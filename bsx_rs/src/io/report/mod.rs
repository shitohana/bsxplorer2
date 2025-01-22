use itertools::Itertools;
use polars::prelude::*;
use std::error::Error;

pub mod fasta_reader;
pub mod schema;
pub mod read;

mod report_read_utils {
    use super::*;
    use crate::bsx_batch::BsxBatch;
    use crate::io::report::schema::ReportTypeSchema;
    use crate::region::RegionCoordinates;
    use std::io::{BufRead, Seek};
    use num::Unsigned;
    use polars::export::num::PrimInt;
    use crate::io::report::read::{ContextData, ReadQueueItem};
    use crate::utils;

    pub(crate) fn get_context_data<R>(
        reader: &mut fasta_reader::FastaCoverageReader<R, u64>,
        item: &ReadQueueItem,
        report_schema: &ReportTypeSchema,
    ) -> Result<ContextData<u64>, Box<dyn Error>>
    where
        R: BufRead + Seek,
    {
        let (data, is_end) = item;
        let last_position =
            utils::last_position(data, report_schema.chr_col(), report_schema.position_col())?;
        let chr = last_position.clone().chr().to_string();
        let chr_coverage = *reader
            .coverage()
            .get(chr.to_string())
            .expect("Chromosome not found in index");

        // To ensure we capture all contexts
        let sequence_overhead = 3;
        let fetch_start = if chr_coverage.read() == 0 || chr_coverage.read() < sequence_overhead {
            0
        } else {
            chr_coverage.read() - sequence_overhead
        };
        let fetch_end = if last_position.position() as u64 + sequence_overhead
            > chr_coverage.total()
            || *is_end
        {
            chr_coverage.total()
        } else {
            last_position.position() as u64 + sequence_overhead
        };
        let fetch_region =
            RegionCoordinates::new(chr.to_string(), fetch_start, fetch_end);

        let sequence = reader.inner_mut().fetch_region(fetch_region.clone())?;
        let context_data = ContextData::from_sequence(&sequence, fetch_region.start_gpos() + 1)
            .filter(|pos| pos > chr_coverage.read())
            .filter(|pos| {
                if !*is_end {
                    pos <= last_position.position()
                } else {
                    true
                }
            });
        if !*is_end {
            reader
                .coverage_mut()
                .shift_to(fetch_region.chr(), last_position.position())?;
        } else {
            reader
                .coverage_mut()
                .shift_to(fetch_region.chr(), chr_coverage.total())?;
        }

        Ok(context_data)
    }

    const JOIN_ARGS: JoinArgs = JoinArgs {
        how: JoinType::Left,
        validation: JoinValidation::OneToOne,
        suffix: None,
        slice: None,
        join_nulls: false,
        coalesce: JoinCoalesce::CoalesceColumns,
        maintain_order: MaintainOrderJoin::Left,
    };

    pub(crate) fn align_data_with_context<N: PrimInt + Unsigned>(
        data_frame: &DataFrame,
        context_data: ContextData<N>,
    ) -> PolarsResult<DataFrame> {
        let data_join_columns = [BsxBatch::pos_col()];
        let context_join_columns = [ContextData::<N>::position_col()];
        let chr = utils::first_position(data_frame, BsxBatch::chr_col(), BsxBatch::pos_col())?
            .chr()
            .to_string();

        let context_df = context_data.into_dataframe()?;
        let mut context_df_lazy = context_df
            .lazy()
            .cast(
                PlHashMap::from_iter(
                    data_join_columns
                        .iter()
                        .cloned()
                        .map(|name| (name, data_frame.schema().get(name).unwrap().clone())),
                ),
                true,
            )
            .with_column(lit(chr).alias("chr"));
        let drop_columns = context_df_lazy
            .collect_schema()?
            .iter_names()
            .filter(|name| {
                !data_join_columns.contains(&name.as_str()) && data_frame.column(name).is_ok()
            })
            .cloned()
            .collect::<Vec<_>>();

        context_df_lazy = utils::decode_context(context_df_lazy, "context", "context");
        context_df_lazy = utils::decode_strand(context_df_lazy, "strand", "strand");

        context_df_lazy
            .collect()?
            .join(
                &data_frame.drop_many(drop_columns),
                context_join_columns,
                data_join_columns,
                JOIN_ARGS,
            )?
            .lazy()
            .cast(BsxBatch::hashmap(), true)
            .collect()?
            .select(BsxBatch::col_names().iter().cloned())
    }
}
