use crate::data_structs::batch::encoded::EncodedBsxBatch;
use crate::data_structs::batch::traits::{BsxBatchMethods, BsxColNames, BsxTypeTag};
use itertools::Itertools;
use polars::frame::DataFrame;
use polars::prelude::{ChunkedArray, Column, IntoSeries};

pub fn merge_replicates<B: BsxBatchMethods + BsxTypeTag>(
    mut batches: Vec<B>,
    count_agg: fn(
        Vec<&ChunkedArray<B::CountType>>,
    ) -> ChunkedArray<B::CountType>,
    density_agg: fn(
        Vec<&ChunkedArray<B::DensityType>>,
    ) -> ChunkedArray<B::DensityType>,
) -> anyhow::Result<B>
where
    ChunkedArray<B::ChrType>: IntoSeries,
    ChunkedArray<B::PosType>: IntoSeries,
    ChunkedArray<B::StrandType>: IntoSeries,
    ChunkedArray<B::ContextType>: IntoSeries,
    ChunkedArray<B::CountType>: IntoSeries,
    ChunkedArray<B::DensityType>: IntoSeries,
{
    if batches.is_empty() {
        anyhow::bail!("batches cannot be empty");
    } else if batches.len() == 1 {
        return Ok(batches.pop().unwrap());
    } else {
        let chr_col = batches[0].chr();
        let pos_col = batches[0].position();
        let strand_col = batches[0].strand();
        let context_col = batches[0].context();

        let count_m_col = count_agg(
            batches
                .iter()
                .map(|b| b.count_m())
                .collect_vec(),
        );
        let count_total_col = count_agg(
            batches
                .iter()
                .map(|b| b.count_total())
                .collect_vec(),
        );
        let density_col = density_agg(
            batches
                .iter()
                .map(|b| b.density())
                .collect_vec(),
        );
        let df = DataFrame::from_iter([
            Column::new(EncodedBsxBatch::CHR_NAME.into(), chr_col.to_owned()),
            Column::new(EncodedBsxBatch::POS_NAME.into(), pos_col.to_owned()),
            Column::new(EncodedBsxBatch::STRAND_NAME.into(), strand_col.to_owned()),
            Column::new(EncodedBsxBatch::CONTEXT_NAME.into(), context_col.to_owned()),
            Column::new(EncodedBsxBatch::COUNT_M_NAME.into(), count_m_col),
            Column::new(EncodedBsxBatch::COUNT_TOTAL_NAME.into(), count_total_col),
            Column::new(EncodedBsxBatch::DENSITY_NAME.into(), density_col),
        ]);

        let batch = unsafe { B::new_unchecked(df) };
        Ok(batch)
    }
}
