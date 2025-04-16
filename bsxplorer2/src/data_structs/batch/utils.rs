use crate::data_structs::batch::traits::{colnames, BsxBatchMethods};
use itertools::Itertools;
use polars::frame::DataFrame;
use polars::prelude::Column;

pub fn merge_replicates<B: BsxBatchMethods>(
    mut batches: Vec<B>,
    count_agg: fn(Vec<&Column>) -> Column,
    density_agg: fn(Vec<&Column>) -> Column,
) -> anyhow::Result<B> {
    if batches.is_empty() {
        anyhow::bail!("batches cannot be empty");
    } else if batches.len() == 1 {
        return Ok(batches.pop().unwrap());
    } else {
        let chr_col = batches[0]
            .data()
            .column(colnames::CHR_NAME)?;
        let pos_col = batches[0]
            .data()
            .column(colnames::CHR_NAME)?;
        let strand_col = batches[0]
            .data()
            .column(colnames::STRAND_NAME)?;
        let context_col = batches[0]
            .data()
            .column(colnames::CONTEXT_NAME)?;

        let count_m_col = count_agg(
            batches
                .iter()
                .map(|b| {
                    b.data()
                        .column(colnames::COUNT_M_NAME)
                        .unwrap()
                })
                .collect_vec(),
        );
        let count_total_col = count_agg(
            batches
                .iter()
                .map(|b| {
                    b.data()
                        .column(colnames::COUNT_TOTAL_NAME)
                        .unwrap()
                })
                .collect_vec(),
        );
        let density_col = density_agg(
            batches
                .iter()
                .map(|b| {
                    b.data()
                        .column(colnames::DENSITY_NAME)
                        .unwrap()
                })
                .collect_vec(),
        );
        let df = DataFrame::from_iter([
            chr_col.to_owned(),
            pos_col.to_owned(),
            strand_col.to_owned(),
            context_col.to_owned(),
            count_m_col,
            count_total_col,
            density_col,
        ]);

        let batch = unsafe { B::new_unchecked(df) };
        Ok(batch)
    }
}
