use itertools::Itertools;
use polars::error::PolarsResult;

use super::BsxFileReader;
use crate::data_structs::batch::{
    merge_replicates,
    AggMethod,
    BsxBatch,
};

pub struct MultiBsxFileReader {
    readers: Vec<BsxFileReader>,
}

impl MultiBsxFileReader {
    pub fn n_readers(&self) -> usize {
        self.readers.len()
    }

    pub fn validate(
        &mut self,
        deep: bool,
    ) -> anyhow::Result<()> {
        if self
            .readers
            .iter()
            .map(BsxFileReader::blocks_total)
            .unique()
            .count()
            > 1
        {
            anyhow::bail!("BSX files have different block counts")
        }
        let mut iterators = self.readers.iter_mut().map(|r| r.iter()).collect_vec();
        let mut batch_count = 0;

        while let Some(batches) = iterators
            .iter_mut()
            .map(|i| i.next())
            .collect::<Option<PolarsResult<Vec<_>>>>()
            .transpose()?
        {
            if !batches.iter().map(|b| b.as_contig()).all_equal() {
                anyhow::bail!(
                    "BSX files have different contig lengths (batch {})",
                    batch_count
                )
            }

            if deep
                && !batches
                    .iter()
                    .map(|b| b.position().to_vec_null_aware().unwrap_left())
                    .all_equal()
            {
                anyhow::bail!(
                    "BSX files have different positions (batch {})",
                    batch_count
                )
            }

            batch_count += 1;
        }

        Ok(())
    }

    pub fn get_batch(
        &mut self,
        batch_idx: usize,
    ) -> Option<Vec<PolarsResult<BsxBatch>>> {
        self.readers
            .iter_mut()
            .map(|r| r.get_batch(batch_idx))
            .collect::<Option<Vec<_>>>()
    }

    pub fn get_batch_merged(
        &mut self,
        batch_idx: usize,
        count_agg: AggMethod,
        density_agg: AggMethod,
    ) -> Option<PolarsResult<BsxBatch>> {
        let batches = match self
            .get_batch(batch_idx)?
            .into_iter()
            .collect::<PolarsResult<Vec<_>>>()
        {
            Ok(batches) => batches,
            Err(e) => return Some(Err(e)),
        };
        Some(merge_replicates(
            batches,
            count_agg.get_expr(),
            density_agg.get_expr(),
        ))
    }

    pub fn blocks_total(&self) -> usize {
        self.readers.first().unwrap().blocks_total()
    }

    pub fn iter(&mut self) -> MultiBsxIterator<'_> {
        MultiBsxIterator {
            inner: self
                .readers
                .iter_mut()
                .map(|r| {
                    Box::new(r.iter())
                        as Box<dyn Iterator<Item = PolarsResult<BsxBatch>> + '_>
                })
                .collect_vec(),
        }
    }

    pub fn iter_merged(
        &mut self,
        count_agg: AggMethod,
        density_agg: AggMethod,
    ) -> impl Iterator<Item = PolarsResult<BsxBatch>> + '_ {
        self.iter().map(move |batches| -> PolarsResult<BsxBatch> {
            merge_replicates(batches?, count_agg.get_expr(), density_agg.get_expr())
        })
    }
}

impl FromIterator<BsxFileReader> for MultiBsxFileReader {
    fn from_iter<I: IntoIterator<Item = BsxFileReader>>(iter: I) -> Self {
        let readers = iter.into_iter().collect_vec();
        assert!(!readers.is_empty(), "No readers provided");
        MultiBsxFileReader { readers }
    }
}

pub struct MultiBsxIterator<'a> {
    inner: Vec<Box<dyn Iterator<Item = PolarsResult<BsxBatch>> + 'a>>,
}

impl Iterator for MultiBsxIterator<'_> {
    type Item = PolarsResult<Vec<BsxBatch>>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.iter_mut().map(|r| r.next()).collect()
    }
}
