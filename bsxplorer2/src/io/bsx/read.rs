use std::collections::HashMap;
use std::io::{ Read, Seek };
use std::ops::Range;

use bio::data_structures::interval_tree::IntervalTree;
use itertools::Itertools;
use polars::error::PolarsResult;
use polars::export::arrow::array::Array;
use polars::export::arrow::record_batch::RecordBatchT;
use polars::frame::DataFrame;

use super::ipc::IpcFileReader;
use crate::data_structs::batch::{ BsxBatchMethods, LazyBsxBatch };
use crate::data_structs::batch::EncodedBsxBatch;
use crate::data_structs::coords::{ Contig, GenomicPosition };

type BatchIndexMap = HashMap<String, IntervalTree<GenomicPosition<String, u32>, usize>>;

/// Reader for BSX files
pub struct BsxFileReader<R: Read + Seek> {
    reader: IpcFileReader<R>,
    index: Option<BatchIndexMap>,
}

impl<R: Read + Seek> BsxFileReader<R> {
    /// Creates a new BSX file reader
    pub fn new(handle: R) -> Self {
        Self {
            reader: IpcFileReader::new(handle, None, None),
            index: None,
        }
    }

    pub fn index(&mut self) -> anyhow::Result<&BatchIndexMap> {
        let initialized = self.index.is_some();
        if initialized {
            return Ok(self.index.as_ref().unwrap());
        } else {
            let mut new_index = BatchIndexMap::new();

            for batch_idx in 0..self.reader.blocks_total() {
                let batch = self.get_batch(batch_idx).unwrap()?;
                let contig = batch.as_contig()?;

                new_index
                    .entry(contig.seqname().to_string())
                    .and_modify(|tree| {
                        tree.insert(Range::<_>::from(contig.clone()), batch_idx);
                    })
                    .or_insert_with(|| {
                        let mut tree = IntervalTree::new();
                        tree.insert(Range::<_>::from(contig), batch_idx);
                        tree
                    });
            }
            self.index = Some(new_index.clone());

            Ok(self.index.as_ref().unwrap())
        }
    }

    pub fn find(&mut self, contig: Contig<String, u32>) -> Option<Vec<usize>> {
        let index = self.index().expect("Failed to retrieve index");
        if let Some(tree) = index.get(&contig.seqname()) {
            let batches = tree
                .find(Range::<_>::from(contig))
                .map(|entry| entry.data().clone())
                .collect_vec();
            Some(batches)
        } else {
            None
        }
    }

    pub fn query(
        &mut self,
        contig: Contig<String, u32>
    ) -> Option<anyhow::Result<EncodedBsxBatch>> {
        let batches = self.find(contig.clone())?;
        if batches.is_empty() {
            return None;
        }

        let mut results = Vec::new();
        for batch_idx in batches {
            if let Some(batch) = self.get_batch(batch_idx) {
                results.push(batch);
            }
        }
        let batches = match results.into_iter().collect::<PolarsResult<Vec<_>>>() {
            Ok(mut batches) => {
                let mut first = batches.remove(0);
                for batch in batches.drain(..) {
                    first.extend(&batch);
                }
                Ok(first)
            }
            Err(e) => Err(anyhow::anyhow!(e)),
        };

        match batches {
            Ok(batch) =>
                Some(
                    LazyBsxBatch::from(batch)
                        .filter_pos_gt(contig.start().position().saturating_sub(1))
                        .filter_pos_lt(contig.end().position())
                        .collect()
                ),
            Err(e) => Some(Err(e)),
        }
    }

    /// Processes a record batch into an EncodedBsxBatch
    fn process_record_batch(
        &self,
        batch: PolarsResult<RecordBatchT<Box<dyn Array>>>
    ) -> PolarsResult<EncodedBsxBatch> {
        batch
            .map(|batch| {
                DataFrame::try_from((batch, self.reader.metadata().schema.as_ref())).expect(
                    "Failed to create DataFrame from batch - schema mismatch"
                )
            })
            .map(|df| unsafe { EncodedBsxBatch::new_unchecked(df) })
    }

    /// Retrieves a specific batch by index
    pub fn get_batch(&mut self, batch_idx: usize) -> Option<PolarsResult<EncodedBsxBatch>> {
        self.reader.read_at(batch_idx).map(|res| self.process_record_batch(res))
    }

    /// Returns the total number of blocks in the file
    pub fn blocks_total(&self) -> usize {
        let count = self.reader.blocks_total();
        count
    }
}

impl<R: Read + Seek> Iterator for BsxFileReader<R> {
    type Item = PolarsResult<EncodedBsxBatch>;

    /// Returns the next batch or None when finished
    fn next(&mut self) -> Option<Self::Item> {
        let next = self.reader.next();
        let res = next.map(|res| self.process_record_batch(res));

        if let Some(Ok(data)) = res.as_ref() {
            if data.height() == 0 {
                return None;
            }
        }

        res
    }
}
