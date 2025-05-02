use anyhow::anyhow;
use bio::data_structures::interval_tree::IntervalTree;
use indexmap::IndexSet;
use itertools::Itertools;
use num::{PrimInt, Unsigned};
use polars::error::PolarsResult;
use polars::export::arrow::array::Array;
use polars::export::arrow::record_batch::RecordBatchT;
use polars::frame::DataFrame;
use std::collections::HashMap;
use std::hash::Hash;
use std::io::{Read, Seek};
use std::ops::Range;

use super::ipc::IpcFileReader;
use crate::data_structs::batch::BsxBatchMethods;
use crate::data_structs::batch::{BsxBatchBuilder, EncodedBsxBatch};
use crate::data_structs::coords::{Contig, GenomicPosition};

/// Index for batches in a BSX file.
#[derive(Debug, Clone)]
pub struct BatchIndex<S, P>
where
    S: AsRef<str> + Clone,
    P: Unsigned + PrimInt,
{
    map: HashMap<S, IntervalTree<GenomicPosition<S, P>, usize>>,
    chr_order: IndexSet<S>,
}

impl<S, P> BatchIndex<S, P>
where
    S: AsRef<str> + Clone + Eq + Hash,
    P: Unsigned + PrimInt,
{
    fn new() -> Self {
        Self {
            map: HashMap::new(),
            chr_order: IndexSet::new(),
        }
    }

    /// Insert a contig and its corresponding batch index.
    pub fn insert(&mut self, contig: Contig<S, P>, batch_idx: usize) {
        self.chr_order.insert(contig.seqname().clone());

        self.map
            .entry(contig.seqname())
            .and_modify(|tree| {
                tree.insert(
                    Range::<_>::from(contig.clone()),
                    batch_idx,
                );
            })
            .or_insert_with(|| {
                let mut tree = IntervalTree::new();
                tree.insert(Range::<_>::from(contig), batch_idx);
                tree
            });
    }

    /// Sort a set of contigs according to the chromosome order and start position.
    pub fn sort<I>(&self, contigs: I) -> Vec<Contig<S, P>>
    where
        I: IntoIterator<Item = Contig<S, P>>,
    {
        contigs.into_iter()
            .map(|contig| (self.chr_order.get_index_of(&contig.seqname()).unwrap_or(0), contig))
            .sorted_by(|(left_chr, left_contig), (right_chr, right_contig)| left_chr.cmp(&right_chr).then(left_contig.start().cmp(&right_contig.start())))
            .map(|(_, contig)| contig)
            .collect::<Vec<_>>()
    }

    /// Find the batch indices that overlap with a given contig.
    pub fn find(
        &self,
        contig: &Contig<S, P>,
    ) -> Option<Vec<usize>> {
        if let Some(tree) = self.map.get(&contig.seqname()) {
            let batches = tree
                .find(Range::<_>::from(contig.clone()))
                .map(|entry| entry.data().clone())
                .collect_vec();
            Some(batches)
        } else {
            None
        }
    }

    /// Returns the chromosome order.
    pub fn chr_order(&self) -> &IndexSet<S> {
        &self.chr_order
    }

    /// Returns the underlying map.
    pub fn map(&self) -> &HashMap<S, IntervalTree<GenomicPosition<S, P>, usize>> {
        &self.map
    }
}

/// Reader for BSX files
pub struct BsxFileReader<R: Read + Seek> {
    ipc_reader: IpcFileReader<R>,
    index: Option<BatchIndex<String, u32>>,
}

impl<R: Read + Seek> BsxFileReader<R> {
    /// Creates a new BSX file reader
    pub fn new(handle: R) -> Self {
        Self {
            ipc_reader: IpcFileReader::new(handle, None, None),
            index: None,
        }
    }

    /// Indexes the BSX file and returns the index.
    pub fn index(&mut self) -> anyhow::Result<&BatchIndex<String, u32>> {
        let initialized = self.index.is_some();
        if initialized {
            Ok(self.index.as_ref().unwrap())
        } else {
            let mut new_index = BatchIndex::new();

            for batch_idx in 0..self.ipc_reader.blocks_total() {
                let batch = self
                    .get_batch(batch_idx)
                    .expect("Batch index out of bounds")?;
                // We unwrap as we expect that there is NO empty batches in bsx file
                let contig = batch.as_contig()?.unwrap();

                new_index.insert(contig, batch_idx);
            }
            self.index = Some(new_index);

            Ok(self.index.as_ref().unwrap())
        }
    }

    /// Queries the BSX file for a given contig.
    pub fn query(
        &mut self,
        contig: &Contig<String, u32>,
    ) -> anyhow::Result<Option<EncodedBsxBatch>> {
        let batch_indices = self.index()?
            .find(contig)
            .ok_or(anyhow!("No batches found for contig: {}", contig))?;
        if batch_indices.is_empty() {
            return Ok(None);
        }

        let mut batches = batch_indices
            .into_iter()
            .map(|idx| {
                self.get_batch(idx)
                    .ok_or(anyhow!("Batch with index {} not found", idx))
            })
            .collect::<anyhow::Result<PolarsResult<Vec<_>>>>()??;
        batches.sort_by_key(|b| {
            b.start_pos()
                .expect("Unexpected no data")
        });

        let res = BsxBatchBuilder::concat(batches)?
            .lazy()
            .filter_pos_gt(contig.start() - 1)
            .filter_pos_lt(contig.end())
            .collect()?;
        Ok(Some(res))
    }

    /// Processes a record batch into an EncodedBsxBatch
    fn process_record_batch(
        &self,
        batch: PolarsResult<RecordBatchT<Box<dyn Array>>>,
    ) -> PolarsResult<EncodedBsxBatch> {
        batch
            .map(|batch| {
                DataFrame::try_from((
                    batch,
                    self.ipc_reader.metadata().schema.as_ref(),
                ))
                .expect(
                    "Failed to create DataFrame from batch - schema mismatch",
                )
            })
            .map(|df| unsafe { EncodedBsxBatch::new_unchecked(df) })
    }

    /// Retrieves a specific batch by index
    pub fn get_batch(
        &mut self,
        batch_idx: usize,
    ) -> Option<PolarsResult<EncodedBsxBatch>> {
        self.ipc_reader
            .read_at(batch_idx)
            .map(|res| self.process_record_batch(res))
    }

    /// Returns the total number of blocks in the file
    pub fn blocks_total(&self) -> usize {
        let count = self.ipc_reader.blocks_total();
        count
    }

    /// Retrieves the next batch
    fn _next(&mut self) -> Option<PolarsResult<EncodedBsxBatch>> {
        let next = self.ipc_reader.next();
        let res = next.map(|res| self.process_record_batch(res));

        if let Some(Ok(data)) = res.as_ref() {
            if data.height() == 0 {
                return None;
            }
        }

        res
    }

    /// Creates an iterator for the BSX file.
    fn iter(&mut self) -> BsxFileIterator<R> {
        BsxFileIterator {
            reader: self,
        }
    }
}

impl<R: Read + Seek> Iterator for BsxFileReader<R> {
    type Item = PolarsResult<EncodedBsxBatch>;

    /// Returns the next batch or None when finished
    fn next(&mut self) -> Option<Self::Item> {
        self._next()
    }
}

/// Iterator for BSX files.
pub struct BsxFileIterator<'a, R: Read + Seek> {
    reader: &'a mut BsxFileReader<R>,
}

impl<'a, R: Read + Seek> Iterator for BsxFileIterator<'a, R> {
    type Item = PolarsResult<EncodedBsxBatch>;

    /// Returns the next batch or None when finished
    fn next(&mut self) -> Option<Self::Item> {
        self.reader._next()
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.reader.ipc_reader.set_current_block(self.reader.ipc_reader.current_block() + n);
        self.reader._next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.reader.ipc_reader.blocks_total(), Some(self.reader.ipc_reader.blocks_total()))
    }

    fn count(self) -> usize
    where
        Self: Sized,
    {
        self.reader.ipc_reader.blocks_total() - self.reader.ipc_reader.current_block()
    }

    fn last(self) -> Option<Self::Item>
    where
        Self: Sized,
    {
        self.reader.get_batch(self.reader.blocks_total() - 1)
    }
}
