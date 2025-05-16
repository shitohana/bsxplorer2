use std::io::{Read, Seek};

use anyhow::anyhow;
use polars::error::PolarsResult;
use polars::export::arrow::array::Array;
use polars::export::arrow::record_batch::RecordBatchT;
use polars::frame::DataFrame;

use super::ipc::IpcFileReader;
use super::BatchIndex;
use crate::data_structs::batch::{BsxBatch, BsxBatchBuilder};
use crate::data_structs::coords::Contig;
use crate::data_structs::typedef::BsxSmallStr;

/// Reader for BSX files
pub struct BsxFileReader<R: Read + Seek> {
    ipc_reader: IpcFileReader<R>,
    index:      Option<BatchIndex<BsxSmallStr, u32>>,
}

impl<R: Read + Seek> BsxFileReader<R> {
    /// Creates a new BSX file reader
    pub fn new(handle: R) -> Self {
        Self {
            ipc_reader: IpcFileReader::new(handle, None, None),
            index:      None,
        }
    }

    pub fn from_file_and_index(
        report: R,
        index: &mut R,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            ipc_reader: IpcFileReader::new(report, None, None),
            index:      Some(BatchIndex::from_file(index)?),
        })
    }

    /// Indexes the BSX file and returns the index.
    pub fn index(&mut self) -> anyhow::Result<&BatchIndex<BsxSmallStr, u32>> {
        let initialized = self.index.is_some();
        if initialized {
            Ok(self.index.as_ref().unwrap())
        }
        else {
            let mut new_index = BatchIndex::new();

            for batch_idx in 0..self.ipc_reader.blocks_total() {
                let batch = self
                    .get_batch(batch_idx)
                    .expect("Batch index out of bounds")?;
                // We unwrap as we expect that there is NO empty batches in bsx
                // file
                let contig = batch.as_contig()
                    .ok_or(anyhow!("Batch with no data encountered in file"))?;

                new_index.insert(contig, batch_idx);
            }
            self.index = Some(new_index);

            Ok(self.index.as_ref().unwrap())
        }
    }

    /// Queries the BSX file for a given contig.
    pub fn query(
        &mut self,
        contig: &Contig<BsxSmallStr, u32>,
    ) -> anyhow::Result<Option<BsxBatch>> {
        let batch_indices = self
            .index()?
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
        batches.sort_by_key(|b| b.first_pos());

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
    ) -> PolarsResult<BsxBatch> {
        batch
            .map(|batch| {
                DataFrame::try_from((
                    batch,
                    self.ipc_reader
                        .metadata()
                        .schema
                        .as_ref(),
                ))
                .expect(
                    "Failed to create DataFrame from batch - schema mismatch",
                )
            })
            .map(|df| unsafe { BsxBatch::new_unchecked(df) })
    }

    /// Retrieves a specific batch by index
    pub fn get_batch(
        &mut self,
        batch_idx: usize,
    ) -> Option<PolarsResult<BsxBatch>> {
        self.ipc_reader
            .read_at(batch_idx)
            .map(|res| self.process_record_batch(res))
    }

    /// Returns the total number of blocks in the file
    pub fn blocks_total(&self) -> usize {
        self.ipc_reader.blocks_total()
    }

    /// Retrieves the next batch
    fn _next(&mut self) -> Option<PolarsResult<BsxBatch>> {
        let next = self.ipc_reader.next();
        let res = next.map(|res| self.process_record_batch(res));

        if let Some(Ok(data)) = res.as_ref() {
            if data.len() == 0 {
                return None;
            }
        }

        res
    }

    /// Creates an iterator for the BSX file.
    pub fn iter(&mut self) -> BsxFileIterator<R> {
        BsxFileIterator { reader: self }
    }

    pub fn set_index(
        &mut self,
        index: Option<BatchIndex<BsxSmallStr, u32>>,
    ) {
        self.index = index;
    }
}

impl<R: Read + Seek> Iterator for BsxFileReader<R> {
    type Item = PolarsResult<BsxBatch>;

    /// Returns the next batch or None when finished
    fn next(&mut self) -> Option<Self::Item> {
        self._next()
    }
}

/// Iterator for BSX files.
pub struct BsxFileIterator<'a, R: Read + Seek> {
    reader: &'a mut BsxFileReader<R>,
}

impl<R: Read + Seek> Iterator for BsxFileIterator<'_, R> {
    type Item = PolarsResult<BsxBatch>;

    /// Returns the next batch or None when finished
    fn next(&mut self) -> Option<Self::Item> {
        self.reader._next()
    }

    fn nth(
        &mut self,
        n: usize,
    ) -> Option<Self::Item> {
        self.reader
            .ipc_reader
            .set_current_block(self.reader.ipc_reader.current_block() + n);
        self.reader._next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (
            self.reader.ipc_reader.blocks_total(),
            Some(self.reader.ipc_reader.blocks_total()),
        )
    }

    fn count(self) -> usize
    where
        Self: Sized, {
        self.reader.ipc_reader.blocks_total()
            - self.reader.ipc_reader.current_block()
    }

    fn last(self) -> Option<Self::Item>
    where
        Self: Sized, {
        self.reader
            .get_batch(self.reader.blocks_total() - 1)
    }
}
