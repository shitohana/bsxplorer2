use std::collections::{BTreeMap, HashMap};
use std::io::{Read, Seek};

use log::{debug, trace};
use polars::error::PolarsResult;
use polars::export::arrow::array::Array;
use polars::export::arrow::record_batch::RecordBatchT;
use polars::frame::DataFrame;

use crate::data_structs::batch::BsxBatchMethods;
use crate::data_structs::batch::EncodedBsxBatch;
use super::ipc::IpcFileReader;

/// Maps chromosome names to positions and batch indices
/// Structure: chromosome -> (start_position -> batch_index)
pub type BSXIndex = HashMap<String, BTreeMap<u64, usize>>;

/// Reader for BSX files
pub struct BsxFileReader<R: Read + Seek> {
    reader: IpcFileReader<R>,
}

impl<R: Read + Seek> BsxFileReader<R> {
    /// Creates a new BSX file reader
    pub fn new(handle: R) -> Self {
        Self {
            reader: IpcFileReader::new(handle, None, None),
        }
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
                    self.reader.metadata().schema.as_ref(),
                ))
                .expect(
                    "Failed to create DataFrame from batch - schema mismatch",
                )
            })
            .map(|df| {
                unsafe { EncodedBsxBatch::new_unchecked(df) }
            })
    }

    /// Retrieves a specific batch by index
    pub fn get_batch(
        &mut self,
        batch_idx: usize,
    ) -> Option<PolarsResult<EncodedBsxBatch>> {
        self.reader
            .read_at(batch_idx)
            .map(|res| self.process_record_batch(res))
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
