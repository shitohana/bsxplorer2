use std::{collections::VecDeque, marker::PhantomData};

use polars::{frame::DataFrame, io::mmap::MmapBytesReader, prelude::DataType};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::data_structs::batch::{builder::BsxBatchBuilder, decoded::BsxBatch, encoded::EncodedBsxBatch, traits::BsxBatchMethods};

use super::{read::OwnedBatchedCsvReader, schema::ReportTypeSchema};

struct ReportReader<R: MmapBytesReader + 'static, B: BsxBatchMethods> {
    csv_reader: OwnedBatchedCsvReader<R>,
    report_schema: ReportTypeSchema,
    batch_per_read: usize,
    cache: VecDeque<B>,
    chr_dtype: Option<DataType>
}

impl<R: MmapBytesReader + 'static> ReportReader<R, BsxBatch> {
    fn process_batches(&self, batches: Vec<DataFrame>) -> VecDeque<BsxBatch> {
        batches.into_par_iter()
            .map(
                |b| {
                    b.partition_by_stable([self.report_schema.chr_col()], true)
                        .expect("Couldn't partition batch by chromosome")
                })
            .flatten()
            .map(|b| {
                BsxBatchBuilder::default()
                    .with_report_type(self.report_schema)
                    .with_check_single_chr(false)
                    .build_decoded(b)
                    .expect("Failed to build batch")
            })
            .collect()
    }
}

impl<R: MmapBytesReader + 'static> Iterator for ReportReader<R, BsxBatch> {
    type Item = BsxBatch;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(batch) = self.cache.pop_front() {
            Some(batch)
        } else {
            match self.csv_reader.next_batches(self.batch_per_read) {
                Ok(Some(new_batches)) => {
                    let mut processed = self.process_batches(new_batches);
                    self.cache.append(&mut processed);
                    self.next()
                },
                Ok(None) => {None},
                Err(e) => panic!("{e}")
            }
        }
    }
}

impl<R: MmapBytesReader + 'static> ReportReader<R, EncodedBsxBatch> {
    fn process_batches(&self, batches: Vec<DataFrame>) -> VecDeque<EncodedBsxBatch> {
        let chr_dtype = self.chr_dtype.clone()
            .expect("chr dtype must be specified to read into EncodedBsxBatch");

        batches.into_par_iter()
            .map(
                |b| {
                    b.partition_by_stable([self.report_schema.chr_col()], true)
                        .expect("Couldn't partition batch by chromosome")
                })
            .flatten()
            .map(|b| {
                BsxBatchBuilder::default()
                    .with_report_type(self.report_schema)
                    .with_check_single_chr(false)
                    .build_encoded(b, chr_dtype.clone())
                    .expect("Failed to build batch")
            })
            .collect()
    }
}

impl<R: MmapBytesReader + 'static> Iterator for ReportReader<R, EncodedBsxBatch> {
    type Item = EncodedBsxBatch;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(batch) = self.cache.pop_front() {
            Some(batch)
        } else {
            match self.csv_reader.next_batches(self.batch_per_read) {
                Ok(Some(new_batches)) => {
                    let mut processed = self.process_batches(new_batches);
                    self.cache.append(&mut processed);
                    self.next()
                },
                Ok(None) => {None},
                Err(e) => panic!("{e}")
            }
        }
    }
}

struct BsxBatchCache<B: BsxBatchMethods> {
    cache: Option<B>
}
