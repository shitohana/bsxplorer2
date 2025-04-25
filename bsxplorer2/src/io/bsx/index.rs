use core::panic;
use bio_types::annot::pos;
use itertools::Itertools;
use memmap2::Mmap;
use num::{ToPrimitive, Unsigned};
use polars::prelude::*;
use polars_arrow::array::TryExtend;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Debug;
use std::io::{BufRead, Cursor, Read, Seek};
use std::marker::PhantomData;
use std::ops::Bound::{Excluded, Included};
use std::os::fd::AsRawFd;
use anyhow::anyhow;
use super::ipc::IpcFileReader;
use crate::data_structs::batch::{
    BatchType, BsxBatchBuilder, BsxBatchMethods, BsxTypeTag,
    EncodedBsxBatch, LazyBsxBatch,
};
use crate::io::report::{BedGraphRow, BismarkRow, CgMapRow, CoverageRow, ReportRow, ReportTypeSchema};
use crate::utils::get_categorical_dtype;

trait GenomicQuery {
    fn query<C: AsRef<str>, N: Unsigned>(
        &self,
        chr: C,
        start: N,
        end: N,
    ) -> impl BsxBatchMethods + BsxTypeTag;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
struct ChunkLabel {
    chr_idx: usize,
    start: usize,
    end: usize,
}

impl PartialOrd for ChunkLabel {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<std::cmp::Ordering> {
        if self.chr_idx == other.chr_idx {
            if self.end >= other.end {
                Some(std::cmp::Ordering::Greater)
            } else if self.start <= other.start {
                Some(std::cmp::Ordering::Less)
            } else {
                Some(std::cmp::Ordering::Equal)
            }
        } else {
            Some(self.chr_idx.cmp(&other.chr_idx))
        }
    }
}

impl Ord for ChunkLabel {
    fn cmp(
        &self,
        other: &Self,
    ) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl ChunkLabel {
    fn new(
        chr_idx: usize,
        start: usize,
        end: usize,
    ) -> Self {
        Self {
            chr_idx,
            start,
            end,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ChrIndex(HashMap<String, usize>);

impl ChrIndex {
    fn new() -> Self {
        Self(HashMap::new())
    }

    fn index(
        &mut self,
        chr: String,
    ) -> usize {
        let l = self.0.len();
        self.0
            .entry(chr)
            .or_insert_with(|| l)
            .clone()
    }

    fn keys(&self) -> Vec<String> {
        self.0.keys().cloned().collect()
    }
}

trait IndexedHandle<'a> {
    type Value: Clone + Debug + Serialize + Deserialize<'a>;
    type Output: BsxBatchMethods;

    fn get_tree(&self) -> &BTreeMap<ChunkLabel, Self::Value>;
    fn get_tree_mut(&mut self) -> &mut BTreeMap<ChunkLabel, Self::Value>;

    fn get_chr_index(&self) -> &ChrIndex;
    fn get_chr_index_mut(&mut self) -> &mut ChrIndex;

    fn range<C: AsRef<str>, N: Unsigned + ToPrimitive + Copy>(
        &mut self,
        chr: C,
        start: N,
        end: N,
    ) -> anyhow::Result<Option<Self::Output>>;

    fn serialize_tree<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> Result<S::Ok, S::Error> {
        serde::Serialize::serialize(self.get_tree(), serializer)
    }

    fn serialize_chr_index<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> Result<S::Ok, S::Error> {
        serde::Serialize::serialize(self.get_chr_index(), serializer)
    }

    fn insert_record<C: AsRef<str>>(
        &mut self,
        chr: C,
        start: usize,
        end: usize,
        value: Self::Value,
    ) {
        let chr_idx = self
            .get_chr_index_mut()
            .index(chr.as_ref().to_string());
        let label = ChunkLabel::new(chr_idx, start, end);
        self.get_tree_mut().insert(label, value);
    }

    fn remove_record<C: AsRef<str>>(
        &mut self,
        chr: C,
        start: usize,
        end: usize,
    ) -> Option<Self::Value> {
        let chr_idx = self
            .get_chr_index_mut()
            .index(chr.as_ref().to_string());
        let label = ChunkLabel::new(chr_idx, start, end);
        self.get_tree_mut().remove(&label)
    }

    /// Query the chr:[start; end) region.
    fn range_records<C: AsRef<str>, N: Unsigned + ToPrimitive>(
        &mut self,
        chr: C,
        start: N,
        end: N,
    ) -> Vec<Self::Value> {
        let chr_idx = self
            .get_chr_index_mut()
            .index(chr.as_ref().to_string());
        let start = start.to_usize().unwrap();
        let end = end.to_usize().unwrap();

        let range = (
            Included(ChunkLabel::new(chr_idx, start, start)),
            Excluded(ChunkLabel::new(chr_idx, end, end)),
        );

        let result = self
            .get_tree()
            .range(range)
            .map(|(k, v)| v.clone())
            .collect::<Vec<_>>();
        result
    }
}

struct IndexedBsxReader<'a, R: Read + Seek> {
    reader: IpcFileReader<R>,
    chr_index: ChrIndex,
    batch_index: BTreeMap<ChunkLabel, usize>,
    _lifetime: &'a (),
}

impl<'a, R: Read + Seek> IndexedBsxReader<'a, R> {
    pub fn try_new(handle: R) -> anyhow::Result<Self> {
        let mut reader = IpcFileReader::new(handle, None, None);
        let mut batch_index = BTreeMap::new();
        let mut chr_index = ChrIndex::new();
        for idx in 0..reader.blocks_total() {
            let data = reader.read_df_at(idx)?;
            let batch = BsxBatchBuilder::no_checks().build::<EncodedBsxBatch>(data)?;

            let chr = batch.chr_val()?;
            let chr_idx = chr_index.index(chr.to_string());
            let start = batch.start_pos().ok_or(anyhow!("no data"))?;
            let end = batch.end_pos().ok_or(anyhow!("no data"))?;

            let label = ChunkLabel::new(chr_idx, start as usize, end as usize);
            batch_index.insert(label, idx);
        }
        reader = reader.reopen();
        Ok(Self {
            reader,
            chr_index,
            batch_index,
             _lifetime: &(),
        })
    }
}

impl<'a, R: Read + Seek> IndexedBsxReader<'a, R> {}

impl<'a, R: Read + Seek> IndexedHandle<'a> for IndexedBsxReader<'a, R> {
    type Value = usize;
    type Output = EncodedBsxBatch;

    fn get_tree(&self) -> &BTreeMap<ChunkLabel, Self::Value> {
        &self.batch_index
    }

    fn get_tree_mut(&mut self) -> &mut BTreeMap<ChunkLabel, Self::Value> {
        &mut self.batch_index
    }

    fn get_chr_index(&self) -> &ChrIndex {
        &self.chr_index
    }

    fn get_chr_index_mut(&mut self) -> &mut ChrIndex {
        &mut self.chr_index
    }

    fn range<C: AsRef<str>, N: Unsigned + ToPrimitive + Copy>(
        &mut self,
        chr: C,
        start: N,
        end: N,
    ) -> anyhow::Result<Option<Self::Output>> {
        let batch_indices = self.range_records(chr, start, end);
        if batch_indices.is_empty() {
            return Ok(None);
        }

        let batches = batch_indices
            .into_iter()
            .map(|i| {
                self.reader
                    .read_at(i)
                    .unwrap_or_else(|| {
                        panic!("Batch exists in the index, but not in file")
                    })
            })
            .collect_vec();
        let mut df = DataFrame::empty_with_arrow_schema(
            self.reader.metadata().schema.as_ref(),
        );
        df.try_extend(batches)?;

        let batch = BsxBatchBuilder::no_checks()
            .with_rechunk(true)
            .with_chr_dtype(Some(df.column("chr")?.dtype().clone()))
            .build::<EncodedBsxBatch>(df)?;
        let batch_filtered = LazyBsxBatch::from(batch)
            .filter_pos_gt(
                start
                    .to_u64()
                    .unwrap()
                    .saturating_sub(1),
            )
            .filter_pos_lt(end.to_u64().unwrap())
            .collect()?;

        Ok(Some(batch_filtered))
    }
}

struct IndexedReportReader<'a, B: BsxBatchMethods> {
    mmap: Mmap,
    report_type: ReportTypeSchema,
    chr_index: ChrIndex,
    /// Values are offsets in the file to the start (inclusive)
    /// and end (exclusive) offsets of chunk's data
    batch_index: BTreeMap<ChunkLabel, usize>,
    _lifetime: &'a (),
    _batch: PhantomData<B>,
}

impl<'a, B: BsxBatchMethods> IndexedReportReader<'a, B> {
    pub fn new<F: memmap2::MmapAsRawDesc>(
        handle: F,
        report_type: ReportTypeSchema,
        block_size: usize,
    ) -> Self {
        if report_type.need_align() {
            unimplemented!("Region reader is not yet implemented for report types, which need alignment.")
        }
        let mmap = unsafe { Mmap::map(handle).unwrap() };
        todo!()
    }
}

impl<'a, B: BsxBatchMethods> IndexedHandle<'_> for IndexedReportReader<'a, B> {
    type Value = usize;
    type Output = B;

    fn get_tree(&self) -> &BTreeMap<ChunkLabel, Self::Value> {
        &self.batch_index
    }

    fn get_tree_mut(&mut self) -> &mut BTreeMap<ChunkLabel, Self::Value> {
        &mut self.batch_index
    }

    fn get_chr_index(&self) -> &ChrIndex {
        &self.chr_index
    }

    fn get_chr_index_mut(&mut self) -> &mut ChrIndex {
        &mut self.chr_index
    }

    fn range<C: AsRef<str>, N: Unsigned + ToPrimitive + Copy>(
        &mut self,
        chr: C,
        start: N,
        end: N,
    ) -> anyhow::Result<Option<Self::Output>> {
        let batch_indices = self.range_records(chr.as_ref(), start, end);
        if batch_indices.is_empty() {
            return Ok(None);
        };
        let start_offset = *batch_indices.first().unwrap();

        let cursor = Cursor::new(&self.mmap[start_offset..(start_offset + CHUNK_SIZE)]);
        let mut res = self
            .report_type
            .read_options()
            .into_reader_with_file_handle(cursor)
            .finish()?
            .partition_by_stable([self.report_type.chr_col()], true)?
            .into_iter()
            .map(|df| {
                BsxBatchBuilder::all_checks()
                    .with_check_single_chr(false)
                    .with_chr_dtype({
                        if matches!(B::type_enum(), BatchType::Encoded) {
                            Some(get_categorical_dtype(self.chr_index.keys()))
                        } else {
                            None
                        }
                    })
                    .with_report_type(self.report_type)
                    .build::<B>(df).unwrap()
            })
            .filter(|batch| batch.chr_val().unwrap() == chr.as_ref())
            .collect_vec();

        if res.len() == 0 { anyhow::bail!("Index maps to wrong location. Chromosome mismatch") }
        let batch = res.pop().unwrap();

        let batch_filtered = LazyBsxBatch::from(batch)
            .filter_pos_gt(
                start
                    .to_u64()
                    .unwrap()
                    .saturating_sub(1),
            )
            .filter_pos_lt(end.to_u64().unwrap())
            .collect()?;

        Ok(Some(batch_filtered))
    }
}

/// 64 kb
const CHUNK_SIZE: usize = 64 * 1024;
fn index_buffer<B: AsRef<[u8]>>(buffer: B, report_type: ReportTypeSchema) -> anyhow::Result<BTreeMap<ChunkLabel, usize>> {
    todo!()
}

/// Extracts the label from the buffer based on the report type.
/// 
/// # Arguments
/// * `buffer` - The buffer containing the data (single CSV row).
/// * `report_type` - The type of report (Bismark, BedGraph, CgMap, Coverage).
fn extract_label(
    buffer: &[u8],
    report_type: &ReportTypeSchema,
) -> anyhow::Result<(String, usize)> {
    let record = report_type.header_record();
    let (chr, pos) = match report_type {
        ReportTypeSchema::Bismark => {
            let row: BismarkRow = record.deserialize(None)?;
            (row.get_chr(), row.get_pos())
        },
        ReportTypeSchema::BedGraph => {
            let row: BedGraphRow = record.deserialize(None)?;
            (row.get_chr(), row.get_pos())
        },
        ReportTypeSchema::CgMap => {
            let row: CgMapRow = record.deserialize(None)?;
            (row.get_chr(), row.get_pos())
        },
        ReportTypeSchema::Coverage => {
            let row: CoverageRow = record.deserialize(None)?;
            (row.get_chr(), row.get_pos())
        },
    };
    Ok((chr, pos))
}