use super::report_read_utils::{align_data_with_context, get_context_data};
use crate::data_structs::bsx_batch::BsxBatch;
use crate::data_structs::region::GenomicPosition;
use crate::io::report::fasta_reader::{FastaCoverageReader, FastaReader};
use crate::io::report::schema::ReportTypeSchema;
#[cfg(feature = "python")]
use crate::utils::wrap_box_result;
use crate::utils::{first_position, last_position};
use itertools::Itertools;
use log::debug;
use num::{PrimInt, Unsigned};
use polars::df;
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::io::mmap::MmapBytesReader;
use polars::io::RowIndex;
use polars::prelude::{BatchedCsvReader, CsvReader, Schema, SchemaRef};
#[cfg(feature = "python")]
use pyo3::exceptions::PyIOError;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3_polars::error::PyPolarsErr;
#[cfg(feature = "python")]
use pyo3_polars::PyDataFrame;
use serde::Serialize;
use std::error::Error;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::path::PathBuf;
use std::sync::mpsc::{Receiver, SyncSender};
use std::sync::{mpsc, Arc};
use std::thread;
use std::thread::JoinHandle;
/*
TODO:
    Add decompression
    Add expand_user
*/

pub struct ReportReaderBuilder {
    pub report_type: ReportTypeSchema,
    pub rechunk: bool,
    pub n_threads: Option<usize>,
    pub low_memory: bool,
    pub n_rows: Option<usize>,
    pub row_index: Option<RowIndex>,
    pub chunk_size: usize,
    pub skip_rows_after_header: usize,
    pub fasta_path: Option<PathBuf>,
    pub fai_path: Option<PathBuf>,
    pub batch_per_read: usize,
    pub batch_size: usize,
}

impl ReportReaderBuilder {
    pub fn new(report_type: ReportTypeSchema) -> ReportReaderBuilder {
        ReportReaderBuilder {
            report_type,
            ..Default::default()
        }
    }

    pub fn with_rechunk(mut self, rechunk: bool) -> Self {
        self.rechunk = rechunk;
        self
    }
    pub fn with_fasta(mut self, fasta_path: PathBuf, fai_path: PathBuf) -> Self {
        self.fasta_path = Some(fasta_path);
        self.fai_path = Some(fai_path);
        self
    }
    pub fn with_n_threads(mut self, n_threads: usize) -> Self {
        self.n_threads = Some(n_threads);
        self
    }
    pub fn with_skip_rows_after_header(mut self, skip_rows_after_header: usize) -> Self {
        self.skip_rows_after_header = skip_rows_after_header;
        self
    }
    pub fn with_low_memory(mut self, low_memory: bool) -> Self {
        self.low_memory = low_memory;
        self
    }
    pub fn with_n_rows(mut self, n_rows: usize) -> Self {
        self.n_rows = Some(n_rows);
        self
    }
    pub fn with_report_type(mut self, report_type: ReportTypeSchema) -> Self {
        self.report_type = report_type;
        self
    }
    pub fn with_row_index(mut self, row_index: RowIndex) -> Self {
        self.row_index = Some(row_index);
        self
    }
    pub fn with_chunk_size(mut self, chunk_size: usize) -> Self {
        self.chunk_size = chunk_size;
        self
    }
    pub fn with_batch_per_read(mut self, n: usize) -> Self {
        self.batch_per_read = n;
        self
    }
    pub fn with_batch_size(mut self, n: usize) -> Self {
        self.batch_size = n;
        self
    }

    pub fn try_finish<F>(self, handle: F) -> Result<ReportReader, Box<dyn Error>>
    where
        F: Read + Seek + MmapBytesReader + 'static,
    {
        let csv_reader = self
            .report_type
            .read_options()
            .with_low_memory(self.low_memory)
            .with_n_rows(self.n_rows)
            .with_skip_rows_after_header(self.skip_rows_after_header)
            .with_row_index(self.row_index.clone())
            .with_chunk_size(self.batch_size)
            .with_low_memory(self.low_memory)
            .with_n_threads(self.n_threads)
            .into_reader_with_file_handle(handle);

        let indexed_reader = match (self.fasta_path, self.fai_path) {
            (Some(fasta_path), Some(fai_path)) => {
                let reader = FastaReader::try_from_handle(
                    BufReader::new(File::open(&fasta_path)?),
                    BufReader::new(File::open(&fai_path)?),
                )?;
                Some(FastaCoverageReader::from(reader))
            }
            (Some(_), None) => todo!("Add auto FASTA indexing?"),
            (None, Some(_)) => return Err(Box::from("No FASTA path given")),
            (None, None) => None,
        };

        if self.report_type.need_align() && indexed_reader.is_none() {
            return Err(Box::from(format!(
                "FASTA path must be provided for report type: {:?}",
                self.report_type
            )));
        }

        Ok(ReportReader::new(
            csv_reader,
            indexed_reader,
            self.report_type,
            self.batch_per_read,
            self.chunk_size,
        ))
    }
}

impl Default for ReportReaderBuilder {
    fn default() -> Self {
        ReportReaderBuilder {
            report_type: ReportTypeSchema::Bismark,
            rechunk: false,
            n_threads: None,
            low_memory: false,
            n_rows: None,
            row_index: None,
            chunk_size: 10_000,
            skip_rows_after_header: 0,
            fasta_path: None,
            fai_path: None,
            batch_per_read: 16,
            batch_size: 2 << 20,
        }
    }
}

pub struct OwnedBatchedCsvReader<F>
where
    F: Read + Seek + MmapBytesReader,
{
    #[allow(dead_code)]
    // this exists because we need to keep ownership
    pub schema: SchemaRef,
    pub batched_reader: BatchedCsvReader<'static>,
    // keep ownership
    pub _reader: CsvReader<F>,
}

impl<F> OwnedBatchedCsvReader<F>
where
    F: Read + Seek + MmapBytesReader + 'static,
{
    pub(crate) fn new(mut reader: CsvReader<F>, schema: Arc<Schema>) -> Self {
        let batched_reader = reader
            .batched_borrowed()
            .expect("Could not create batched CSV reader.");
        let batched_reader: BatchedCsvReader<'static> =
            unsafe { std::mem::transmute(batched_reader) };
        Self {
            batched_reader,
            _reader: reader,
            schema,
        }
    }
}

impl<F> OwnedBatchedCsvReader<F>
where
    F: Read + Seek + MmapBytesReader + 'static,
{
    pub fn next_batches(&mut self, n: usize) -> PolarsResult<Option<Vec<DataFrame>>> {
        self.batched_reader.next_batches(n)
    }
}

pub struct ContextData<N>
where
    N: PrimInt + Unsigned + Serialize + Display,
{
    positions: Vec<N>,
    contexts: Vec<Option<bool>>,
    strands: Vec<bool>,
}

impl<N: PrimInt + Unsigned + Serialize + Display> ContextData<N> {
    pub fn len(&self) -> usize {
        self.contexts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.contexts.is_empty()
    }

    pub(crate) fn new() -> Self {
        ContextData {
            positions: Vec::new(),
            contexts: Vec::new(),
            strands: Vec::new(),
        }
    }

    pub(crate) fn filter<F: Fn(N) -> bool>(&self, predicate: F) -> Self {
        let mut new_self = Self::new();
        for (idx, pos) in self.positions.iter().enumerate() {
            if predicate(*pos) {
                new_self.add_row(*pos, self.contexts[idx], self.strands[idx])
            }
        }
        new_self
    }

    #[allow(dead_code)]
    pub(crate) fn col_names() -> &'static [&'static str] {
        &["position", "context", "strand"]
    }

    pub(crate) fn position_col() -> &'static str {
        "position"
    }

    pub(crate) fn with_capacity(capacity: usize) -> Self {
        ContextData {
            positions: Vec::with_capacity(capacity),
            contexts: Vec::with_capacity(capacity),
            strands: Vec::with_capacity(capacity),
        }
    }

    pub(crate) fn add_row(&mut self, position: N, context: Option<bool>, strand: bool) {
        self.positions.push(position);
        self.strands.push(strand);
        self.contexts.push(context);
    }

    pub fn from_sequence(seq: &[u8], start: GenomicPosition<N>) -> Self {
        let start_pos = start.position();
        let fw_bound: usize = seq.len() - 2;
        let rv_bound: usize = 2;

        let mut new = Self::with_capacity(seq.len());

        let ascii_seq = seq.to_ascii_uppercase();

        'seq_iter: for (index, nuc) in ascii_seq.iter().enumerate() {
            let forward = match nuc {
                b'C' => true,
                b'G' => false,
                _ => continue 'seq_iter,
            };
            let context = if forward {
                if index >= fw_bound {
                    continue 'seq_iter;
                };

                if ascii_seq[index + 1] == b'G' {
                    Some(true)
                } else if ascii_seq[index + 2] == b'G' {
                    Some(false)
                } else {
                    None
                }
            } else {
                if index <= rv_bound {
                    continue 'seq_iter;
                };

                if ascii_seq[index - 1] == b'C' {
                    Some(true)
                } else if ascii_seq[index - 2] == b'C' {
                    Some(false)
                } else {
                    None
                }
            };

            new.add_row(start_pos + N::from(index).unwrap(), context, forward);
        }

        new.shrink_to_fit();
        new
    }

    pub(crate) fn shrink_to_fit(&mut self) {
        self.positions.shrink_to_fit();
        self.contexts.shrink_to_fit();
        self.strands.shrink_to_fit();
    }

    pub fn into_dataframe(self) -> PolarsResult<DataFrame> {
        df![
            "position" => self.positions.iter().map(|x| x.to_u64().unwrap()).collect_vec(),
            "context" => self.contexts,
            "strand" => self.strands
        ]
    }

    pub fn positions(&self) -> &Vec<N> {
        &self.positions
    }

    pub fn contexts(&self) -> &Vec<Option<bool>> {
        &self.contexts
    }

    pub fn strands(&self) -> &Vec<bool> {
        &self.strands
    }
}

/// Data itself and marker if the batch is the last
pub(in crate::io::report) type ReadQueueItem = (DataFrame, bool);

#[cfg_attr(feature = "python", pyclass)]
pub struct ReportReader {
    _join_handle: JoinHandle<()>,
    receiver: Receiver<ReadQueueItem>,
    report_schema: ReportTypeSchema,
    fasta_reader: Option<FastaCoverageReader<BufReader<File>, u64>>,
    chunk_size: usize,
    batch_cache: Option<DataFrame>,
}

impl ReportReader {
    pub(crate) fn new<F>(
        reader: CsvReader<F>,
        fasta_reader: Option<FastaCoverageReader<BufReader<File>, u64>>,
        report_schema: ReportTypeSchema,
        batch_per_read: usize,
        chunk_size: usize,
    ) -> Self
    where
        F: MmapBytesReader + 'static,
    {
        let (sender, receiver) = mpsc::sync_channel(batch_per_read);
        let join_handle =
            thread::spawn(move || reader_thread(reader, sender, report_schema, batch_per_read));
        // Create struct
        Self {
            receiver,
            chunk_size,
            report_schema,
            fasta_reader,
            _join_handle: join_handle,
            batch_cache: None,
        }
    }

    /// Item must has single chromosome
    fn extend_cache(&mut self, item: ReadQueueItem) -> Result<(), Box<dyn Error>> {
        let context_data = if let Some(reader) = self.fasta_reader.as_mut() {
            let data = get_context_data(reader, &item, &self.report_schema)?;
            Some(data)
        } else {
            None
        };

        let data_bsx = self.report_schema.bsx_mutate()(item.0)?;

        let result = if let Some(context_data) = context_data {
            align_data_with_context(&data_bsx, context_data)?
        } else {
            data_bsx
        };

        if let Some(cache) = self.batch_cache.as_mut() {
            cache.extend(&result)?;
        } else {
            self.batch_cache = Some(result);
        }
        Ok(())
    }

    /// It is caller responsibility to check, if cache size is enough
    fn get_chunk(&mut self) -> Result<DataFrame, Box<dyn Error>> {
        if let Some(cache) = self.batch_cache.as_mut() {
            let (mut output, mut rest) = cache.split_at(self.chunk_size as i64);
            let chr_col = output.column("chr")?.as_series().unwrap().rechunk();
            let start_chr = chr_col.first().as_any_value().str_value().to_string();
            let end_chr = chr_col
                .get(self.chunk_size - 1)
                .unwrap()
                .str_value()
                .to_string();

            if start_chr != end_chr {
                let mut partitioned = output.partition_by_stable(["chr"], true)?;
                partitioned.reverse();
                output = partitioned.pop().unwrap();

                if let Some(new_rest) = partitioned
                    .into_iter()
                    .reduce(|acc, df| acc.vstack(&df).unwrap())
                {
                    rest = new_rest.vstack(&rest)?;
                }
            }

            self.batch_cache = Some(rest);
            Ok(output)
        } else {
            Err(Box::from("No batch cache found"))
        }
    }
}

impl Iterator for ReportReader {
    type Item = BsxBatch;

    fn next(&mut self) -> Option<Self::Item> {
        if self
            .batch_cache
            .as_ref()
            .map(DataFrame::height)
            .unwrap_or(0)
            >= self.chunk_size
        {
            let df = self.get_chunk().expect("Could not get chunk");
            Some(unsafe { BsxBatch::new_unchecked(df) })
        } else {
            match self.receiver.recv() {
                Ok(data) => {
                    debug!("Received data with {} rows", data.0.height());
                    self.extend_cache(data)
                        .expect("Failed extending batch cache");
                    self.next()
                }
                Err(_) => {
                    if let Some(cache) = self.batch_cache.take() {
                        debug!("Batch cache is emptied");
                        let res = unsafe { BsxBatch::new_unchecked(cache.clone()) };
                        Some(res)
                    } else {
                        None
                    }
                }
            }
        }
    }
}

#[cfg(feature = "python")]
#[pymethods]
impl ReportReader {
    #[new]
    #[pyo3(signature = (
        report_path,
        report_type,
        /,
        rechunk=false,
        n_threads=None,
        low_memory=false,
        n_rows=None,
        chunk_size=10000,
        skip=0,
        fasta_path=None,
        fai_path=None,
        batch_per_read=16,
        batch_size=2 << 20
    ))]
    pub fn new_py(
        report_path: String,
        report_type: ReportTypeSchema,
        rechunk: bool,
        n_threads: Option<usize>,
        low_memory: bool,
        n_rows: Option<usize>,
        chunk_size: usize,
        skip: usize,
        fasta_path: Option<PathBuf>,
        fai_path: Option<PathBuf>,
        batch_per_read: usize,
        batch_size: usize,
    ) -> PyResult<Self> {
        let builder = ReportReaderBuilder {
            report_type,
            rechunk,
            n_threads,
            low_memory,
            n_rows,
            row_index: None,
            chunk_size,
            skip_rows_after_header: skip,
            fasta_path,
            fai_path,
            batch_per_read,
            batch_size,
        };
        let file = BufReader::new(File::open(report_path)?);
        wrap_box_result!(PyIOError, builder.try_finish(file))
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<BsxBatch> {
        slf.next()
    }
}

fn reader_thread<R>(
    reader: CsvReader<R>,
    send_channel: SyncSender<ReadQueueItem>,
    report_schema: ReportTypeSchema,
    batch_per_read: usize,
) where
    R: MmapBytesReader + 'static,
{
    let mut owned_batched = OwnedBatchedCsvReader::new(reader, Arc::from(report_schema.schema()));
    let chr_col = report_schema.chr_col();
    let pos_col = report_schema.position_col();

    let mut cached_batch: Option<DataFrame> = None;

    loop {
        let incoming_batches = owned_batched
            .next_batches(batch_per_read)
            .expect("Error while reading batches");
        if let Some(batches_vec) = incoming_batches {
            // Flatten batches and partition by chr
            let batches_vec = batches_vec
                .into_iter()
                .flat_map(|batch| partition_batch(batch, chr_col, pos_col).unwrap())
                .collect_vec();

            for batch in batches_vec.into_iter() {
                if let Some(cached_batch_data) = cached_batch.take() {
                    // Compare with cached
                    let cached_chr = first_position(&cached_batch_data, chr_col, pos_col)
                        .unwrap()
                        .chr()
                        .to_owned();
                    let new_chr = first_position(&batch, chr_col, pos_col)
                        .unwrap()
                        .chr()
                        .to_owned();
                    let item = if cached_chr != new_chr {
                        (cached_batch_data, true)
                    } else {
                        (cached_batch_data, false)
                    };
                    // Return cached
                    send_channel
                        .send(item)
                        .expect("Could not send data to main thread.");
                }
                // Update cached
                cached_batch = Some(batch);
            }
        } else {
            // Release cached
            if let Some(cached_batch_data) = cached_batch.take() {
                let item = (cached_batch_data, true);
                send_channel
                    .send(item)
                    .expect("Could not send data to main thread.");
            }
            break;
        }
    }
}

/// It is caller responsibility to validate schema
fn partition_batch(
    batch: DataFrame,
    chr_col: &str,
    pos_col: &str,
) -> Result<Vec<DataFrame>, Box<dyn Error>> {
    let first_pos = first_position(&batch, chr_col, pos_col)?;
    let last_pos = last_position(&batch, chr_col, pos_col)?;

    if first_pos.chr() != last_pos.chr() {
        Ok(batch.partition_by_stable([chr_col], true)?)
    } else {
        Ok(vec![batch])
    }
}
