use itertools::Itertools;
use polars::prelude::*;
use std::error::Error;

pub mod bsx_batch;
pub mod fasta_reader;
pub mod schema;

pub mod reader {
    use super::*;
    use crate::io::report::bsx_batch::BsxBatch;
    use crate::io::report::report_batch_utils::{
        align_data_with_context, first_position, get_context_data, last_position,
    };
    use crate::io::report::schema::ReportTypeSchema;
    use crate::region::GenomicPosition;
    use log::debug;
    use polars::io::mmap::MmapBytesReader;
    use polars::io::RowIndex;

    use crate::io::report::fasta_reader::{FastaCoverageReader, FastaReader};
    use std::error::Error;
    use std::fs::File;
    use std::io::{BufReader, Read, Seek};
    use std::path::PathBuf;
    use std::sync::mpsc;
    use std::sync::mpsc::{Receiver, SyncSender};
    use std::thread;
    use std::thread::JoinHandle;
    use polars::export::num::{PrimInt, Unsigned};

    pub struct ReportReaderBuilder {
        report_type: ReportTypeSchema,
        rechunk: bool,
        n_threads: Option<usize>,
        low_memory: bool,
        n_rows: Option<usize>,
        row_index: Option<RowIndex>,
        chunk_size: usize,
        skip_rows_after_header: usize,
        fasta_path: Option<PathBuf>,
        fai_path: Option<PathBuf>,
        batch_per_read: usize,
        batch_size: usize,
    }

    pub fn new(report_type: ReportTypeSchema) -> ReportReaderBuilder {
        ReportReaderBuilder {
            report_type,
            ..Default::default()
        }
    }

    impl ReportReaderBuilder {
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
                (Some(_), None) => todo!(),
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
        N: PrimInt + Unsigned {
        positions: Vec<N>,
        contexts: Vec<Option<bool>>,
        strands: Vec<bool>,
    }

    impl<N: num::PrimInt + num::Unsigned> ContextData<N> {
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
    pub(super) type ReadQueueItem = (DataFrame, bool);

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
        /// TODO maybe change it
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

    fn reader_thread<R>(
        reader: CsvReader<R>,
        send_channel: SyncSender<ReadQueueItem>,
        report_schema: ReportTypeSchema,
        batch_per_read: usize,
    ) where
        R: MmapBytesReader + 'static,
    {
        let mut owned_batched =
            OwnedBatchedCsvReader::new(reader, Arc::from(report_schema.schema()));
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
}

pub mod report_batch_utils {
    use super::*;
    use crate::io::report::bsx_batch::BsxBatch;
    use crate::io::report::reader::{ContextData, ReadQueueItem};
    use crate::io::report::schema::ReportTypeSchema;
    use crate::region::{GenomicPosition, RegionCoordinates};
    use std::io::{BufRead, Seek};
    use num::Unsigned;
    use polars::export::num::PrimInt;

    pub(crate) fn get_context_data<R>(
        reader: &mut fasta_reader::FastaCoverageReader<R, u64>,
        item: &ReadQueueItem,
        report_schema: &ReportTypeSchema,
    ) -> Result<ContextData<u64>, Box<dyn Error>>
    where
        R: BufRead + Seek,
    {
        let (data, is_end) = item;
        let last_position =
            last_position(data, report_schema.chr_col(), report_schema.position_col())?;
        let chr = last_position.clone().chr().to_string();
        let chr_coverage = *reader
            .coverage()
            .get(chr.to_string())
            .expect("Chromosome not found in index");

        // To ensure we capture all contexts
        let sequence_overhead = 3;
        let fetch_start = if chr_coverage.read() == 0 || chr_coverage.read() < sequence_overhead {
            0
        } else {
            chr_coverage.read() - sequence_overhead
        };
        let fetch_end = if last_position.position() as u64 + sequence_overhead
            > chr_coverage.total()
            || *is_end
        {
            chr_coverage.total()
        } else {
            last_position.position() as u64 + sequence_overhead
        };
        let fetch_region =
            RegionCoordinates::new(chr.to_string(), fetch_start, fetch_end);

        let sequence = reader.inner_mut().fetch_region(fetch_region.clone())?;
        let context_data = ContextData::from_sequence(&sequence, fetch_region.start_gpos() + 1)
            .filter(|pos| pos > chr_coverage.read())
            .filter(|pos| {
                if !*is_end {
                    pos <= last_position.position()
                } else {
                    true
                }
            });
        if !*is_end {
            reader
                .coverage_mut()
                .shift_to(fetch_region.chr(), last_position.position())?;
        } else {
            reader
                .coverage_mut()
                .shift_to(fetch_region.chr(), chr_coverage.total())?;
        }

        Ok(context_data)
    }

    const JOIN_ARGS: JoinArgs = JoinArgs {
        how: JoinType::Left,
        validation: JoinValidation::OneToOne,
        suffix: None,
        slice: None,
        join_nulls: false,
        coalesce: JoinCoalesce::CoalesceColumns,
        maintain_order: MaintainOrderJoin::Left,
    };

    pub(crate) fn align_data_with_context<N: PrimInt + Unsigned>(
        data_frame: &DataFrame,
        context_data: ContextData<N>,
    ) -> PolarsResult<DataFrame> {
        let data_join_columns = [BsxBatch::pos_col()];
        let context_join_columns = [ContextData::<N>::position_col()];
        let chr = first_position(data_frame, BsxBatch::chr_col(), BsxBatch::pos_col())?
            .chr()
            .to_string();

        let context_df = context_data.into_dataframe()?;
        let mut context_df_lazy = context_df
            .lazy()
            .cast(
                PlHashMap::from_iter(
                    data_join_columns
                        .iter()
                        .cloned()
                        .map(|name| (name, data_frame.schema().get(name).unwrap().clone())),
                ),
                true,
            )
            .with_column(lit(chr).alias("chr"));
        let drop_columns = context_df_lazy
            .collect_schema()?
            .iter_names()
            .filter(|name| {
                !data_join_columns.contains(&name.as_str()) && data_frame.column(name).is_ok()
            })
            .cloned()
            .collect::<Vec<_>>();

        context_df_lazy = decode_context(context_df_lazy, "context", "context");
        context_df_lazy = decode_strand(context_df_lazy, "strand", "strand");

        context_df_lazy
            .collect()?
            .join(
                &data_frame.drop_many(drop_columns),
                context_join_columns,
                data_join_columns,
                JOIN_ARGS,
            )?
            .lazy()
            .cast(BsxBatch::hashmap(), true)
            .collect()?
            .select(BsxBatch::col_names().iter().cloned())
    }

    pub(crate) fn schema_from_arrays(names: &[&str], dtypes: &[DataType]) -> Schema {
        Schema::from_iter(names.iter().cloned().map_into().zip(dtypes.iter().cloned()))
    }

    pub(crate) fn hashmap_from_arrays<'a>(
        names: &[&'a str],
        dtypes: &[DataType],
    ) -> PlHashMap<&'a str, DataType> {
        PlHashMap::from_iter(names.iter().cloned().map_into().zip(dtypes.iter().cloned()))
    }

    pub(crate) fn first_position(
        data: &DataFrame,
        chr_col: &str,
        pos_col: &str,
    ) -> PolarsResult<GenomicPosition<u64>> {
        let chr = data
            .column(chr_col)?
            .as_series()
            .unwrap()
            .first()
            .as_any_value()
            .cast(&DataType::String)
            .str_value()
            .to_string();
        let pos = data
            .column(pos_col)?
            .cast(&DataType::UInt64)?
            .u64()?
            .first()
            .unwrap();
        Ok(GenomicPosition::new(chr, pos))
    }

    pub(crate) fn last_position(
        data: &DataFrame,
        chr_col: &str,
        pos_col: &str,
    ) -> PolarsResult<GenomicPosition<u64>> {
        let chr = data
            .column(chr_col)?
            .as_series()
            .unwrap()
            .last()
            .as_any_value()
            .cast(&DataType::String)
            .str_value()
            .to_string();
        let pos = data
            .column(pos_col)?
            .cast(&DataType::UInt64)?
            .u64()?
            .last()
            .unwrap();
        Ok(GenomicPosition::new(chr, pos))
    }

    pub fn encode_strand(lazy_frame: LazyFrame, strand_col: &str) -> LazyFrame {
        lazy_frame.with_column(
            when(col(strand_col).eq(lit("+")))
                .then(lit(true))
                .when(col(strand_col).eq(lit("-")))
                .then(lit(false))
                .otherwise(lit(NULL))
                .cast(DataType::Boolean)
                .alias("strand"),
        )
    }

    pub fn encode_context(lazy_frame: LazyFrame, context_col: &str) -> LazyFrame {
        lazy_frame.with_column(
            when(col(context_col).eq(lit("CG")))
                .then(lit(true))
                .when(col(context_col).eq(lit("CHG")))
                .then(lit(false))
                .otherwise(lit(NULL))
                .cast(DataType::Boolean)
                .alias(context_col),
        )
    }

    pub fn decode_strand(lazy_frame: LazyFrame, strand_col: &str, result_name: &str) -> LazyFrame {
        lazy_frame.with_column(
            when(col(strand_col).eq(lit(true)))
                .then(lit("+"))
                .when(col(strand_col).eq(lit(false)))
                .then(lit("-"))
                .otherwise(lit("."))
                .cast(DataType::String)
                .alias(result_name),
        )
    }

    pub fn decode_context(
        lazy_frame: LazyFrame,
        context_col: &str,
        result_name: &str,
    ) -> LazyFrame {
        lazy_frame.with_column(
            when(col(context_col).eq(lit(true)))
                .then(lit("CG"))
                .when(col(context_col).eq(lit(false)))
                .then(lit("CHG"))
                .otherwise(lit("CHH"))
                .cast(DataType::String)
                .alias(result_name),
        )
    }
}
