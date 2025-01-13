use std::error::Error;
use itertools::Itertools;
use polars::prelude::*;
use crate::io::report::reader::ReportReader;

pub mod schema;
pub mod bsx_batch;


pub mod reader {
    use super::*;
    use crate::io::report::bsx_batch::BsxBatch;
    use crate::io::report::report_batch_utils::{align_data_with_context, first_position, get_context_data, last_position};
    use crate::region::{GenomicPosition, RegionCoordinates};
    use bio::io::fasta::{IndexedReader, Sequence};
    use polars::io::mmap::MmapBytesReader;
    use polars::io::RowIndex;
    use std::error::Error;
    use std::fs::File;
    use std::io::{Read, Seek};
    use std::ops::{Add, Sub};
    use std::path::PathBuf;
    use std::sync::mpsc;
    use std::sync::mpsc::Receiver;
    use std::thread::JoinHandle;
    use std::{io, thread};
    use std::borrow::Cow;
    use log::{debug, info};
    use crate::io::report::schema;
    use crate::io::report::schema::ReportTypeSchema;

    pub struct ReportReaderBuilder {
        report_type: schema::ReportTypeSchema,
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
        sequence_buffer_size: Option<usize>,
    }
    
    pub fn new(report_type: schema::ReportTypeSchema) -> ReportReaderBuilder {
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
        pub fn with_report_type(mut self, report_type: schema::ReportTypeSchema) -> Self {
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
        pub fn with_sequence_buffer_size(mut self, size: Option<usize>) -> Self {
            self.sequence_buffer_size = size;
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
                (Some(fasta_path), Some(fai_path)) => Some(FastaCoverageReader::new(
                    fasta_path,
                    fai_path,
                    self.sequence_buffer_size.unwrap_or(1_000_000),
                )),
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
                sequence_buffer_size: None,
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

    pub struct FastaCoverageReader {
        reader: IndexedReader<File>,
        buffer: Vec<u8>,
        seq_buffer_size: usize,
        metadata: Vec<(String, Sequence)>,
        last_position: GenomicPosition,
    }

    impl FastaCoverageReader {
        pub fn new(fa_path: PathBuf, fai_path: PathBuf, sequence_buffer_size: usize) -> Self {
            let reader = IndexedReader::new(
                File::open(&fa_path)
                    .unwrap_or_else(|_| panic!("Could not open FASTA file {}", fa_path.display())),
                File::open(&fai_path).unwrap_or_else(|_| {
                    panic!("Could not open FASTA index file {}", fai_path.display())
                }),
            )
            .expect("Could not create Indexed reader");

            Self::from_indexed_reader(reader, sequence_buffer_size)
        }

        pub(crate) fn from_indexed_reader(
            reader: IndexedReader<File>,
            sequence_buffer_size: usize,
        ) -> Self {
            let seq_metadata = Self::get_seq_metadata(&reader);
            assert!(!seq_metadata.is_empty(), "Index is empty!");

            let (first_chr, _) = seq_metadata.first().unwrap();
            let sequence_buffer = Vec::with_capacity(sequence_buffer_size);
            let last_position = GenomicPosition::new(first_chr.to_string(), 0);

            Self {
                reader,
                buffer: sequence_buffer,
                seq_buffer_size: sequence_buffer_size,
                metadata: seq_metadata,
                last_position,
            }
        }

        pub fn buffer(&self) -> &Vec<u8> {
            &self.buffer
        }

        pub fn seq_buffer_size(&self) -> usize {
            self.seq_buffer_size
        }

        pub fn metadata(&self) -> &Vec<(String, Sequence)> {
            &self.metadata
        }

        pub fn last_position(&self) -> &GenomicPosition {
            &self.last_position
        }

        pub fn get_chr_metadata(&self, chr: &str) -> Option<&Sequence> {
            self.metadata
                .iter()
                .find(|(name, _)| name == chr)
                .map(|(_, sequence)| sequence)
        }

        pub fn get_seq_region(&mut self, region: &RegionCoordinates) -> io::Result<Vec<u8>> {
            let mut sequence = Vec::new();

            self.reader
                .fetch(region.chr(), region.start.into(), region.end.into())?;
            self.reader.read(&mut sequence)?;
            Ok(sequence)
        }

        pub fn get_seq_until(&mut self, pos: &GenomicPosition) -> io::Result<Vec<u8>> {
            match self.check_chr(pos.chr()) {
                Some(false) => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Last chromosome reading not finished",
                    ))
                }
                Some(true) => {
                    self.last_position = GenomicPosition::new(pos.chr().parse().unwrap(), 0);
                }
                _ => {}
            }

            let region = (self.tell_pos() >> pos.clone()).unwrap();
            let (start_i, end_i) = loop {
                if let Some((start, stop)) = self.find_buffer_slice(&region) {
                    break (start, stop);
                } else {
                    let new_end = self.get_shifted_end();
                    self.fill_buffer(new_end)?;
                }
            };
            debug!("Fetching region {}", region);
            
            if let Some(sequence) = self.buffer.get(start_i..end_i) {
                let res = Ok(sequence.to_vec());
                self.drain_buffer(end_i);
                res
            } else {
                Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Error reading sequence",
                ))
            }
        }

        pub fn get_seq_leftover(&mut self) -> io::Result<Vec<u8>> {
            let last_position = self.get_chr_metadata(self.last_position.chr()).unwrap().len;
            self.get_seq_until(&GenomicPosition::new(
                self.last_position.chr().to_string(),
                last_position.try_into().unwrap(),
            ))
        }

        pub fn tell(&self) -> usize {
            <usize as Sub>::sub(
                self.last_position.position().try_into().unwrap(),
                self.buffer.len(),
            ) + 1
        }

        pub fn tell_pos(&self) -> GenomicPosition {
            GenomicPosition::new(
                self.last_position.chr().to_string(),
                self.tell().try_into().unwrap(),
            )
        }
        
        pub fn chr_finished(&self) -> bool {
            self.get_chr_metadata(self.last_position.chr()).unwrap().len == self.last_position().position() as u64
        }

        fn get_seq_metadata(reader: &IndexedReader<File>) -> Vec<(String, Sequence)> {
            reader
                .index
                .sequences()
                .iter()
                .map(|seq| (seq.name.clone(), seq.clone()))
                .collect_vec()
        }

        fn get_shifted_end(&mut self) -> GenomicPosition {
            let mut res_position = <u64 as Add>::add(
                self.seq_buffer_size.try_into().unwrap(),
                self.last_position.position().into(),
            );
            let current_chr = self.last_position.chr();
            let (_, chr_data) = self
                .metadata
                .iter()
                .find(|(name, _)| name == current_chr)
                .unwrap();

            if chr_data.len.lt(&res_position) {
                res_position = chr_data.len;
            }
            GenomicPosition::new(current_chr.to_string(), res_position.try_into().unwrap())
        }

        pub(crate) fn set_chr(&mut self, chr: &str) -> Result<(), Box<dyn Error>> {
            if self.get_chr_metadata(chr).is_some() {
                self.last_position = GenomicPosition::new(chr.to_string(), 0);
                self.buffer.clear();
                Ok(())
            } else {
                Err(Box::from(format!("Chr {} not found in index", chr)))
            }
        }
        
        /// Does not check position outof chromosome bounds
        pub(crate) fn set_position(&mut self, position: u64) {
            let new_last_position = GenomicPosition::new(
                self.last_position.chr().to_string(),
                position.try_into().unwrap(),
            );
            self.last_position = new_last_position;
        }

        fn find_buffer_slice(&self, region: &RegionCoordinates) -> Option<(usize, usize)> {
            // Check same chromosome
            assert_eq!(region.chr(), self.last_position.chr());

            let last_pos_usize: usize = self.last_position.position().try_into().unwrap();
            let buffer_start = self.tell();
            assert!(
                buffer_start.le(&region.start().try_into().unwrap()),
                "Sequence has already been read till {}, but requested position is {}",
                buffer_start,
                region.start()
            );

            if last_pos_usize.gt(&region.end().try_into().unwrap()) {
                let region_start_idx =
                    <usize as Sub>::sub(region.start().try_into().unwrap(), buffer_start);
                let region_end_idx =
                    <usize as Sub>::sub(region.end().try_into().unwrap(), buffer_start);
                Some((region_start_idx, region_end_idx))
            } else {
                None
            }
        }

        fn drain_buffer(&mut self, until: usize) {
            let (_, new_buf) = self.buffer.split_at(until);
            self.buffer = Vec::from(new_buf)
        }

        fn fill_buffer(&mut self, new_pos: GenomicPosition) -> io::Result<()> {
            // Check same chromosome
            assert_eq!(new_pos.chr(), self.last_position.chr());

            let read_region = (self.last_position.clone() >> new_pos).unwrap();

            let temp_buf = self.get_seq_region(&read_region)?;

            self.buffer.extend(temp_buf);
            self.last_position = read_region.end_gpos();
            Ok(())
        }

        fn check_chr(&self, chr: &str) -> Option<bool> {
            if self.last_position.chr() == chr {
                None
            } else if self.last_position.position()
                == self
                    .metadata
                    .iter()
                    .find(|(name, _)| name == chr)
                    .unwrap()
                    .1
                    .len as u32
            {
                Some(true)
            } else {
                Some(false)
            }
        }
    }

    pub struct ContextData {
        positions: Vec<u32>,
        contexts: Vec<Option<bool>>,
        strands: Vec<bool>,
    }

    impl ContextData {
        pub fn len(&self) -> usize {
            self.contexts.len()
        }
        
        pub(crate) fn new() -> Self {
            ContextData {
                positions: Vec::new(),
                contexts: Vec::new(),
                strands: Vec::new(),
            }
        }
        
        pub(crate) fn col_names() -> &'static[&'static str] {
            &[
                "position",
                "context",
                "strand",
            ]
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

        pub(crate) fn add_row(&mut self, position: u32, context: Option<bool>, strand: bool) {
            self.positions.push(position);
            self.strands.push(strand);
            self.contexts.push(context);
        }

        pub fn from_sequence(seq: &[u8], start: GenomicPosition) -> Self {
            let start_pos = start.position() - 1;
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

                new.add_row((start_pos + index as u32), context, forward);
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
                "position" => self.positions,
                "context" => self.contexts,
                "strand" => self.strands
            ]
        }

        pub fn positions(&self) -> &Vec<u32> {
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
        join_handle: JoinHandle<()>,
        receiver: Receiver<ReadQueueItem>,
        report_schema: schema::ReportTypeSchema,
        fasta_reader: Option<FastaCoverageReader>,
        chunk_size: usize,
        batch_cache: DataFrame,
    }

    impl ReportReader {
        pub(crate) fn new<F>(
            reader: CsvReader<F>,
            fasta_reader: Option<FastaCoverageReader>,
            report_schema: schema::ReportTypeSchema,
            batch_per_read: usize,
            chunk_size: usize,
        ) -> Self
        where
            F: MmapBytesReader + 'static,
        {
            let (join_handle, receiver) =
                Self::init_reader_thread(reader, &report_schema, batch_per_read);
            info!("Reader thread initialized");
            let batch_cache = DataFrame::empty_with_schema(&BsxBatch::schema());
            // Create struct
            Self {
                receiver,
                chunk_size,
                report_schema,
                fasta_reader,
                join_handle,
                batch_cache,
            }
        }

        fn init_reader_thread<F>(
            reader: CsvReader<F>,
            report_schema: &ReportTypeSchema,
            batch_per_read: usize,
        ) -> (JoinHandle<()>, Receiver<ReadQueueItem>)
        where
            F: MmapBytesReader + 'static,
        {
            let (sc, receiver) = mpsc::sync_channel(batch_per_read);

            // Clone schema for thread
            let report_schema_clone = report_schema.clone();
            // Start reading thread
            let join_handle = thread::spawn(move || {
                let mut owned_batched =
                    OwnedBatchedCsvReader::new(reader, Arc::from(report_schema_clone.schema()));
                
                loop {
                    let incoming_batches = owned_batched.next_batches(batch_per_read);
                    match incoming_batches {
                        Ok(Some(data)) => {
                            data.into_iter().for_each(|frame| {
                                let first_pos = first_position(
                                    &frame,
                                    report_schema_clone.chr_col(),
                                    report_schema_clone.position_col(),
                                )
                                    .expect("failed to read first position");
                                let last_pos = last_position(
                                    &frame,
                                    report_schema_clone.chr_col(),
                                    report_schema_clone.position_col(),
                                )
                                    .expect("failed to read last position");

                                if (first_pos >> last_pos).is_some() {
                                    sc.send((frame, false))
                                        .expect("Could not send data to main thread.")
                                } else {
                                    let mut partitioned = frame
                                        .partition_by_stable([report_schema_clone.chr_col()], true)
                                        .expect("Could not get partition by chr");
                                    while let Some(batch) = partitioned.pop() {
                                        if partitioned.is_empty() {
                                            sc.send((batch, false))
                                                .expect("Could not send data to main thread.")
                                        } else {
                                            sc.send((batch, true))
                                                .expect("Could not send data to main thread.")
                                        }
                                    }
                                }
                            });
                        },
                        Ok(None) => {break}
                        Err(e) => {panic!("{}", e)}
                    }
                }
            });

            (join_handle, receiver)
        }
        
        /// Item must has single chromosome
        fn extend_cache(&mut self, item: ReadQueueItem) -> Result<(), Box<dyn Error>> {
            let context_data = if let Some(reader) = self.fasta_reader.as_mut() {
                Some(get_context_data(reader, &item, &self.report_schema)?)
            } else {
                None
            };
            
            let data_bsx = (self.report_schema.bsx_mutate())(item.0)?;
            
            let result = if let Some(context_data) = context_data {
                align_data_with_context(&data_bsx, context_data)?
            } else { 
                data_bsx 
            };
            
            self.batch_cache.extend(&result)?;
            Ok(())
        }
    }

    impl Iterator for ReportReader {
        type Item = BsxBatch;

        fn next(&mut self) -> Option<Self::Item> {
            if self.batch_cache.height() >= self.chunk_size {
                let chr_col = self.batch_cache
                    .column("chr").unwrap()
                    .as_series()?.rechunk();
                
                let start_chr = chr_col.first()
                    .as_any_value().str_value()
                    .to_string();
                let end_chr = chr_col
                    .get(self.chunk_size - 1).unwrap()
                    .str_value()
                    .to_string();
                
                let split_at = if start_chr == end_chr {
                    self.chunk_size as i64
                } else {
                    chr_col.iter()
                        .position(|v| v.str_value() != start_chr.clone())
                        .unwrap() as i64
                };

                let (output, rest) = self.batch_cache.split_at(split_at);
                self.batch_cache = rest;
                Some(BsxBatch::new(output))
                
            } else {
                match self.receiver.recv() {
                    Ok(data) => { 
                        debug!("Received data with {} rows", data.0.height());
                        self.extend_cache(data)
                            .expect("Failed extending batch cache");
                        self.next()
                    }
                    Err(_) => {
                        if self.batch_cache.is_empty() { None } else { 
                            let res = BsxBatch::new(self.batch_cache.clone());
                            self.batch_cache = DataFrame::empty();
                            Some(res)
                        }
                    }
                }
            }
        }
    }
    
}



pub mod report_batch_utils {
    use crate::io::report::bsx_batch::BsxBatch;
    use crate::io::report::reader::{ContextData, FastaCoverageReader, ReadQueueItem};
    use crate::io::report::schema::ReportTypeSchema;
    use super::*;
    use crate::region::GenomicPosition;
    
    pub(crate) fn get_context_data(
        reader: &mut FastaCoverageReader, 
        item: &ReadQueueItem,
        report_schema: &ReportTypeSchema,
    )  -> PolarsResult<ContextData> {
        let (data, is_end) = item;
        let last_position = last_position(
            data, 
            report_schema.chr_col(), 
            report_schema.position_col()
        )? + 1;

        if reader.chr_finished() {
            reader.set_chr(last_position.chr()).
                expect("Could not set chr in fasta reader.");
        }

        let sequence = if *is_end {
            // Read sequence till the end of the chromosome
            let seq = reader.get_seq_leftover()?;
            // Set position explicitly to mark chromosome as read
            let chr_length = reader.get_chr_metadata(reader.last_position().chr()).unwrap().len;
            reader.set_position(chr_length);
            seq
        } else {
            // Read region
            reader.get_seq_until(&last_position)?
        };

        let first_position = last_position - sequence.len() as u32;
        let context_data = ContextData::from_sequence(&sequence, first_position);
        Ok(context_data)
    }
    
    const JOIN_ARGS: JoinArgs = JoinArgs {
        how: JoinType::Left,
        validation: JoinValidation::OneToOne,
        suffix: None,
        slice: None,
        join_nulls: false,
        coalesce: JoinCoalesce::CoalesceColumns,
        maintain_order: MaintainOrderJoin::Left
    };
    
    pub(crate) fn align_data_with_context(
        data_frame: &DataFrame, 
        context_data: ContextData,
    ) -> PolarsResult<DataFrame> {
        let data_join_columns = [BsxBatch::pos_col()];
        let context_join_columns = [ContextData::position_col()];
        let chr = first_position(data_frame, BsxBatch::chr_col(), BsxBatch::pos_col())?.chr().to_string();
        
        let context_df = context_data.into_dataframe()?;
        let mut context_df_lazy = context_df.lazy()
            .cast(PlHashMap::from_iter(
                data_join_columns.iter().cloned().map(|name| {(name, data_frame.schema().get(name).unwrap().clone())})
            ), true)
            .with_column(lit(chr).alias("chr"));
        let drop_columns = context_df_lazy.collect_schema()?
            .iter_names()
            .filter(|name| !data_join_columns.contains(&name.as_str()) && data_frame.column(name).is_ok())
            .cloned().collect::<Vec<_>>();
        
        context_df_lazy = decode_context(context_df_lazy, "context", "context");
        context_df_lazy = decode_strand(context_df_lazy, "strand", "strand");

        context_df_lazy.collect()?
            .join(
                &data_frame.drop_many(drop_columns),
                context_join_columns,
                data_join_columns,
                JOIN_ARGS
            )?.lazy()
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
    ) -> PolarsResult<GenomicPosition> {
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
        Ok(GenomicPosition::new(chr, pos as u32))
    }

    pub(crate) fn last_position(
        data: &DataFrame,
        chr_col: &str,
        pos_col: &str,
    ) -> PolarsResult<GenomicPosition> {
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
            .last()
            .unwrap();
        Ok(GenomicPosition::new(chr, pos as u32))
    }
    
    pub fn encode_strand(lazy_frame: LazyFrame, strand_col: &str) -> LazyFrame {
        lazy_frame.with_column(
            when(col(strand_col).eq(lit("+")))
                .then(lit(true))
                .when(col(strand_col).eq(lit("-")))
                .then(lit(false))
                .otherwise(lit(NULL))
                .cast(DataType::Boolean),
        )
    }

    pub fn encode_context(lazy_frame: LazyFrame, context_col: &str) -> LazyFrame {
        lazy_frame.with_column(
            when(col(context_col).eq(lit("CG")))
                .then(lit(true))
                .when(col(context_col).eq(lit("CHG")))
                .then(lit(false))
                .otherwise(lit(NULL))
                .cast(DataType::Boolean),
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

    pub fn decode_context(lazy_frame: LazyFrame, context_col: &str, result_name: &str) -> LazyFrame {
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
