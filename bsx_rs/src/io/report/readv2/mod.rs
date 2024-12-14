use crate::region::RegionCoordinates;
use crate::ubatch2::{report_type_to_bsx, BSXBatch, BSXCol, UBatch};
use crate::utils::types::BSXResult;
use bio::io::fasta::IndexedReader;
use log::*;
use polars::frame::DataFrame;
use polars::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::error::Error;
use std::fs::File;
use std::mem::ManuallyDrop;
use std::path::{Path, PathBuf};
use std::sync::mpmc::{Receiver, Sender};
use std::sync::{mpmc, Mutex};
use std::thread;
use std::thread::{sleep, spawn, JoinHandle};
use itertools::Itertools;
use polars::io::mmap::MmapBytesReader;
use polars_arrow::array::Utf8ViewArray;
use batch_inner::WorkBatch;
use batch_res::ReadFinalBatch;
use crate::io::report::types::ReportType;

// TODO:    1. Fix UUID
//          2. Make batches concat with one another if length < chunk size
//          4. Make better logging



mod batch_res;
mod batch_inner;

#[derive(Debug, PartialOrd, Ord, Eq, PartialEq, Clone, Copy)]
enum BatchStatus {
    /// Initial raw data, no checks have been provided
    InitialRaw = 0,
    /// Data have been checked to have obligatory columns
    ColumnsChecked = 1,
    /// Data have been split to chromosome unique batches
    PartitionedChr = 2,
    /// Data have been aligned with the reference
    Aligned = 3,
    /// Chromosome column have been casted to categorical datatype
    ChrCasted = 4,
    /// All checks have been passed. Data has been converted to [WorkBatch]
    FullChecked = 5,
}

/// Union of [DataFrame] and [ReadFinalBatch] to be able to go into the
/// same queue
union ReadQueueData {
    df: ManuallyDrop<DataFrame>,
    bsx: ManuallyDrop<ReadFinalBatch>,
}

impl ReadQueueData {
    fn new_df(df: DataFrame) -> ReadQueueData {
        Self {
            df: ManuallyDrop::new(df),
        }
    }

    fn new_bsx(bsx: ReadFinalBatch) -> ReadQueueData {
        Self {
            bsx: ManuallyDrop::new(bsx),
        }
    }
}

#[derive(Clone)]
struct ParserThread {
    fasta_reader: Arc<Mutex<IndexedReader<File>>>,
    res_heap: Arc<Mutex<BinaryHeap<ReadFinalBatch>>>,
    cache_heap: Arc<Mutex<BinaryHeap<WorkBatch>>>,
    receiver: Arc<Receiver<WorkBatch>>,
    sender: Arc<Sender<WorkBatch>>,
    chr_type: DataType,
}

impl ParserThread {
    fn new(
        fasta_reader: IndexedReader<File>,
        receiver: Arc<Receiver<WorkBatch>>,
        sender: Arc<Sender<WorkBatch>>,
        chr_type: DataType,
    ) -> Self {
        Self {
            fasta_reader: Arc::new(Mutex::new(fasta_reader)),
            res_heap: Arc::new(Mutex::new(BinaryHeap::new())),
            cache_heap: Arc::new(Mutex::new(BinaryHeap::new())),
            receiver,
            sender,
            chr_type,
        }
    }

    fn run(&mut self) -> JoinHandle<()> {
        let mut cloned_self = self.clone();
        spawn(move || {
            info!("ParserThread started");
            loop {
                if !cloned_self.receiver.is_empty() {
                    let data = cloned_self.receiver.recv().unwrap();
                    cloned_self.cache_heap.lock().unwrap().push(data);
                }
                if cloned_self.cache_heap.lock().unwrap().len() > 0 {

                    let wbatch = cloned_self.check({
                            cloned_self.cache_heap.lock().unwrap().pop().unwrap()
                        })
                        .expect("Step 1 failed. Checks not passed");
                    
                    cloned_self.single_chr(wbatch)
                        .expect("Step 2 failed. Failed to partition chromosomes")
                        .into_iter()
                        .map(
                            |b| cloned_self.align(b).expect("Step 3 failed. Failed to align")
                        ).collect_vec().into_iter()
                        .map(
                            |b| cloned_self.cast_chromosomes(b).expect("Step 4 failed. Failed to cast chromosome col")
                        ).map(
                            |b| cloned_self.final_check(b).expect("Step 5 failed. Final check not passed")
                        )
                        .for_each(
                            |mut b| {
                                let un = b.data_mut();
                                let data = unsafe { ManuallyDrop::take(&mut un.bsx) };
                                cloned_self.res_heap.lock().unwrap().push(data);
                            }
                        );
                }
            }
        })
    }

    /// [BatchStatus::ChrCasted] ==> [BatchStatus::FullChecked]
    ///
    /// Performs transoframtion to [ReadFinalBatch] and performs checks:
    /// -   [BSXBatch::check_chr_unique]
    /// -   [BSXBatch::check_position_sorted]
    /// -   [BSXBatch::check_types]
    /// -   [BSXBatch::check_cols]
    ///
    /// Puts the final item back to the channel
    fn final_check(&self, mut item: WorkBatch) -> BSXResult<WorkBatch> {
        let data = unsafe { ManuallyDrop::take(&mut item.data_mut().df) };
        let formatted = report_type_to_bsx(&item.report_type(), data)?;
        let record_batch = ReadFinalBatch::from_df(formatted);

        record_batch.check_cols()?;
        record_batch.check_types()?;
        record_batch.check_chr_unique()?;
        record_batch.check_position_sorted()?;

        let item = WorkBatch::new(
            ReadQueueData::new_bsx(record_batch),
            BatchStatus::FullChecked,
            item.report_type(),
        );
        debug!("Final check passed for {}", item.uuid());
        Ok(item)
    }

    /// [BatchStatus::Aligned] ==> [BatchStatus::ChrCasted]
    /// Casts chromosome column to [DataType::Categorical] with preset
    /// possible chromosomes and their indices
    fn cast_chromosomes(&self, mut item: WorkBatch) -> BSXResult<WorkBatch> {
        let df = unsafe { ManuallyDrop::take(&mut item.data_mut().df) };
        let data = df
            .lazy().cast(PlHashMap::from_iter(vec![(BSXCol::Chr.as_str(), self.chr_type.clone())]), true)
            .collect()?
            .to_owned();
        let item = WorkBatch::new(
            ReadQueueData::new_df(data),
            BatchStatus::ChrCasted,
            item.report_type(),
        );
        debug!("Cast chromosomes passed for {}", item.uuid());
        Ok(item)
    }

    /// [BatchStatus::PartitionedChr] ==> [BatchStatus::Aligned]
    /// Aligns data to the reference genome cytosines
    fn align(&mut self, batch: WorkBatch) -> BSXResult<WorkBatch> {
        let report_type = batch.report_type();
        let coordinates = RegionCoordinates {
            start: unsafe { &batch.data().df }
                .column("position")?
                .u32()?
                .first()
                .unwrap(),
            end: unsafe { &batch.data().df }
                .column("position")?
                .u32()?
                .last()
                .unwrap(),
            chr: unsafe { &batch.data().df }
                .column("chr")?
                .str()?
                .first()
                .unwrap()
                .to_owned(),
        };

        let context_df = self.get_context_df(coordinates.clone())?;
        let mut old_data = unsafe { ManuallyDrop::take(&mut batch.into_data().df) };
        
        old_data = old_data.drop(BSXCol::Context.as_str())?;
        old_data = old_data.drop(BSXCol::Strand.as_str())?;
        
        let new = context_df
            .lazy()
            .join(
                old_data.lazy(),
                [BSXCol::Position.into()],
                [BSXCol::Position.into()],
                JoinArgs::new(JoinType::Left),
            )
            .with_column(
                Expr::from(BSXCol::Chr).fill_null(lit(coordinates.chr)),
            )
            .collect()?;

        let item = WorkBatch::new(
            ReadQueueData::new_df(new),
            BatchStatus::Aligned,
            report_type,
        );
        debug!("Aligned data for {}", item.uuid());
        Ok(item)
    }

    /// [BatchStatus::ColumnsChecked] ==> [BatchStatus::PartitionedChr]
    /// Data is being partitioned by chromosome
    fn single_chr(&self, mut item: WorkBatch) -> BSXResult<Vec<WorkBatch>> {
        let item = unsafe { ManuallyDrop::take(&mut item.data_mut().df) }
            .partition_by_stable(["chr"], true)?
            .into_iter()
            .map(|data| {
                let item = WorkBatch::new(
                    ReadQueueData::new_df(data),
                    BatchStatus::PartitionedChr,
                    item.report_type(),
                );
                debug!("Single chr for {}", item.uuid());
                item
            }).collect_vec();
        Ok(item)
    }

    /// [BatchStatus::InitialRaw] ==> [BatchStatus::ColumnsChecked]
    fn check(&self, mut item: WorkBatch) -> BSXResult<WorkBatch> {
        let mut df = unsafe { ManuallyDrop::take(&mut item.data_mut().df) };
        if !df.schema().try_get("chr")?.is_string() {
            return Err(Box::from("Chr col missing"));
        }
        if !df.schema().try_get("position")?.is_numeric() {
            return Err(Box::from("Position must be numeric"));
        }
        
        df = df.lazy()
            .cast(PlHashMap::from_iter(vec![(BSXCol::Position.as_str(), BSXCol::Position.data_type())]), true)
            .collect()?;
        let new = WorkBatch::new(
            ReadQueueData::new_df(df), BatchStatus::ColumnsChecked, item.report_type(),
        );
        debug!("Check passed for {}", item.uuid());
        Ok(new)
    }

    /// Use FASTA [IndexedReader] to retrieve region sequence and convert it
    /// to methylation contexts table
    fn get_context_df(&mut self, region_coordinates: RegionCoordinates) -> BSXResult<DataFrame> {
        // Get reader ref
        let mut reader = self.fasta_reader.lock().unwrap();

        // Get chromosome length from index
        let chr_len = match reader
            .index
            .sequences()
            .iter()
            .find(|s| s.name == *region_coordinates.chr)
        {
            Some(chr_data) => chr_data.len,
            None => {
                return Err(Box::<dyn Error>::from(format!(
                    "chr {} not found in index",
                    region_coordinates.chr
                )))
            }
        };
        let expanded_coords = region_coordinates.clone().expand(2, chr_len as u32);
        // Read sequence
        let mut seq = Vec::new();
        reader.fetch(
            region_coordinates.chr.as_str(),
            expanded_coords.start as u64,
            expanded_coords.end as u64,
        )?;
        reader.read(&mut seq)?;

        // Convert to DataFrame
        let (positions, contexts, strands) =
            Self::parse_cytosines(String::from_utf8(seq)?, expanded_coords);
        Ok(DataFrame::from_iter([
            Column::new("position".into(), positions),
            Column::new("context".into(), contexts),
            Column::new("strand".into(), strands),
        ]))
    }
    fn parse_cytosines(
        seq: String,
        region_coordinates: RegionCoordinates,
    ) -> (Vec<u32>, Vec<Option<bool>>, Vec<bool>) {
        let length = region_coordinates.length() as usize;
        let (fw_bound, rv_bound) = (2, length - 2);

        let mut positions = Vec::<u32>::with_capacity(length);
        let mut contexts = Vec::<Option<bool>>::with_capacity(length);
        let mut strands = Vec::<bool>::with_capacity(length);
        let uppercased = seq.to_ascii_uppercase();
        let ascii_seq = uppercased.as_bytes();

        'seq_iter: for (index, nuc) in ascii_seq.iter().enumerate() {
            let forward = match nuc {
                b'C' => true,
                b'G' => false,
                _ => continue 'seq_iter,
            };

            let skip_cond = |strand: bool, index: usize| {
                let cmp = index.cmp(if strand { &fw_bound } else { &rv_bound });
                if !strand {
                    cmp.reverse()
                } else {
                    cmp
                }
            };

            if skip_cond(forward, index) == Ordering::Less {
                continue 'seq_iter;
            } else {
                let target_nuc = if forward { b'G' } else { b'C' };
                let k: isize = if forward { 1 } else { -1 };

                let context = if ascii_seq[((index as isize + k) as usize)] == target_nuc {
                    Some(true)
                } else if ascii_seq[(index as isize + 2 * k) as usize] == target_nuc {
                    Some(false)
                } else {
                    None
                };
                positions.push(region_coordinates.start + index as u32);
                contexts.push(context);
                strands.push(forward);
            }
        }
        positions.shrink_to_fit();
        contexts.shrink_to_fit();
        strands.shrink_to_fit();
        (positions, contexts, strands)
    }
}

impl Iterator for ParserThread {
    type Item = ReadFinalBatch;
    fn next(&mut self) -> Option<Self::Item> {
        if self.res_heap.lock().unwrap().len() == 0 && !self.receiver.is_empty() {
            sleep(std::time::Duration::from_millis(10));
            self.next()
        } else if self.res_heap.lock().unwrap().len() == 0 {
            return None;
        } else {
            self.res_heap.lock().unwrap().pop()
        }
    }
}

#[derive(Clone)]
struct ReaderThread {
    batched_reader: Arc<Mutex<crate::io::report::read::OwnedBatchedCsvReader>>,
    batch_per_read: usize,
    sender: Arc<Sender<WorkBatch>>,
    report_type: ReportType
}

impl ReaderThread {
    pub fn new(
        mut csv_reader: CsvReader<Box<dyn MmapBytesReader>>, 
        batch_per_read: usize, 
        sender: Arc<Sender<WorkBatch>>, 
        report_type: ReportType
    ) -> Self {
        let batched_reader = csv_reader.batched_borrowed().unwrap();
        let batched_reader: BatchedCsvReader<'static> =
            unsafe { std::mem::transmute(batched_reader) };
        let schema = SchemaRef::from(report_type.get_schema());

        let batched_owned = crate::io::report::read::OwnedBatchedCsvReader {
            schema,
            batched_reader,
            _reader: csv_reader,
        };
        Self { 
            batched_reader: Arc::from(Mutex::new(batched_owned)), 
            batch_per_read, 
            sender,
            report_type,
        }
    }
    
    pub fn run(&mut self) -> JoinHandle<()> {
        let cloned_self = self.clone();
        spawn(move || {
            let mut reader = cloned_self.batched_reader.lock().unwrap();
            
            info!("ReaderThread started");
            
            while let Some(batches) = reader.next_batches(cloned_self.batch_per_read).unwrap() {
                for df in batches {
                    let data = UBatch::from_report_type(df, &cloned_self.report_type).unwrap();
                    let work_batch = WorkBatch::new(
                        ReadQueueData::new_df(data.into()),
                        BatchStatus::InitialRaw,
                        cloned_self.report_type,
                    );
                    cloned_self.sender.send(work_batch).unwrap();
                }
            }
        })
    }

    pub fn batch_per_read(&self) -> usize {
        self.batch_per_read
    }

    pub fn report_type(&self) -> ReportType {
        self.report_type
    }
}

pub struct ReportReader {
    reader: ReaderThread,
    parser: ParserThread,
    active_batch: Option<ReadFinalBatch>,
    chunk_size: usize,
    join_handles: Vec<JoinHandle<()>>,
}

impl Iterator for ReportReader {
    type Item = ReadFinalBatch;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(batch) = self.active_batch.take() {
            let (ret, left) = batch.split(self.chunk_size as i64);
            self.active_batch = Some(left);
            Some(ret)
        } else {
            match self.parser.next() {
                Some(batch) => {
                    self.active_batch = Some(batch);
                    self.next()
                },
                None => {
                    for handle in self.join_handles.drain(..) {
                        handle.join().unwrap();
                    }
                    None
                },
            }
        }
    }
}

pub struct ReportReaderBuilder {
    report_type: ReportType,
    report_path: PathBuf,
    fasta_path: PathBuf,
    fai_path: PathBuf,
    batch_per_read: usize,
    batch_size: usize,
    chunk_size: usize,
    low_memory: bool,
    n_threads: Option<usize>,
    rechunk: bool
}

impl Default for ReportReaderBuilder {
    fn default() -> Self {
        Self {
            report_type: ReportType::BISMARK,
            report_path: PathBuf::new(),
            fasta_path: PathBuf::new(),
            fai_path: PathBuf::new(),
            batch_per_read: 16,
            batch_size: 2 << 20,
            chunk_size: 10_000,
            low_memory: false,
            n_threads: None,
            rechunk: true
        }
    }
}

impl ReportReaderBuilder {

    pub fn with_report_type(self, report_type: ReportType) -> Self {
        Self { report_type, ..self }
    }
    pub fn with_report_path(self, report_path: PathBuf) -> Self {
        Self { report_path, ..self }
    }
    pub fn with_fasta_path(self, fasta_path: PathBuf) -> Self {
        Self { fasta_path, ..self }
    }
    pub fn with_fai_path(self, fai_path: PathBuf) -> Self {
        Self { fai_path, ..self }
    }
    pub fn with_batch_size(self, batch_size: usize) -> Self {
        Self { batch_size, ..self }
    }
    pub fn with_chunk_size(self, chunk_size: usize) -> Self {
        Self { chunk_size, ..self }
    }
    pub fn with_batch_per_read(self, batch_per_read: usize) -> Self {
        Self { batch_per_read, ..self }
    }
    
    pub fn build(self) -> BSXResult<ReportReader> {
        if !self.report_path.exists() {
            return Err(Box::from(format!("Report path does not exist: {}", self.report_path.display())));
        }
        if !self.fasta_path.exists() {
            return Err(Box::from(format!("FASTA path does not exist: {}", self.report_path.display())));
        }
        if !self.fai_path.exists() {
            return Err(Box::from(format!("FAI path does not exist: {}", self.report_path.display())));
        }
        
        let csv_reader: CsvReader<Box<dyn MmapBytesReader>> = self.report_type.get_reader(
            Box::new(File::open(self.report_path)?),
            Some(self.chunk_size), Some(self.low_memory), self.n_threads, Some(self.rechunk)
        );
        
        let fasta_reader = IndexedReader::new(
            File::open(self.fasta_path)?,
            File::open(self.fai_path)?,
        )?;
        let chr_names = fasta_reader.index
            .sequences().iter()
            .map(|s| Some(s.name.clone()))
            .collect_vec();
        let categories = DataType::Enum(
            Some(Arc::new(RevMapping::build_local({
                Utf8ViewArray::from_slice(chr_names.as_slice())
            }))),
            CategoricalOrdering::Physical,
        );

        let (sender, reciever) = mpmc::sync_channel(self.batch_per_read);
        let mut reader_thread = ReaderThread::new(
            csv_reader, self.batch_per_read, Arc::new(sender.clone()), self.report_type,
        );
        let mut parser_thread = ParserThread::new(
            fasta_reader, Arc::new(reciever), Arc::new(sender), categories 
        );
        
        let join_handle = vec![reader_thread.run(), parser_thread.run()];
        
        Ok( ReportReader {
            reader: reader_thread, parser: parser_thread, active_batch: None, chunk_size: self.chunk_size,
            join_handles: join_handle
        } )
    }
}