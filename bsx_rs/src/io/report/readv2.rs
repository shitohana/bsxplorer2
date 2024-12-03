use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::ops::Deref;
use std::path::PathBuf;
use std::sync::mpsc;
use std::thread;
use std::thread::{JoinHandle};
use bio::io::fasta::IndexedReader;
use itertools::Itertools;
use log::{debug, info};
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::io::mmap::MmapBytesReader;
use polars::prelude::{*};
use crate::io::report::read::parse_cytosines;
use crate::io::report::types::ReportType;
use crate::region::RegionCoordinates;
use crate::ubatch::UniversalBatch;

pub(crate) struct OwnedBatchedCsvReader {
    #[allow(dead_code)]
    // this exists because we need to keep ownership
    pub schema: SchemaRef,
    pub batched_reader: BatchedCsvReader<'static>,
    // keep ownership
    pub _reader: CsvReader<Box<dyn MmapBytesReader>>,
}

impl OwnedBatchedCsvReader {
    pub fn next_batches(&mut self, n: usize) -> PolarsResult<Option<Vec<DataFrame>>> {
        self.batched_reader.next_batches(n)
    }
}

struct IncomingBatch {
    data: DataFrame,
    report_type: ReportType
}

struct ReadQueueItem {
    data: UniversalBatch,
    chr_idx: usize,
}

impl PartialEq<Self> for ReadQueueItem {
    fn eq(&self, other: &Self) -> bool {
        self.chr_idx == other.chr_idx 
            && self.data.first_position() == other.data.first_position()
            && self.data.last_position() == other.data.last_position()
            && self.data.count() == other.data.count()
    }
}

impl PartialOrd for ReadQueueItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.chr_idx.cmp(&other.chr_idx).then(self.data.first_position().cmp(&other.data.first_position())))
    }
}

impl Eq for ReadQueueItem {}

impl Ord for ReadQueueItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

pub struct ReportReaderBuilder {
    batch_per_read: usize,
    max_cache_size: usize,
    chunk_size: usize,
    read_batch_size: usize,
    fasta_file: Option<PathBuf>,
    fasta_index_file: Option<PathBuf>,
    low_memory: bool,
    rechunk: bool,
    n_threads: usize,
}

impl Default for ReportReaderBuilder {
    fn default() -> Self {
        Self {
            batch_per_read: 30,
            max_cache_size: 100,
            chunk_size: 10_000,
            read_batch_size: 100 * 2 << 10,
            fasta_file: None, 
            fasta_index_file: None,
            low_memory: false, 
            rechunk: true,
            n_threads: usize::from(std::thread::available_parallelism().unwrap()),
        }
    }
}

impl ReportReaderBuilder {
    pub fn with_batch_per_read(mut self, n: usize) -> Self {
        self.batch_per_read = n;
        self
    }
    pub fn with_chunk_size(mut self, n: usize) -> Self {
        self.chunk_size = n;
        self
    }
    pub fn with_read_batch_size(mut self, n: usize) -> Self {
        self.read_batch_size = n;
        self
    }
    pub fn with_fasta_file(mut self, path: PathBuf) -> Self {
        self.fasta_file = Some(path);
        self
    }
    pub fn with_fasta_index_file(mut self, path: PathBuf) -> Self {
        self.fasta_index_file = Some(path);
        self
    }
    pub fn with_low_memory(mut self, low: bool) -> Self {
        self.low_memory = low;
        self
    }
    pub fn with_rechunk(mut self, rechunk: bool) -> Self {
        self.rechunk = rechunk;
        self
    }
    pub fn with_n_threads(mut self, n: usize) -> Self {
        self.n_threads = n;
        self
    }
    
    pub fn finish(self, report_file: PathBuf, report_type: ReportType) -> Result<ReportReader, Box<dyn Error>> {
        let fasta_reader = if let (Some(fasta_file), Some(fasta_index)) = (self.fasta_file, self.fasta_index_file) {
            let reader = IndexedReader::new(
                File::open(fasta_file.clone())?,
                File::open(fasta_index.clone())?,
            )?;
            info!("Opened Fasta reader from {fasta_file:?} with index {fasta_index:?}");
            Some(reader)
        } else {None};
        
        let mut reader: CsvReader<Box<dyn MmapBytesReader>> = report_type.clone().get_reader(
            Box::new(BufReader::new(File::open(report_file.clone())?)),
            Some(self.read_batch_size), Some(self.low_memory), Some(self.n_threads), Some(self.rechunk)
        );
        let schema = SchemaRef::from(report_type.clone().get_schema());
        let batched_reader = reader.batched_borrowed().unwrap();
        let batched_reader: BatchedCsvReader<'static> =
            unsafe { std::mem::transmute(batched_reader) };
        let mut batched_owned = OwnedBatchedCsvReader {
                schema,
                batched_reader,
                _reader: reader,
            };
        
        let (sc, rc) = mpsc::sync_channel::<IncomingBatch>(self.max_cache_size);
        
        let reading_thread = thread::spawn(move || -> () {
            while let Some(batches) = batched_owned.next_batches(self.batch_per_read).unwrap() {
                info!("Read {} batches", batches.len());
                batches.into_iter()
                    .map(|data| IncomingBatch { data, report_type: report_type.clone() })
                    .for_each(|batch| sc.send(batch).unwrap());
            };
        });
        
        info!("Created ReportReader on {report_file:?}");
        Ok({
            ReportReader {
                reading_thread,
                rc,
                batch_cache: BinaryHeap::new(),
                fasta_reader,
                read_batch_size: self.read_batch_size,
                batch_per_read: self.batch_per_read,
                chunk_size: self.chunk_size,
                observed_chroms: vec![],
            }
        })
    }
}

pub struct ReportReader {
    reading_thread: JoinHandle<()>,
    
    rc: mpsc::Receiver<IncomingBatch>,
    
    batch_cache: BinaryHeap<Reverse<ReadQueueItem>>,
    fasta_reader: Option<IndexedReader<File>>,
    observed_chroms: Vec<String>,

    // Config
    batch_per_read: usize,
    chunk_size: usize,
    read_batch_size: usize,
}

impl ReportReader {
    /// Formats incoming data to UniversalBatch format, performing checks and
    /// splitting batch data by chromosome. Observed chromosomes are written into
    /// observed_chroms field.
    fn format_data(&mut self, data: IncomingBatch) -> Result<Vec<UniversalBatch>, Box<dyn Error>> {
        debug!("Formatting data with {} rows", data.data.height());
        let converted = UniversalBatch::from_report_type(data.data, &data.report_type)?;
        let is_border = !converted.unique_chr();
        
        let mut universal_checked = if is_border {
            converted.partition_chr()
        } else {
            vec![converted]
        };
        universal_checked.iter_mut().for_each(|batch| {batch.check_sorted().unwrap()});
        universal_checked.iter().for_each(|batch| {self.update_observed_chroms(vec![batch.get_chr()])});
        Ok(universal_checked)
    }
    
    /// Aligns data to reference sequence
    fn align_with_reference(&mut self, batch: UniversalBatch) -> Result<UniversalBatch, Box<dyn Error>> {
        let coordinates = batch.get_region_coordinates();
        debug!("Aligning data for region {coordinates:?}");
        let context_df = self.get_context_df(coordinates)?;
        Ok({
            context_df.lazy()
                .join(
                    DataFrame::from(batch).lazy(),
                    [col("position")], [col("position")],
                    JoinArgs::new(JoinType::Left),
                ).collect()?.into()
        })
    }
    
    /// Make [ReadQueueItem] from incoming batch. Calls [ReportReader::format_data] and 
    /// [ReportReader::align_with_reference] if FASTA reader is specified or when 
    /// report_type is [ReportType::BEDGRAPH] or [ReportType::COVERAGE].
    fn process_data(&mut self, data: IncomingBatch) -> Result<(), Box<dyn Error>> {
        let report_type = data.report_type.clone();
        let formatted = self.format_data(data)?;
        let aligned = if [ReportType::BEDGRAPH, ReportType::COVERAGE].contains(&report_type)
            || self.fasta_reader.is_some() {
            formatted.into_iter().map(|b| self.align_with_reference(b)).collect::<Result<Vec<_>, _>>()?
        } else if [ReportType::BEDGRAPH, ReportType::COVERAGE].contains(&report_type) && 
            self.fasta_reader.is_none() {
            return Err(Box::from("FASTA reader must be specified for BEDGRAPH or COVERAGE reports reading"))
        } else {
            formatted
        };
        
        let chunk_size = self.chunk_size.clone();
        let items = aligned.into_iter()
            .map(|b| b.partition_by_chunks(chunk_size))
            .flatten()
            .map(|b| Reverse(self.make_item(b)))
            .collect_vec();
        items.into_iter().for_each(|item| {self.batch_cache.push(item);});
        Ok(())
    }
    
    fn make_item(&self, batch: UniversalBatch) -> ReadQueueItem {
        let batch_chr = batch.get_chr();
        ReadQueueItem {
            data: batch,
            chr_idx: self.observed_chroms.iter().position(|chr| *chr == batch_chr).unwrap()
        }
    }

    /// Get reference to FASTA [IndexedReader]
    fn get_fasta_reader(&mut self) -> Option<&mut IndexedReader<File>> {
        self.fasta_reader.as_mut()
    }
    
    /// Get observed chroms
    fn get_chroms(&self) -> &Vec<String> {
        &self.observed_chroms
    }
    
    /// Get observed chroms mutable
    fn get_chroms_mut(&mut self) -> &mut Vec<String> {
        self.observed_chroms.as_mut()
    }

    /// Update observed chroms with array of chromosomes
    fn update_observed_chroms(&mut self, new_chroms: Vec<String>) {
        let chroms = self.get_chroms_mut();
        new_chroms.into_iter().for_each(|chrom| {if !chroms.contains(&chrom) {chroms.push(chrom);}})
    }
    
    /// Use FASTA [IndexedReader] to retrieve region sequence and convert it
    /// to methylation contexts table
    fn get_context_df(&mut self, region_coordinates: RegionCoordinates) -> Result<DataFrame, Box<dyn Error>> {
        // Get reader ref
        let reader = match self.get_fasta_reader() {
            Some(reader) => reader,
            None => return Err(Box::<dyn Error>::from("FASTA reader should be present"))
        };

        // Get chromosome length from index
        let chr_len = match reader.index.sequences().iter().find(|s| s.name == *region_coordinates.chr) {
            Some(chr_data) => chr_data.len,
            None => return Err(Box::<dyn Error>::from(format!("chr {} not found in index", region_coordinates.chr)))
        };

        // Expand read bounds if possible
        let start = if region_coordinates.start > 2 { region_coordinates.start - 2 } else { region_coordinates.start };
        let end = if region_coordinates.end + 2 <= chr_len as u32 { region_coordinates.end + 2} else { chr_len as u32 };

        // Read sequence
        let mut seq = Vec::new();
        reader.fetch(region_coordinates.chr.as_str(), start as u64, end as u64)?;
        reader.read(&mut seq)?;

        // Convert to DataFrame
        let (positions, contexts, strands) = parse_cytosines(String::from_utf8(seq)?, start);
        Ok(DataFrame::from_iter([
            Column::new("position".into(), positions),
            Column::new("context".into(), contexts),
            Column::new("strand".into(), strands)
        ]))
    }
}

impl Iterator for ReportReader {
    type Item = Result<UniversalBatch, Box<dyn Error>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.batch_cache.len() < self.batch_per_read && !self.reading_thread.is_finished() {
            let new_batch = self.rc.recv().unwrap();
            if self.process_data(new_batch).is_err() {
                return Some(Err(Box::from("Error processing batch")));
            };
            return self.next();
        };
        
        if !self.reading_thread.is_finished() || !self.batch_cache.is_empty() {
            let next = self.batch_cache.pop().unwrap().0.data;
            
            if next.count() != self.chunk_size {
                let peek = self.batch_cache.peek().unwrap();
            } 
            
            Some(Ok(next))
        } else {
            None
        }
    }
}
