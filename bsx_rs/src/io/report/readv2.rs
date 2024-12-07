use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::sync::mpsc;
use std::thread;
use std::thread::{JoinHandle};
use bio::io::fasta::IndexedReader;
use log::{debug, info};
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::io::mmap::MmapBytesReader;
use polars::prelude::{*};
use polars_arrow::array::Utf8ViewArray;
use crate::io::report::read::parse_cytosines;
use crate::io::report::types::ReportType;
use crate::region::RegionCoordinates;
use crate::ubatch2::{BSXBatch, BSXCol};

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

impl BSXBatch for IncomingBatch {
    fn from_df(data_frame: DataFrame) -> Self {
        Self {
            data: data_frame,
            report_type: ReportType::BISMARK
        }
    }

    fn get_data(&self) -> &DataFrame {
        &self.data
    }

    fn get_data_mut(&mut self) -> &mut DataFrame {
        &mut self.data
    }

    fn _set_data(&mut self, data_frame: DataFrame) {
        debug!("Incoming batch for {}", data_frame);
        self.data = data_frame;
    }
}

/// Unlike [UBatch], which already expects data to be casted to
/// proper types, [IncomingBatch::new] performs casting and [DataType::Categorical]
/// creation for [BSXCol::Chr]
impl IncomingBatch {
    fn new(data_frame: DataFrame, report_type: &ReportType, chr_names: &Vec<String>) -> Self {

        let categories = DataType::Enum(
            Some({
                let mut cat_builder =
                    CategoricalChunkedBuilder::new(Default::default(), 0, Default::default());
                for chr_name in chr_names {
                    let _ = &cat_builder.append(Some(chr_name.as_str()));
                }
                cat_builder.finish().get_rev_map().clone()
            }),
            CategoricalOrdering::Physical,
        );

        assert!(categories.contains_categoricals());
        
        let obj = Self::from_report_type(data_frame, report_type).unwrap();
        
        if let Err(e) = obj.check_cols() {
            panic!("Error while checking cols: {}", e);
        }
        let data: DataFrame = obj.into();
        IncomingBatch::from_cast(data.lazy(), chr_names).expect("Error casting")
            .with_report_type(report_type.clone())
    }

    fn with_report_type(mut self, report_type: ReportType) -> Self {
        self.report_type = report_type;
        self
    }
}
impl Into<DataFrame> for IncomingBatch {
    fn into(self) -> DataFrame {
        self.data
    }
}

/// Wrapper for [DataFrame], which implements [BSXBatch].
/// Differences from other implementation are:
///     1. Check for single chromosome data is performed
///     2. Check if data is sorted is performed
///     3. Compare methods are implemented
pub struct ReadBatch {
    data: DataFrame,
}

impl ReadBatch {
    pub fn new(data: DataFrame) -> Self {
        let mut obj = Self::from_df(data);

        if let Err(e) = obj.check_cols() {
            panic!("{:?}", e)
        }
        if let Err(e) = obj.check_types() {
            panic!("{:?}", e)
        }
        if !obj.check_chr_unique() {
            panic!("Chromosomes differ!")
        };
        if !obj.check_position_sorted() {
            obj.sort_positions();
        }
        obj
    }
}

impl From<ReadBatch> for DataFrame {
    fn from(item: ReadBatch) -> Self {
        item.data
    }
}


impl BSXBatch for ReadBatch {
    fn from_df(data_frame: DataFrame) -> Self {
        let obj = Self { data: data_frame };
        obj
    }
    fn get_data(&self) -> &DataFrame {
        &self.data
    }
    fn get_data_mut(&mut self) -> &mut DataFrame {
        &mut self.data
    }

    fn _set_data(&mut self, data_frame: DataFrame) {
        self.data = data_frame;
    }
}

impl PartialEq<Self> for ReadBatch {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl Eq for ReadBatch {}

impl PartialOrd for ReadBatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.get_chr_idx().cmp(&other.get_chr_idx()).then(self.get_first_pos().cmp(&other.get_last_pos())))
    }
}

impl Ord for ReadBatch {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

pub struct ReportReaderBuilder {
    batch_per_read: usize,
    max_cache_size: usize,
    chunk_size: usize,
    read_batch_size: usize,
    low_memory: bool,
    rechunk: bool,
    n_threads: usize,
    align: bool
}

impl Default for ReportReaderBuilder {
    fn default() -> Self {
        Self {
            batch_per_read: 30,
            max_cache_size: 100,
            chunk_size: 10_000,
            read_batch_size: 100 * 2 << 10,
            low_memory: false,
            rechunk: true,
            n_threads: usize::from(thread::available_parallelism().unwrap()),
            align: false
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
    pub fn with_align(mut self, align: bool) -> Self {
        self.align = align;
        self
    }
    
    pub fn finish(
        self,
        report_file: PathBuf,
        report_type: ReportType,
        fasta_file: PathBuf,
        fasta_index_file: PathBuf,
    ) -> Result<ReportReader, Box<dyn Error>> {
        let fasta_reader = IndexedReader::new(
            File::open(fasta_file.clone())?,
            File::open(fasta_index_file.clone())?,
        )?;
        info!("Opened Fasta reader from {fasta_file:?} with index {fasta_index_file:?}");
        let chr_revmap = Arc::new(RevMapping::build_local(
            Utf8ViewArray::arr_from_iter(fasta_reader.index.sequences().iter().map(|seq| seq.name.as_str()))
        ));

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
            let mapping_clone = chr_revmap.clone();
            let report_type_clone = report_type.clone();
            while let Some(batches) = batched_owned.next_batches(self.batch_per_read).unwrap() {
                info!("Read {} batches", batches.len());
                batches.into_iter()
                    .map(|data| IncomingBatch::new(data, &report_type_clone, &mapping_clone))
                    .for_each(|batch| sc.send(batch).unwrap());
            };
        });
        
        info!("Created ReportReader on {report_file:?}");
        Ok({
            ReportReader {
                reading_thread,
                rc,
                batch_cache: BinaryHeap::new(),
                align: self.align,
                fasta_reader,
                read_batch_size: self.read_batch_size,
                batch_per_read: self.batch_per_read,
                chunk_size: self.chunk_size,
            }
        })
    }
}

pub struct ReportReader {
    reading_thread: JoinHandle<()>,
    
    rc: mpsc::Receiver<IncomingBatch>,

    batch_cache: BinaryHeap<Reverse<ReadBatch>>,
    fasta_reader: IndexedReader<File>,

    // Config
    batch_per_read: usize,
    chunk_size: usize,
    read_batch_size: usize,
    align: bool
}

impl ReportReader {
    /// Formats incoming data to [ReadBatch] format, performing checks and
    /// splitting batch data by chromosome.
    fn format_data(&mut self, data: IncomingBatch) -> Result<Vec<ReadBatch>, Box<dyn Error>> {
        debug!("Formatting data with {} rows", data.data.height());
        let converted = if data.check_chr_unique() {
            vec![ReadBatch::new(data.into())]
        } else {
            data.partition_by_bsx(BSXCol::Chr, true).into_iter().map(|b| ReadBatch::new(b.into())).collect()
        };

        Ok(converted)
    }
    
    /// Aligns data to reference sequence
    fn align_with_reference(&mut self, batch: ReadBatch) -> Result<ReadBatch, Box<dyn Error>> {
        let coordinates = batch.get_region_coordinates().unwrap();
        debug!("Aligning data for region {coordinates:?}");
        let context_df = self.get_context_df(coordinates)?;
        let new = context_df.lazy()
            .join(
                DataFrame::from(batch).lazy(),
                [col("position")], [col("position")],
                JoinArgs::new(JoinType::Left),
        ).collect()?;
        Ok(ReadBatch::from_df(new))
    }
    
    /// Make [ReadBatch] from incoming batch. Calls [ReportReader::format_data] and
    /// [ReportReader::align_with_reference] if FASTA reader is specified or when 
    /// report_type is [ReportType::BEDGRAPH] or [ReportType::COVERAGE].
    fn process_data(&mut self, data: IncomingBatch) -> Result<(), Box<dyn Error>> {
        let report_type = data.report_type.clone();
        let formatted = self.format_data(data)?;
        let aligned = if [ReportType::BEDGRAPH, ReportType::COVERAGE].contains(&report_type)
            || self.align {
            formatted.into_iter().map(|b| self.align_with_reference(b)).collect::<Result<Vec<_>, _>>()?
        } else {
            formatted
        };
        aligned.into_iter().for_each(|item| {self.batch_cache.push(Reverse(item));});
        Ok(())
    }

    /// Use FASTA [IndexedReader] to retrieve region sequence and convert it
    /// to methylation contexts table
    fn get_context_df(&mut self, region_coordinates: RegionCoordinates) -> Result<DataFrame, Box<dyn Error>> {
        // Get reader ref
        let reader = &mut self.fasta_reader;

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
    type Item = Result<ReadBatch, Box<dyn Error>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.batch_cache.len() < self.batch_per_read && !self.reading_thread.is_finished() {
            let new_batch = self.rc.recv().unwrap();
            if self.process_data(new_batch).is_err() {
                return Some(Err(Box::from("Error processing batch")));
            };
            return self.next();
        };
        
        if !self.reading_thread.is_finished() || !self.batch_cache.is_empty() {
            let current = self.batch_cache.pop().unwrap().0;

            let res = match current.length() as isize - self.chunk_size as isize {
                0 => current,
                ..0 => {
                    match self.batch_cache.pop() {
                        Some(val) if val.0.get_chr_idx() == current.get_chr_idx() => {
                            let next = val.0;
                            let (out, back) = current.vstack(&next).split(self.chunk_size as i64);
                            self.batch_cache.push(Reverse(back));
                            out
                        },
                        _ => current
                    }

                }
                1.. => {
                    let (out, back) = current.split(self.chunk_size as i64);
                    self.batch_cache.push(Reverse(back));
                    out
                }
            };
            
            Some(Ok(res))
        } else {
            None
        }
    }
}
