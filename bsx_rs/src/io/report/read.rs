use std::cmp::{Ordering, Reverse};
use crate::io::report::types::ReportType;
use crate::ubatch::UniversalBatch;
use itertools::Itertools;
use log::{debug, info};
use polars::error::PolarsResult;
use polars::io::mmap::MmapBytesReader;
use polars::prelude::*;
use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};
use std::collections::hash_map::Drain;
use std::fs::File;
use rayon::prelude::*;


pub(crate) struct OwnedBatchedCsvReader {
    #[allow(dead_code)]
    // this exists because we need to keep ownership
    schema: SchemaRef,
    batched_reader: BatchedCsvReader<'static>,
    // keep ownership
    _reader: CsvReader<Box<dyn MmapBytesReader>>,
}

impl OwnedBatchedCsvReader {
    pub(crate) fn next_batches(&mut self, n: usize) -> PolarsResult<Option<Vec<DataFrame>>> {
        self.batched_reader.next_batches(n)
    }
}

/// Parse cytosine sequence and extract cytosine methylation contexts and positions.
/// # Returns
/// Vector of \[position, context, strand\]
fn parse_cytosines(seq: String, start_pos: u32) -> (Vec<u32>, Vec<Option<bool>>, Vec<bool>) {
    let start_pos = start_pos - 1;
    let fw_bound: usize = seq.len() - 2;
    let rv_bound: usize = 2;

    let mut positions: Vec<u32> = Vec::new();
    let mut contexts: Vec<Option<bool>> = Vec::new();
    let mut strands: Vec<bool> = Vec::new();
    let uppercased = seq.to_uppercase();
    let ascii_seq = uppercased.as_bytes();

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
        positions.push(start_pos + index as u32);
        contexts.push(context);
        strands.push(forward);
    }
    (positions, contexts, strands)
}

/// Stores information about batch genomic position
/// Tuple of (chromosome name, min position, max position)
pub struct BatchStats(String, u64, u64);

impl BatchStats {
    
}

/// Iterator over methylation report data. Output batches
/// are sorted by genomic position.
///
/// **Input** data is expected to be **sorted**.

pub struct ReportReader {
    /// [ReportType] of the report
    report_type: ReportType,
    /// Initialized [BatchedCsvReader]. Use [ReportType::get_reader] to create binding for
    /// reader and [CsvReader::batched_borrowed] to create Batched reader.
    batched_reader: OwnedBatchedCsvReader,
    /// How many batches will be read into memory. Should be >= num cores.
    batch_per_read: usize,
    // TODO: implement for faster writing
    /// Number of rows in the output batches. All batches will be same specified
    /// length if possible.
    batch_size: usize,
    /// If [ReportType] is [ReportType::COVERAGE] or [ReportType::BEDGRAPH], this option
    /// should be specified with initialized FASTA reader. If other report types -
    /// it is optional to apply additional alignment to reference cytosines.
    fasta_reader: bio::io::fasta::IndexedReader<File>,
    /// Queue of computed batches.
    read_queue: ReadQueue,
}

impl ReportReader {
    pub fn get_report_type(&self) -> &ReportType {
        &self.report_type
    }

    pub fn get_chr_order(&mut self) -> Vec<String> {
        self.fasta_reader.index
            .sequences().iter()
            .map(|s| s.name.clone())
            .collect_vec()
    }

    pub fn new(
        report_type: ReportType,
        mut reader: CsvReader<Box<dyn MmapBytesReader>>,
        batch_per_read: Option<usize>,
        batch_size: Option<usize>,
        fa_path: String,
        fai_path: String,
    ) -> Self {
        let fasta_reader = bio::io::fasta::IndexedReader::new(
            File::open(fa_path.clone()).expect("Failed to open fasta file"),
            File::open(fai_path.clone()).expect("Failed to open index file"),
        ).expect("Failed to read fasta file");
        info!("Opened Fasta reader from {fa_path:?} with index {fai_path:?}");

        let batch_per_read =
            batch_per_read.unwrap_or(std::thread::available_parallelism().unwrap().get());
        let batch_size = batch_size.unwrap_or(10_000);
        let chr_order = fasta_reader.index
            .sequences().iter()
            .map(|s| s.name.clone())
            .collect_vec();
        let read_queue = ReadQueue::new(batch_size, chr_order);

        let schema = SchemaRef::from(report_type.get_schema());
        let batched_reader = reader.batched_borrowed().unwrap();
        let batched_reader: BatchedCsvReader<'static> =
            unsafe { std::mem::transmute(batched_reader) };

        let batched_owned = OwnedBatchedCsvReader {
            schema,
            batched_reader,
            _reader: reader,
        };

        
        Self {
            report_type,
            batched_reader: batched_owned,
            batch_per_read,
            batch_size,
            fasta_reader,
            read_queue
        }
    }

    fn fill_queue(&mut self) -> Option<()> {
        if let Some(batches) = self.batched_reader.next_batches(self.batch_per_read).unwrap() {
            let report_type = self.report_type.clone();
            let mut join_args = JoinArgs::new(JoinType::Left);
            join_args.validation = JoinValidation::OneToOne;
            
            let chrom_batches = batches.into_iter()
                // Split all batches by chromosome
                .map(
                    |batch| UniversalBatch::from_report_type(batch, &report_type).unwrap().partition_chr()
                )
                .flatten()
                .into_grouping_map_by(|batch| batch.get_chr()).reduce(
                    |acc, _chr, new| acc.vstack(&new)
                )
                .into_iter()
                // Sort by chromosome name
                .collect_vec();
            let context_dfs = chrom_batches
                .iter()
                .map(
                    |(chr, batch)| {
                        let mut seq = Vec::new();
                        let chr_len = self.fasta_reader
                            .index.sequences().iter()
                            .find(|s| s.name == *chr)
                            .expect(format!("Sequence {chr} not found in FASTA file").as_str())
                            .len;

                        let mut start = batch.first_position();
                        if start > 2 { start -= 2; }
                        let mut end = batch.last_position();
                        if end + 2 <= chr_len as u32 { end += 2; } else { end = chr_len as u32; }

                        self.fasta_reader.fetch(chr.as_str(), start as u64, end as u64).expect(
                            format!("Failed to fetch region ({}, {}, {})", chr, start, end)
                                .as_str(),
                        );
                        self.fasta_reader
                            .read(&mut seq)
                            .expect("Failed to read fasta sequence");
                        let (positions, contexts, strands) =
                            parse_cytosines(String::from_utf8(seq.clone()).unwrap(), start.clone());
                        let pos_col = Column::new("position".into(), positions);
                        let ctx_col = Column::new("context".into(), contexts);
                        let str_col = Column::new("strand".into(), strands);
                        DataFrame::from_iter(vec![pos_col, ctx_col, str_col])
                        }
                ).collect_vec();
            let res: Vec<UniversalBatch> = 
                chrom_batches.into_par_iter()
                    .zip(context_dfs.into_par_iter())
                    .map(|((chr, batch), context_df)| {
                        let joined = context_df.lazy()
                            .join(
                                batch.get_data().clone().lazy(),
                                [col("position")], [col("position")],
                                join_args.clone()
                            )
                            .with_column(col("chr").fill_null(lit(chr.clone())));

                        let res = UniversalBatch::from(joined.collect().unwrap());
                        res
                    }).collect();
            res.into_iter().sorted_by_cached_key(|b| self.read_queue.get_chr_idx(&b.get_chr())).for_each(|r| {self.read_queue.insert(r).unwrap()});
            Some(())

        } else {
            None
        }
    }
}

impl Iterator for ReportReader {
    type Item = UniversalBatch;

    fn next(&mut self) -> Option<Self::Item> {
        if self.read_queue.is_empty() {
            match self.fill_queue() {
                _ => {},
            };
            self.next();
            debug!("filled queue.");
        }
        // Return values by chromosome order
        self.read_queue.pop()
        // End
    }
}

#[derive(Eq)]
struct ReadQueueItem {
    chr_idx: usize,
    batch: UniversalBatch,
}

impl PartialEq<Self> for ReadQueueItem {
    fn eq(&self, other: &Self) -> bool {
        (self.chr_idx == other.chr_idx) && (self.batch.first_position() == other.batch.first_position())
    }
}

impl PartialOrd<Self> for ReadQueueItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let chr_cmp = self.chr_idx.cmp(&other.chr_idx);
        Some(chr_cmp.then(other.batch.first_position().cmp(&self.batch.first_position())))
    }
}

impl Ord for ReadQueueItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl ReadQueueItem {
    fn new(chr_idx: usize, batch: UniversalBatch) -> ReadQueueItem {
        // let chr_idx = Reverse(chr_idx);
        assert!(
            batch.unique_chr(),
            "Batch has data from multiple chromosomes! {}. It must be single chromosome!", 
            batch.get_data().column("chr").unwrap().unique().unwrap().as_series().unwrap()
        );
        ReadQueueItem { chr_idx, batch }
    }
}

impl From<ReadQueueItem> for UniversalBatch {
    fn from(item: ReadQueueItem) -> Self {
        item.batch
    }
}

struct ReadQueue {
    heap: Vec<ReadQueueItem>,
    assemble_cache: Option<UniversalBatch>,
    chunk_size: usize,
    chr_order: Vec<String>,
}

impl ReadQueue {
    fn get_chr_idx(&self, chr: &String) -> usize {
        self.chr_order.iter().position(|c| c == chr).unwrap()
    }
    
    fn make_item(&self, batch: UniversalBatch) -> ReadQueueItem {
        let item = ReadQueueItem::new(self.get_chr_idx(&batch.get_chr()), batch);
        item
    }

    fn new(chunk_size: usize, chr_order: Vec<String>) -> Self {
        Self {
            heap: Vec::new(),
            assemble_cache: None,
            chr_order, chunk_size,
        }
    }

    fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }
    
    fn insert(&mut self, mut batch: UniversalBatch) -> Result<(), ()> {
        let chr = batch.get_chr();
        
        match self.assemble_cache.take() {
            None => {},
            Some(cache) => {
                if cache.get_chr() != chr {
                    self.heap.push(self.make_item(cache))
                } else {
                    batch.extend(&cache, cache.first_position() > batch.first_position())
                }
            }
        }
        
        let mut partitioned = batch.partition_by_chunks(self.chunk_size);

        if let Some(pos) = partitioned.iter().position(|b| b.get_data().height() != self.chunk_size) {
            self.assemble_cache = Some(partitioned.remove(pos));
        }
        
        partitioned.into_iter().for_each(|b| self.heap.push(self.make_item(b)));
        Ok(())
    }
    
    fn pop(&mut self) -> Option<UniversalBatch> {
        if self.heap.is_empty() {
            return None;
        }
        
        let chr_idxs = self.heap.iter().map(|b| b.chr_idx).unique();
        if chr_idxs.clone().count() == 1 {
            let min_batch = self.heap.iter().position_min_by_key(|b| b.batch.first_position()).unwrap();
            Some(self.heap.remove(min_batch).into())
        } else {
            let min_chr = chr_idxs.min().unwrap();
            let min_batch = self.heap.iter().positions(|b| b.chr_idx == min_chr).into_iter().map(|i| self.heap.get(i).unwrap()).position_min_by_key(|b| b.batch.first_position()).unwrap();
            Some(self.heap.remove(min_batch).into())
        }
    }
}
