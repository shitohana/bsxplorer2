use crate::io::report::types::ReportType;
use itertools::Itertools;
use polars::datatypes::DataType;
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::io::mmap::MmapBytesReader;
use polars::prelude::*;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::collections::VecDeque;
use std::fs::File;
use std::ops::{Add, Div, Sub};

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
fn parse_cytosines(seq: String, start_pos: u64) -> (Vec<u64>, Vec<Option<bool>>, Vec<bool>) {
    let start_pos = start_pos - 1;
    let fw_bound: usize = seq.len() - 2;
    let rv_bound: usize = 2;

    let mut positions: Vec<u64> = Vec::new();
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
            }
            else if ascii_seq[index + 2] == b'G' {
                Some(false)
            }
            else {
                None
            }
        }
        else {
            if index <= rv_bound {
                continue 'seq_iter;
            };
            if ascii_seq[index - 1] == b'C' {
                Some(true)
            }
            else if ascii_seq[index - 2] == b'C' {
                Some(false)
            }
            else {
                None
            }
        };
        positions.push(start_pos + index as u64);
        contexts.push(context);
        strands.push(forward);
    }
    (positions, contexts, strands)
}

/// Stores information about batch genomic position
/// Tuple of (chromosome name, min position, max position)
pub struct BatchStats(String, u64, u64);

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
    /// Number of rows in the output batches. All batches will be same specified
    /// length if possible.
    batch_size: usize,
    /// If [ReportType] is [ReportType::COVERAGE] or [ReportType::BEDGRAPH], this option
    /// should be specified with initialized FASTA reader. If other report types -
    /// it is optional to apply additional alignment to reference cytosines.
    fasta_reader: Option<bio::io::fasta::IndexedReader<File>>,
    /// Queue of computed batches.
    batch_queue: VecDeque<DataFrame>,
}

impl ReportReader {
    pub fn get_report_type(&self) -> &ReportType {
        &self.report_type
    }

    pub fn new(
        report_type: ReportType,
        mut reader: CsvReader<Box<dyn MmapBytesReader>>,
        batch_per_read: Option<usize>,
        batch_size: Option<usize>,
        fa_path: Option<String>,
        fai_path: Option<String>,
    ) -> Self {
        let fasta_reader = match (fa_path, fai_path) {
            (Some(fa_path), Some(fai_path)) => Some(
                bio::io::fasta::IndexedReader::new(
                    File::open(fa_path).expect("Failed to open fasta file"),
                    File::open(fai_path).expect("Failed to open index file"),
                )
                .expect("Failed to read fasta file"),
            ),
            (None, None) => None,
            _ => panic!("Both fasta and fasta index path need to be specified"),
        };

        match (&report_type, &fasta_reader) {
            (ReportType::BEDGRAPH, None) => {
                panic!("FASTA file is required for .bedGraph file reading")
            }
            (ReportType::COVERAGE, None) => {
                panic!("FASTA file is required for .cov file reading")
            }
            _ => {}
        }

        let batch_per_read = batch_per_read.unwrap_or(
            std::thread::available_parallelism()
                .unwrap()
                .get(),
        );
        let batch_queue = VecDeque::new();
        let batch_size = batch_size.unwrap_or(10_000);

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
            batch_queue,
        }
    }

    fn fill_queue(&mut self) -> Option<()> {
        let read_res = self
            .batched_reader
            .next_batches(self.batch_per_read)
            .expect("Failed to read batches");
        if read_res.is_none() {
            return None;
        };
        let mut batches = read_res
            .unwrap()
            .iter()
            // Partition by chromosomes to allow direct sequence reading
            .map(|batch| {
                batch
                    .partition_by(["chr"], true)
                    .unwrap()
            })
            .flatten()
            .map(|batch| {
                self.report_type
                    .to_universal_mutate(batch.lazy())
                    .collect()
                    .unwrap()
            })
            .collect::<Vec<_>>();

        batches = if self.fasta_reader.is_some() {
            let fasta_reader = self.fasta_reader.as_mut().unwrap();
            // Collect stats about batches to prepare fetch args
            let batch_stats: Vec<BatchStats> = batches
                .iter()
                .map(|batch| -> BatchStats {
                    let positions = batch
                        .column("position")
                        .unwrap()
                        .u64()
                        .unwrap();
                    let chr = batch
                        .column("chr")
                        .unwrap()
                        .str()
                        .unwrap()
                        .first()
                        .unwrap();
                    let min = positions.first().unwrap();
                    let max = positions.last().unwrap();
                    BatchStats(String::from(chr), min, max)
                })
                .collect();

            // Read sequences, as bio::io::fasta does not implement synchronized reading
            let sequences: Vec<(Vec<u8>, u64)> = batch_stats
                .iter()
                .map(|stats| {
                    let mut seq = Vec::new();
                    {
                        let chr = stats.0.as_str();
                        let chr_len = fasta_reader
                            .index
                            .sequences()
                            .iter()
                            .find(|s| s.name == chr)
                            .expect(format!("Sequence {chr} not found in FASTA file").as_str())
                            .len;
                        let start = if stats.1 >= 2 { stats.1 - 2 } else { stats.1 };
                        let stop = if stats.2 + 2 <= chr_len {
                            stats.2 + 2
                        }
                        else {
                            chr_len
                        };
                        fasta_reader
                            .fetch(chr, start, stop)
                            .expect(
                                format!("Failed to fetch region ({}, {}, {})", chr, start, stop)
                                    .as_str(),
                            );
                        fasta_reader
                            .read(&mut seq)
                            .expect("Failed to read fasta sequence");
                    }
                    (seq, stats.1)
                })
                .collect();

            // Process sequences and extract cytosine contexts.
            // Convert them into a DataFrame
            let context_df_cols = vec![String::from("context"), String::from("strand")];
            let context_dfs = sequences
                .par_iter()
                .map(|(seq, start)| {
                    let (positions, contexts, strands) =
                        parse_cytosines(String::from_utf8(seq.clone()).unwrap(), start.clone());
                    let pos_col = Column::new("position".into(), positions);
                    let ctx_col = Column::new("context".into(), contexts);
                    let str_col = Column::new("strand".into(), strands);
                    DataFrame::from_iter(vec![pos_col, ctx_col, str_col])
                })
                .collect::<Vec<_>>();
            let keep_cols = batches
                .get(0)
                .unwrap()
                .schema()
                .iter_names()
                .filter_map(|name| {
                    if !context_df_cols.contains(&name.to_string().into()) {
                        Some(name.clone())
                    }
                    else {
                        None
                    }
                })
                .collect::<Vec<_>>();

            // Join raw data to all possible cytosine contexts from reference file.
            let joined: Vec<LazyFrame> = itertools::izip!(batches, context_dfs)
                .map(|(raw, contexts)| {
                    let chr = raw
                        .column("chr")
                        .unwrap()
                        .str()
                        .unwrap()
                        .first()
                        .unwrap();
                    let trimmed = raw.select(keep_cols.clone()).unwrap();
                    let joined = contexts
                        .lazy()
                        .join(
                            trimmed.lazy(),
                            [col("position")],
                            [col("position")],
                            JoinArgs::new(JoinType::Left),
                        )
                        .with_columns([col("chr").fill_null(lit(chr))]);
                    joined
                })
                .collect();

            // Concatenate LazyFrames before splitting into equal-size batches
            let full = concat(joined, UnionArgs::default())
                .unwrap()
                .collect()
                .expect("Could not concat joined data");

            // Split into equal-size batches
            full.partition_by(["chr"], true)
                .unwrap()
        }
        else {
            batches
        };

        for batch in {
            batches
                .iter()
                .sorted_by_cached_key(|batch| {
                    batch
                        .column("position")
                        .unwrap()
                        .u64()
                        .unwrap()
                        .first()
                        .unwrap()
                })
                .collect::<Vec<_>>()
        } {
            self.batch_queue
                .push_back(batch.to_owned());
        }
        Some(())
    }

    pub fn get_chr_names(&mut self) -> Option<Vec<String>> {
        if self.fasta_reader.is_some() {
            let fasta_reader = self.fasta_reader.as_mut().unwrap();
            let names = fasta_reader
                .index
                .sequences()
                .iter()
                .map(|seq| seq.name.clone())
                .collect::<Vec<_>>();
            Some(names)
        }
        else {
            None
        }
    }
}

impl Iterator for ReportReader {
    type Item = DataFrame;

    fn next(&mut self) -> Option<Self::Item> {
        if self.batch_queue.is_empty() {
            self.fill_queue();
        }
        self.batch_queue.pop_front()
    }
}
