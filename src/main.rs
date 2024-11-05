#![feature(btree_cursors)]
#![feature(ascii_char)]

mod bsxipc;
mod genome;
mod python;
mod report;
mod sequence;
mod utils;
mod bsx_reader;
mod bsx_writer;

use bam;
use bio;
use itertools::Itertools;
use polars::io::mmap::MmapBytesReader;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{Read, Seek};
use std::ops::{Div, Sub};
use std::time::Instant;

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
    for (index, nuc) in ascii_seq.iter().enumerate() {
        let forward = match nuc {
            b'C' => true,
            b'G' => false,
            _ => { continue }
        };
        let context = if forward {
            if index >= fw_bound {continue };
            if ascii_seq[index + 1] == b'G' { Some(true) } 
            else if ascii_seq[index + 2] == b'G' { Some(false) }
            else { None }
        } else {
            if index <= rv_bound { continue };
            if ascii_seq[index - 1] == b'C' { Some(true) }
            else if ascii_seq[index - 2] == b'C' { Some(false) }
            else { None }
        };
        positions.push(start_pos + index as u64);
        contexts.push(context);
        strands.push(forward);
    };
    (positions, contexts, strands)
}

/// Stores information about batch genomic position
/// Tuple of (chromosome name, min position, max position)
struct BatchStats(String, u64, u64);

enum ReportType {
    BISMARK,
    CGMAP,
    BEDGRAPH,
    COVERAGE
}

impl ReportType {
    /// Get polars [Schema] of the report type
    fn get_schema(&self) -> Schema {
        match self {
            ReportType::BISMARK => Schema::from(PlIndexMap::from_iter([
                (PlSmallStr::from("chr"), DataType::String),
                (PlSmallStr::from("position"), DataType::UInt64),
                (PlSmallStr::from("strand"), DataType::String),
                (PlSmallStr::from("count_m"), DataType::UInt32),
                (PlSmallStr::from("count_um"), DataType::UInt32),
                (PlSmallStr::from("context"), DataType::String),
                (PlSmallStr::from("trinuc"), DataType::String),
            ])),
            ReportType::CGMAP => Schema::from(PlIndexMap::from_iter([
                (PlSmallStr::from("chr"), DataType::String),
                (PlSmallStr::from("nuc"), DataType::String),
                (PlSmallStr::from("position"), DataType::UInt64),
                (PlSmallStr::from("context"), DataType::String),
                (PlSmallStr::from("dinuc"), DataType::String),
                (PlSmallStr::from("count_m"), DataType::UInt32),
                (PlSmallStr::from("count_total"), DataType::UInt32),
            ])),
            ReportType::BEDGRAPH => Schema::from(PlIndexMap::from_iter([
                (PlSmallStr::from("chr"), DataType::String),
                (PlSmallStr::from("start"), DataType::UInt64),
                (PlSmallStr::from("end"), DataType::UInt64),
                (PlSmallStr::from("density"), DataType::Float64)
            ])),
            ReportType::COVERAGE => Schema::from(PlIndexMap::from_iter([
                (PlSmallStr::from("chr"), DataType::String),
                (PlSmallStr::from("start"), DataType::UInt64),
                (PlSmallStr::from("end"), DataType::UInt64),
                (PlSmallStr::from("density"), DataType::Float64),
                (PlSmallStr::from("count_m"), DataType::UInt32),
                (PlSmallStr::from("count_um"), DataType::UInt32)
            ]))
        }
    }

    /// Get preconfigured CSV reader for specified report type;
    ///
    /// # Defaults
    /// chunk_size = 10_000
    /// low_memory = false
    /// rechunk = true
    fn get_reader<R: MmapBytesReader>(
        &self,
        handle: R,
        chunk_size: Option<usize>,
        low_memory: Option<bool>,
        n_threads: Option<usize>,
        rechunk: Option<bool>,
    ) -> CsvReader<R> {
        let mut read_options = CsvReadOptions::default()
            // As no yet supported formats have headers
            .with_has_header(false)
            .with_chunk_size(chunk_size.unwrap_or(10_000))
            .with_low_memory(low_memory.unwrap_or(false))
            .with_n_threads(n_threads)
            .with_schema(Some(SchemaRef::from(self.get_schema())))
            .with_rechunk(rechunk.unwrap_or(true))
            .with_parse_options({
                CsvParseOptions::default()
                    .with_separator(b'\t')
                    .with_try_parse_dates(false)
                    .with_quote_char(Some(b'#'))
            });
        match self {
            ReportType::BEDGRAPH => {
                read_options = read_options.with_skip_rows(1);
            },
            _ => {}
        };
        read_options.into_reader_with_file_handle(handle)
    }

    /// Modify raw [LazyFrame] to Universal BSX schema
    ///
    /// # Returns
    ///
    /// [LazyFrame] with columns renamed/mutated to universal schema
    fn to_universal_mutate(&self, lazy_frame: LazyFrame) -> LazyFrame {
        let context_encoder = when(col("context").eq(lit("CG")))
            .then(lit(true))
            .when(col("context").eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean);

        match self {
            ReportType::BISMARK => {
                lazy_frame
                    .with_column((col("count_m") + col("count_um")).alias("count_total"))
                    .with_columns([
                        (col("count_m") / col("count_total")).alias("density"),
                        col("strand").eq(lit("+")).alias("strand"),
                        col("count_m") / col("count_total").alias("context")
                    ])
            },
            ReportType::CGMAP => {
                lazy_frame
                    .with_columns([
                        col("nuc").eq(lit("C")).alias("strand"),
                        context_encoder.alias("context")
                    ])
            },
            ReportType::BEDGRAPH => {
                lazy_frame
                    .rename(["start"], ["position"], true)
                    .drop(["end"])
                    .with_columns([
                        lit(NULL).alias("count_m"),
                        lit(NULL).alias("count_total"),
                        col("density").div(lit(100)).alias("density"),
                ])
            },
            ReportType::COVERAGE => {
                lazy_frame
                    .rename(["start"], ["position"], true)
                    .drop(["end"])
            },
        }

    }
}

/// Iterator over methylation report data
struct ReportReader<'a> {
    /// [ReportType] of the report
    report_type: ReportType,
    /// Initialized [BatchedCsvReader]. Use [ReportType::get_reader] to create binding for
    /// reader and [CsvReader::batched_borrowed] to create Batched reader.
    batched_reader: BatchedCsvReader<'a>,
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

impl<'a> ReportReader<'a> {
    fn new(
        report_type: ReportType,
        mut batched_csv_reader: BatchedCsvReader<'a>,
        batch_per_read: Option<usize>,
        batch_size: Option<usize>,
        fa_path: Option<&str>,
        fai_path: Option<&str>,
    ) -> Self {
        let fasta_reader = match (fa_path, fai_path) {
            (Some(fa_path), Some(fai_path)) => {
                Some(
                    bio::io::fasta::IndexedReader::new(
                        File::open(fa_path).expect("Failed to open fasta file"),
                        File::open(fai_path).expect("Failed to open index file"),
                    ).expect("Failed to read fasta file")
                )
            },
            (None, None) => None,
            _ => panic!("Both fasta and fasta index path need to be specified"),
        };

        match (&report_type, &fasta_reader) {
            (ReportType::BEDGRAPH, None) => { panic!("FASTA file is required for .bedGraph file reading") },
            (ReportType::COVERAGE, None) => { panic!("FASTA file is required for .cov file reading") },
            _ => {}
        }

        let batch_per_read = batch_per_read.unwrap_or(
            std::thread::available_parallelism().unwrap().get()
        );
        let batch_queue = VecDeque::new();
        let batch_size = batch_size.unwrap_or(10_000);
        Self {
            report_type,
            batched_reader: batched_csv_reader,
            batch_per_read,
            batch_size,
            fasta_reader,
            batch_queue
        }
    }

    fn fill_queue(&mut self) -> Option<()> {
        let read_res = self.batched_reader.next_batches(
            self.batch_per_read
        ).expect("Failed to read batches");
        if read_res.is_none() { return None; };

        let mut batches = read_res.unwrap().iter()
            // Partition by chromosomes to allow direct sequence reading
            .map(|batch| {
                batch.partition_by(["chr"], true).unwrap()
            })
            .flatten()
            .map(|batch| {
                self.report_type
                    .to_universal_mutate(batch.lazy())
                    .collect().unwrap()
            })
            .collect::<Vec<_>>();

        batches = if self.fasta_reader.is_some() {
            let mut fasta_reader = self.fasta_reader.as_mut().unwrap();
            // Collect stats about batches to prepare fetch args
            let batch_stats: Vec<BatchStats> = batches.iter()
                .map(|batch| -> BatchStats {
                    let positions = batch.column("position").unwrap().u64().unwrap();
                    let chr = batch.column("chr").unwrap().str().unwrap().first().unwrap();
                    let min = positions.first().unwrap();
                    let max = positions.last().unwrap();
                    BatchStats(String::from(chr), min, max)
                })
                .collect();

            // Read sequences, as bio::io::fasta does not implement synchronized reading
            let sequences: Vec<(Vec<u8>, u64)> = batch_stats.iter().map(|stats| {
                let mut seq = Vec::new();
                fasta_reader.fetch(stats.0.as_str(), stats.1 - 2, stats.2 + 2)
                    .expect(
                        format!("Failed to fetch region ({}, {}, {})", stats.0, stats.1 - 2, stats.2 + 2).as_str(),
                    );
                fasta_reader.read(&mut seq).expect("Failed to read fasta sequence");
                (seq, stats.1)
            }).collect();

            // Process sequences and extract cytosine contexts.
            // Convert them into a DataFrame
            let context_df_cols = vec![String::from("context"), String::from("strand")];
            let context_dfs = sequences.par_iter().map(|(seq, start)| {
                let (positions, contexts, strands) = parse_cytosines(
                    String::from_utf8(seq.clone()).unwrap(), start.clone()
                );
                let pos_col = Column::new("position".into(), positions);
                let ctx_col = Column::new("context".into(), contexts);
                let str_col = Column::new("strand".into(), strands);
                DataFrame::from_iter(vec![pos_col, ctx_col, str_col])
            }).collect::<Vec<_>>();
            let keep_cols = batches.get(0).unwrap()
                .schema()
                .iter_names().filter_map(|name| {
                if !context_df_cols.contains(&name.to_string().into()) {
                    Some(name.clone())
                } else { None }
            }).collect::<Vec<_>>();

            // Join raw data to all possible cytosine contexts from reference file.
            let joined: Vec<LazyFrame> = itertools::izip!(batches, context_dfs).map(|(raw, contexts)| {
                let chr = raw.column("chr").unwrap().str().unwrap().first().unwrap();
                let trimmed = raw.select(keep_cols.clone()).unwrap();
                let joined = contexts.lazy().join(
                    trimmed.lazy(),
                    [col("position")],
                    [col("position")],
                    JoinArgs::new(JoinType::Left)
                )
                    .with_columns([
                        col("chr").fill_null(lit(chr)),
                    ]);
                joined
            }).collect();

            // Concatenate LazyFrames before splitting into equal-size batches
            let full = concat(joined, UnionArgs::default()).unwrap()
                .collect()
                .expect("Could not concat joined data");

            // Split into equal-size batches
            full
                .partition_by(["chr"], true)
                .unwrap()
                .iter().map(
                    |batch| {
                        batch
                            .with_row_index("batch".into(), None).unwrap()
                            .lazy()
                            .with_column(col("batch").div(lit(self.batch_size as u32)).sub(lit(0.5)).cast(DataType::UInt32))
                            .collect().unwrap()
                            .partition_by(["batch"], false).unwrap()
                    }
                )
                .flatten()
                .collect()
        } else {batches};

        for batch in {
            batches.iter()
                .sorted_by_cached_key(|batch| {
                    batch.column("position").unwrap().u64().unwrap().first().unwrap()
                }).collect::<Vec<_>>()
        } {
            self.batch_queue.push_back(batch.to_owned());
        }
        Some(())
    }
}

impl<'a> Iterator for ReportReader<'a> {
    type Item = DataFrame;

    fn next(&mut self) -> Option<Self::Item> {
        if self.batch_queue.is_empty() {
            self.fill_queue();
        }
        self.batch_queue.pop_front()
    }
}


fn main() {
    // Some params
    let chunk_size: usize = 100_000;

    let start = Instant::now();
    
    // Initialize paths
    let fasta_path: String = "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa".into();
    let index_path: String = "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa.fai".into();
    let bedgraph_path: String = "/Users/shitohana/Documents/CX_reports/A_thaliana.bedGraph".into();
    let ipc_path: String = "/Users/shitohana/Desktop/RustProjects/bsxplorer2/bedgraph.ipc".into();

    // Initialize FASTA reader
    // let mut fasta_reader = bio::io::fasta::IndexedReader::new(
    //     File::open(fasta_path.clone()).expect("Failed to open fasta file"),
    //     File::open(index_path.clone()).expect("Failed to open index file"),
    // ).expect("Failed to read fasta file");

    let report_type = ReportType::BEDGRAPH;
    // Initialize CSV reader
    let mut reader = report_type.get_reader(
        File::open(bedgraph_path).expect("Failed to open fasta file"),
        Some(chunk_size),
        Some(false),
        None,
        Some(true)
    );
    let batched = reader.batched_borrowed()
        .expect("reader could not be batched");
    
    // // Initialize writer
    // let binding = fasta_reader.index
    //     .sequences();
    // let chr_names = binding.iter()
    //     .map(|x| {x.name.as_str()})
    //     .collect::<Vec<_>>();
    // println!("chr_names {:?}", chr_names);
    // let universal_schema = get_universal_schema(chr_names);
    // let universal_names = universal_schema.iter_names().map(|name| name.as_str().into()).collect::<Vec<String>>();
    // let mut writer = IpcWriter::new(File::create(ipc_path).unwrap())
    //     .batched(&universal_schema).expect("Failed to create IpcWriter");

    let iterator = ReportReader::new(
        report_type, batched, None, None, Some(fasta_path.as_str()), Some(index_path.as_str())
    );

    for batch in iterator {
        println!("batch {:?}", batch);
    }

    let elapsed = start.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
}