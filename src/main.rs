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

use std::fs::File;
use std::ops::{Div, Sub};
use std::time::Instant;
use polars::prelude::*;
use rayon::prelude::*;
use polars_arrow::io::ipc::read::{FileReader, read_file_metadata};
use bio;
use bam;
use itertools::Itertools;
use polars::frame::row::Row;
use polars::prelude::Expr::Column;
use crate::report::get_universal_schema;

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

struct BatchStats(String, u64, u64);


fn main() {
    // Some params
    let chunk_size: usize = 100_000;
    let batch_per_read: usize = 16;
    let rbatch_size: usize = 10_000;
    
    let start = Instant::now();
    
    // Initialize paths
    let fasta_path: String = "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa".into();
    let index_path: String = "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa.fai".into();
    let bedgraph_path: String = "/Users/shitohana/Documents/CX_reports/A_thaliana.bedGraph".into();
    let ipc_path: String = "/Users/shitohana/Desktop/RustProjects/bsxplorer2/bedgraph.ipc".into();
    
    // Initialize FASTA reader
    let mut fasta_reader = bio::io::fasta::IndexedReader::new(
        File::open(&fasta_path).expect("Failed to open fasta file"),
        File::open(index_path).expect("Failed to open index file"),
    ).expect("Failed to read fasta file");
    
    // Initialize CSV reader
    let schema = Schema::from(PlIndexMap::from_iter([
        (PlSmallStr::from("chr"), DataType::String),
        (PlSmallStr::from("start"), DataType::UInt64),
        (PlSmallStr::from("end"), DataType::UInt64),
        (PlSmallStr::from("density"), DataType::Float64)
    ]));
    let mut reader = {
        let parse_options = CsvParseOptions::default()
            .with_separator(b'\t')
            .with_comment_prefix(Some("#"));

        CsvReadOptions::default()
            .with_has_header(false)
            .with_chunk_size(chunk_size)
            .with_rechunk(true)
            .with_parse_options(parse_options)
            .with_skip_rows(1)
            .with_schema(Some(SchemaRef::from(schema.clone())))
            .try_into_reader_with_file_path(Some(bedgraph_path.into()))
            .expect("Failed to create reader")
    };
    let mut batched = reader.batched_borrowed()
        .expect("reader could not be batched");
    
    // Initialize writer
    let binding = fasta_reader.index
        .sequences();
    let chr_names = binding.iter()
        .map(|x| {x.name.as_str()})
        .collect::<Vec<_>>();
    println!("chr_names {:?}", chr_names);
    let universal_schema = get_universal_schema(chr_names);
    let universal_names = universal_schema.iter_names().map(|name| name.as_str().into()).collect::<Vec<String>>();
    let mut writer = IpcWriter::new(File::create(ipc_path).unwrap())
        .batched(&universal_schema).expect("Failed to create IpcWriter");

    // Fill cols
    let fill_exprs = [
        lit(NULL).alias("count_m"),
        lit(NULL).alias("count_total"),
    ];

    // Iterate over CSV
    while let Ok(Some(batches)) = batched.next_batches(batch_per_read) {
        // Partition batches by chromosome to allow direct sequence fetching
        let batches = batches.iter()
            .map(|batch| {batch.partition_by(["chr"], true).unwrap()})
            .flatten()
            .collect::<Vec<_>>();
        
        // Collect stats about batches to prepare fetch args
        let batch_stats: Vec<BatchStats> = batches.iter()
            .map(|batch| -> BatchStats {
                let positions = batch.column("start").unwrap().u64().unwrap();
                let chr = batch.column("chr").unwrap().str().unwrap().first().unwrap();
                let min = positions.first().unwrap();
                let max = positions.last().unwrap();
                BatchStats(String::from(chr), min, max)
            })
            .collect();
        
        // Read sequences, as bio::io::fasta does not implement synchronized reading
        let sequences: Vec<(Vec<u8>, u64)> = batch_stats.iter().map(|stats| {
            let mut seq = Vec::new();
            fasta_reader.fetch(stats.0.as_str(), stats.1 - 2, stats.2 + 2).expect(
                format!("Failed to fetch region ({}, {}, {})", stats.0, stats.1 - 2, stats.2 + 2).as_str(),
            );
            fasta_reader.read(&mut seq).expect("Failed to read fasta sequence");
            (seq, stats.1)
        }).collect();
        
        // Process sequences and extract cytosine contexts.
        // Convert them into a DataFrame
        let context_dfs = sequences.par_iter().map(|(seq, start)| {
            let (positions, contexts, strands) = parse_cytosines(
                String::from_utf8(seq.clone()).unwrap(), start.clone()
            );
            let pos_col = polars::frame::column::Column::new("position".into(), positions);
            let ctx_col = polars::frame::column::Column::new("contexts".into(), contexts);
            let str_col = polars::frame::column::Column::new("strand".into(), strands);
            DataFrame::from_iter(vec![pos_col, ctx_col, str_col])
        }).collect::<Vec<_>>();
        
        // Join raw data to all possible cytosine contexts from reference file.
        let joined: Vec<LazyFrame> = itertools::izip!(batches, context_dfs).map(|(raw, contexts)| {
            let chr = raw.column("chr").unwrap().str().unwrap().first().unwrap();
            let trimmed = raw.select(["start", "density"]).unwrap();
            let joined = contexts.lazy().join(
                trimmed.lazy(),
                [col("position")],
                [col("start")],
                JoinArgs::new(JoinType::Left)
            )
                .with_columns([
                    lit(chr).alias("chr"),
                    col("density").fill_null(lit(0)).div(lit(100))
                ]);
            joined
        }).collect();
        // Concatenate LazyFrames before splitting into equal-size batches
        let full = concat(joined, UnionArgs::default()).unwrap()
            .with_columns(fill_exprs.clone())
            .collect()
            .expect("Could not concat joined data");
        println!("full {:?}", full);
        // Split into equal-size batches
        full
            .partition_by(["chr"], true)
            .unwrap()
            .iter().map(
            |batch| {
                batch
                    .with_row_index("batch".into(), None).unwrap()
                    .lazy()
                    .with_column(col("batch").div(lit(rbatch_size as u32)).sub(lit(0.5)).cast(DataType::UInt32))
                    .collect().unwrap()
                    .partition_by(["batch"], false).unwrap()
            }
            )
            .flatten()
            .for_each(|mut batch| { 
                let batch = batch.align_chunks().to_owned();
                
            });
    }
    let elapsed = start.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
}