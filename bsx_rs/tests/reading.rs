#![feature(path_add_extension)]
extern crate pretty_env_logger;
#[macro_use] extern crate log;

use rayon::iter::ParallelIterator;
use rayon::iter::IndexedParallelIterator;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::ops::Div;
use std::os;
use std::path::PathBuf;
use hashbrown::HashSet;
use itertools::{Itertools};
use num::{PrimInt};
use polars::prelude::*;
use rand::distributions::Distribution;
use rayon::prelude::*;
use bsx_rs::io::report::bsx_batch::BsxBatch;
use bsx_rs::io::report::fasta_reader::FastaReader;
use bsx_rs::io::report::reader::{ContextData, ReportReader, ReportReaderBuilder};
use bsx_rs::io::report::report_batch_utils::{decode_context, decode_strand};
use bsx_rs::io::report::schema::ReportTypeSchema;
use bsx_rs::region::{GenomicPosition};

const MEAN_COVERAGE: f64 = 40.0;
const STD_COVERAGE: f64 = 30.0;
const MAX_CHROMS: usize = 10;
const MAX_CHR_LEN: usize = 1_000_000;


#[derive(Debug, Clone, Eq, PartialEq)]
struct ReportMetadata {
    cytosines_per_chr: HashMap<String, HashSet<u64>>,
}


impl ReportMetadata {
    fn new() -> Self {
        Self {
            cytosines_per_chr: HashMap::new(),
        }
    }
    
    fn get_mut(&mut self, name: &str) -> Option<&mut HashSet<u64>> {
        self.cytosines_per_chr.get_mut(name)
    }

    fn append_chr(&mut self, name: String, cytosines: &[u64]) {
        self.cytosines_per_chr.insert(name, HashSet::from_iter(cytosines.iter().cloned()));
    }
    
    fn extend_positions(&mut self, name: String, cytosines: &[u64]) {
        self.cytosines_per_chr.get_mut(&name).unwrap().extend(cytosines);
    }

    fn cytoines_total(&self) -> usize {
        self.cytosines_per_chr.values().map(HashSet::len).sum()
    }
}

fn create_sample_report<W>(
    report_type: ReportTypeSchema,
    handle: W,
    fasta_sink: W
) -> ReportMetadata
where
    W: Write, 
{
    let context_prob: HashMap<Option<bool>, f64> = HashMap::from_iter([
        (Some(true), 0.4),   // CG
        (Some(false), 0.15), // CHG
        (None, 0.05),        // CHH
    ]);
    /*
    Steps of modeling cytosine report

    1. Generate number of chromosomes. THey will be named as "chr%d".
    2. Generate lengths of chromosomes, length < 2^24.
    3. Generate DNA sequence as uniform distribution of
       nucleotides
    4. Extract context data from the sequence.
    5. Model counts data as Binomial distribution. Total count of
       reads is distributed normally.
    6. Create BsxBatch, which then can be converted to any report
       type
    */

    let chr_number = rand::random::<u8>() % MAX_CHROMS as u8 + 1;

    let chr_lengths = (0..chr_number)
        .map(|_| rand::random::<u32>() % MAX_CHR_LEN as u32 + 1).collect_vec();

    info!("Number of chroms: {}", chr_number);
    
    let mut metadata = ReportMetadata::new();

    let mut writer = CsvWriter::new(handle)
        .with_separator(b'\t')
        .include_header(matches!(report_type, ReportTypeSchema::BedGraph))
        .include_bom(false)
        .with_float_scientific(Some(false))
        .batched(&report_type.schema())
        .unwrap();
    let mut fasta_writer = bio::io::fasta::Writer::new(fasta_sink);

    for (chr_i, chr_length) in chr_lengths.iter().enumerate() {
        let chr_name = format!("chr{}", chr_i);

        let sequence = generate_sequence(*chr_length as usize);
        info!("Writing sequence {} with length {}", chr_name, sequence.len());
        fasta_writer.write(chr_name.as_str(), None, &sequence).unwrap();
        
        let context_data = ContextData::from_sequence(
            &sequence, 
            GenomicPosition::new(chr_name.clone(), 1)
        );
        assert!(context_data.positions().is_sorted());
        let context_len = context_data.len();
        metadata.append_chr(chr_name.to_string(), &context_data.positions().iter().map(|v| *v as u64).collect_vec());
        
        let counts_total_col = generate_counts_total::<u32>(MEAN_COVERAGE, STD_COVERAGE, context_len);
        let counts_m_col = generate_counts_m::<u32>(&counts_total_col, context_data.contexts(), &context_prob);

        let result_bsx = context_data_to_bsx(
            context_data, counts_m_col, counts_total_col, chr_name.as_str()
        ).unwrap();
        let mut result_casted = report_type.report_mutate_from_bsx(result_bsx).unwrap().select(report_type.col_names().iter().cloned()).unwrap();
        result_casted.align_chunks_par();
        writer.write_batch(&result_casted).unwrap()
    }
    writer.finish().expect("Failed to finish writing CSV file");
    fasta_writer.flush().expect("Failed to flush fasta file");
    info!("Report written");
    metadata
}

fn generate_sequence(chr_length: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();

    let nuc_mapping: HashMap<i32, u8> = HashMap::from_iter([
        (0, b'A'),
        (1, b'C'),
        (2, b'G'),
        (3, b'T'),
    ]);
    
    rand::distributions::Uniform::new_inclusive(0, 3)
        .sample_iter(&mut rng).take(chr_length)
        .map(|n| nuc_mapping[&n])
        .collect_vec()
}


fn generate_counts_total<N>(mean: f64, std: f64, length: usize) -> Vec<N>
where N: PrimInt {
    let counts_dist = rand_distr::Normal::new(mean, std).unwrap();
    let mut rng = rand::thread_rng();

    counts_dist.sample_iter(&mut rng)
        .take(length)
        .map(|v| N::from(v.abs().floor()).unwrap())
        .collect_vec()
}

fn context_data_to_bsx<N>(context_data: ContextData, counts_m: Vec<N>, counts_total: Vec<N>, chr_name: &str) -> PolarsResult<DataFrame> 
where N: PrimInt {
    let mut result_df = context_data.into_dataframe()?;
    result_df.with_column(Series::new("count_m".into(), counts_m.iter().map(|v| v.to_u64().unwrap()).collect_vec()))?;
    result_df.with_column(Series::new("count_total".into(), counts_total.iter().map(|v| v.to_u64().unwrap()).collect_vec()))?;

    let mut result_lazy = result_df.lazy();
    result_lazy = decode_context(result_lazy, "context", "context");
    result_lazy = decode_strand(result_lazy, "strand", "strand");

    result_lazy = result_lazy
        .with_columns([
            lit(chr_name).alias("chr"),
            (col("count_m").cast(DataType::Float64).div(col("count_total"))).alias("density"),
        ])
        .cast(BsxBatch::hashmap(), true);
    result_lazy.collect()
}

fn generate_counts_m<N>(counts_total: &[N], contexts: &[Option<bool>], encoded_context_prob: &HashMap<Option<bool>, f64>) -> Vec<N>
where 
    N: PrimInt + Sync + Send,
{
    counts_total.par_iter().zip(contexts.par_iter())
        .map(|(count_total, context)| {
            let context_prob = encoded_context_prob[context];
            let dist = rand_distr::Binomial::new(count_total.clone().to_u64().unwrap(), context_prob).unwrap();
            N::from(dist.sample(&mut rand::thread_rng())).unwrap()
        }).collect()
}

fn test_reading(reader: ReportReader, metadata: ReportMetadata, chunk_size: usize) {
    let mut read_metadata = ReportMetadata::new();
    metadata.cytosines_per_chr.keys().for_each(|k|
        read_metadata.append_chr(k.clone(), &[])
    );
    let mut shortened_batch_num = 0;
    for batch in reader {
        let batch_position = batch.first_position().unwrap();
        let batch_height = batch.data().height();
        // Check chr unique
        // ----------------
        assert_eq!(batch.data().column("chr").unwrap().as_series().unwrap().unique().unwrap().len(), 1);
        // Check position sorted
        // ---------------------
        assert!(
            batch.data()
                .column("position").unwrap()
                .as_series().unwrap()
                .is_sorted(SortOptions::default().with_order_descending(false)).unwrap(),
            "Got unsorted batch {}",
            { 
                let mut prev = 0;
                let mut idx = 0;
                for (i, val) in batch.data().column("position").unwrap().u64().unwrap().into_iter().enumerate() {
                    if let Some(v) = val {
                        if v >= prev {
                            prev = v;
                        } else {
                            idx = i;
                            break
                        }
                    }
                }
                
                batch.data().slice((idx - 2) as i64, 5)
            }
        );
        // Check position unique
        // ---------------------
        assert_eq!(batch.data().column("position").unwrap().as_series().unwrap().unique().unwrap().len(), batch_height);
        
        if batch.data().height() != chunk_size {
            shortened_batch_num += 1;
        }
        
        read_metadata.extend_positions(
            batch_position.chr().to_string(), 
            &batch.data().column("position").unwrap().u64().unwrap().iter().map(|v| v.unwrap()).collect_vec()
        );
    }
    
    assert_eq!(shortened_batch_num, read_metadata.cytosines_per_chr.len());
    
    for chr in metadata.cytosines_per_chr.keys() {
        let reference = metadata.cytosines_per_chr.get(chr).unwrap();
        let real = read_metadata.cytosines_per_chr.get(chr).unwrap();
        let difference = reference.difference(real).cloned().sorted().collect_vec();
        
        if !difference.is_empty() {
            println!("{:?}", reference.iter().sorted().zip(real.iter().sorted()).take(20).collect::<Vec<_>>());;
            panic!("Context positions differ! {:?}", difference);
        }
    }
}

struct SampleSetup {
    report_file: tempfile::NamedTempFile,
    fasta_file: tempfile::NamedTempFile,
    metadata: ReportMetadata,
}

fn test_report_type(report_type_schema: ReportTypeSchema, chunk_size: usize)  {
    let tempfile = tempfile::NamedTempFile::new().unwrap();
    let tempfile_fasta = tempfile::NamedTempFile::new().unwrap();
    let tempfile_fai = tempfile::NamedTempFile::new().unwrap();
    
    info!("Report saved as {:?}", tempfile.path());
    info!("FASTA saved as {:?}", tempfile_fasta.path());
    let sink = tempfile.reopen().unwrap();
    let fasta_sink = tempfile_fasta.reopen().unwrap();
    let metadata = create_sample_report(report_type_schema, sink, fasta_sink);

    {
        rust_htslib::faidx::build(tempfile_fasta.path()).expect("Failed to build .fai");
        let index_path = tempfile_fasta.path().with_added_extension("fai");
        info!("Build index {:?}", index_path);
        
        let mut fai_handle = File::open(index_path.clone()).unwrap();
        let mut index: Vec<u8> = Vec::new();
        fai_handle.read_to_end(&mut index).unwrap();
        std::fs::remove_file(index_path.clone()).unwrap();
        info!("Removed index {:?}", index_path);
        let mut fai_sink = tempfile_fai.reopen().unwrap();
        fai_sink.write_all(&index).unwrap();
    }
    info!("Fasta index created as: {:?}", tempfile_fai.path());

    let reader = {
        let mut builder = ReportReaderBuilder::default()
            .with_chunk_size(chunk_size)
            .with_report_type(report_type_schema)
            .with_batch_size(2 << 30);
        
        if report_type_schema.need_align() {
            builder = builder
                .with_fasta(PathBuf::from(tempfile_fasta.path()), PathBuf::from(tempfile_fai.path()));
        }
        
        builder.try_finish(tempfile.into_file()).unwrap()
    };

    test_reading(reader, metadata, chunk_size);
}

#[test]
fn report_reading() {
    /*
    This is done:
        1.  Fix BSX to Bedgraph/Coverage conversion.
            Currently, NaN columns with no reads are not
            dropped
        2.  Find out, why on alignment step the last batch
            of the chromosome is not aligned with the last
            nucleotides
        3.  Find out, why aligned batch can be unsorted by
            position
        4.  Error with reports less than chunk size
        
    To remember: there were several problems with 
    reports alignments:
            1.  bio::io::seq reader is zero-based, so the
                start position must be 0
            2.  Extra number of nucleotides per sequence
                part is now set to 3 to ensure any possible
                context position will be captured
            3.  Many small bugs/mistakes, which I already
                don't remember.  
        
    TODO new:  
        1.  Rewrite what's written so it can be
            understood and be expandable. Add documentation
            and unit tests.
        2.  Rewrite testing strategy, so different ReportTypes
            could operate the same sequence information
        3.  Add seeding to tests to make them reproducible
        4.  Create bulk test with many possible parameters
            for ReportReader
        5.  Create should_panic tests for both FastaReader and
            ReportReader
    */
    
    pretty_env_logger::init();
    let chunk_size = 10_000;
    let report_schemas = [
        ReportTypeSchema::BedGraph,
        ReportTypeSchema::Bismark,
        ReportTypeSchema::Coverage,
        ReportTypeSchema::CgMap,
    ];
    for report_type in report_schemas {
        info!("Testing {:?}", report_type);
        test_report_type(report_type, chunk_size);
    }
    info!("All done!");
}

