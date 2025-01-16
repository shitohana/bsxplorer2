#![feature(path_add_extension)]
extern crate pretty_env_logger;
#[macro_use]
extern crate log;

use bsx_rs::io::report::bsx_batch::BsxBatch;
use bsx_rs::io::report::reader::{ContextData, ReportReader, ReportReaderBuilder};
use bsx_rs::io::report::report_batch_utils::{decode_context, decode_strand};
use bsx_rs::io::report::schema::ReportTypeSchema;
use bsx_rs::region::GenomicPosition;
use bsx_rs::utils::types::{Context, IPCEncodedEnum};
use hashbrown::HashSet;
use itertools::Itertools;
use num::PrimInt;
use polars::prelude::*;
use rand::distributions::uniform::SampleUniform;
use rand::distributions::Distribution;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::ops::{Div, Range};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use tempfile::NamedTempFile;

#[derive(Debug)]
struct TestReportConfig {
    mean_coverage: f64,
    std_coverage: f64,
    chroms_num: Range<usize>,
    chr_len: Range<usize>,
    context_prob: Vec<(Context, f64)>,
    nuc_prob: [f64; 4],
    seed: Option<u64>,
    chunk_size: usize,
}

impl Default for TestReportConfig {
    fn default() -> Self {
        Self {
            mean_coverage: 40f64,
            std_coverage: 30f64,
            chroms_num: 5..10,
            chr_len: 5e5 as usize..3e6 as usize,
            context_prob: vec![
                (Context::CG, 0.4),
                (Context::CHG, 0.15),
                (Context::CHH, 0.05),
            ],
            nuc_prob: [0.25; 4],
            seed: None,
            chunk_size: 10_000,
        }
    }
}

type ReferenceMetadataEntry = DataFrame;

struct ReferenceMetadata {
    data: HashMap<String, ReferenceMetadataEntry>,
    order: Vec<String>,
}

impl ReferenceMetadata {
    pub fn new() -> Self {
        Self {
            data: HashMap::new(),
            order: Vec::new(),
        }
    }

    pub fn insert(&mut self, key: String, value: ReferenceMetadataEntry) {
        self.data.insert(key.clone(), value);
        self.order.push(key);
    }

    pub fn get(&self, key: &str) -> Option<&ReferenceMetadataEntry> {
        self.data.get(key)
    }

    pub fn get_mut(&mut self, key: &str) -> Option<&mut ReferenceMetadataEntry> {
        self.data.get_mut(key)
    }

    pub fn order(&self) -> &Vec<String> {
        &self.order
    }
}

struct TestReportEnv {
    fasta_tempfile: NamedTempFile,
    fai_tempfile: NamedTempFile,
    config: TestReportConfig,
    rng: ChaCha8Rng,
    reference_metadata: ReferenceMetadata,
}

impl TestReportEnv {
    fn new(config: TestReportConfig) -> Self {
        let fasta_tempfile = NamedTempFile::new().unwrap();
        let fai_tempfile = NamedTempFile::new().unwrap();
        let mut reference_metadata = ReferenceMetadata::new();

        let mut rng = if let Some(seed) = config.seed {
            ChaCha8Rng::seed_from_u64(seed)
        } else {
            ChaCha8Rng::from_entropy()
        };

        let fasta_sink = fasta_tempfile.reopen().unwrap();
        let mut fasta_writer = bio::io::fasta::Writer::new(fasta_sink);

        let chr_number = Self::generate_in_range(&mut rng, config.chroms_num.clone());
        let chr_length = (0..chr_number)
            .map(|_| Self::generate_in_range(&mut rng, config.chr_len.clone()))
            .collect_vec();

        for (idx, chr_length) in (0..chr_number).zip(chr_length) {
            let chr_name = format!("chr{}", idx);

            info!("Chrom '{}': {}", chr_name, chr_length);

            let chr_sequence = Self::generate_sequence(&mut rng, chr_length, config.nuc_prob);

            fasta_writer
                .write(chr_name.as_str(), None, &chr_sequence)
                .unwrap();

            let context_data = ContextData::from_sequence(
                &chr_sequence,
                GenomicPosition::new(chr_name.clone(), 1),
            );
            let context_len = context_data.len();
            assert!(context_data.positions().is_sorted());

            let counts_total_col = Self::generate_counts_total::<u32>(
                &mut rng,
                config.mean_coverage,
                config.std_coverage,
                context_len,
            );
            let counts_m_col = Self::generate_counts_m::<u32>(
                &mut rng,
                &counts_total_col,
                context_data.contexts(),
                &config.context_prob,
            );

            let bsx_result = context_data_to_bsx(
                context_data,
                counts_m_col,
                counts_total_col,
                chr_name.as_str(),
            )
            .unwrap();

            reference_metadata.insert(chr_name, bsx_result);
        }

        fasta_writer.flush().unwrap();

        let fai_sink = fai_tempfile.reopen().unwrap();
        Self::create_fai(fasta_tempfile.path(), fai_sink).expect("Could not create FAI");

        info!("Initialized ReportTestingEnv with config: {:?}", config);
        info!("Fasta saved as: {}", fasta_tempfile.path().display());
        info!("Fasta index saved as: {}", fai_tempfile.path().display());

        Self {
            fasta_tempfile,
            fai_tempfile,
            reference_metadata,
            rng,
            config,
        }
    }

    fn test_report_type(&self, report_type: ReportTypeSchema) -> Result<(), Box<dyn Error>> {
        // Write report
        info!("Testing report type {:?}", report_type);
        let report_tempfile = NamedTempFile::new()?;
        self.write_report(report_tempfile.reopen()?, report_type)?;

        let mut read_metadata = CytosinesHashmap::new();
        self.reference_metadata
            .order
            .iter()
            .for_each(|k| read_metadata.append_chr(k.clone(), &[]));

        let reader = {
            let mut builder = ReportReaderBuilder::default()
                .with_chunk_size(self.config.chunk_size)
                .with_report_type(report_type)
                .with_batch_size(2 << 30);

            if report_type.need_align() {
                builder = builder.with_fasta(
                    PathBuf::from(self.fasta_tempfile.path()),
                    PathBuf::from(self.fai_tempfile.path()),
                );
            }

            builder.try_finish(report_tempfile.reopen()?)?
        };

        let mut shortened_batch_num = 0;

        for batch in reader {
            Self::check_chr_unique(&batch);
            Self::check_positions_sorted(&batch);
            Self::check_positions_unique(&batch);

            if batch.data().height() != self.config.chunk_size {
                shortened_batch_num += 1;
            }

            read_metadata.extend_positions(
                batch.first_position()?.chr().to_string(),
                &batch
                    .data()
                    .column("position")?
                    .u64()?
                    .iter()
                    .map(|v| v.unwrap())
                    .collect_vec(),
            );
        }

        assert_eq!(shortened_batch_num, read_metadata.c_positions.len());

        for chr in self.reference_metadata.order() {
            let reference = HashSet::from_iter(
                self.reference_metadata
                    .get(chr)
                    .unwrap()
                    .column("position")?
                    .u64()?
                    .into_iter()
                    .map(|v| v.unwrap()),
            );
            let real = read_metadata.c_positions.get(chr).unwrap();
            let difference = reference.difference(real).cloned().sorted().collect_vec();

            if !difference.is_empty() {
                println!(
                    "{:?}",
                    reference
                        .iter()
                        .sorted()
                        .zip(real.iter().sorted())
                        .take(20)
                        .collect::<Vec<_>>()
                );
                panic!("Context positions differ! {:?}", difference);
            }
        }
        info!("Finished testing report type {:?}. SUCCESS", report_type);
        report_tempfile.close()?;
        Ok(())
    }

    fn check_chr_unique(batch: &BsxBatch) {
        let batch_chroms = batch
            .data()
            .column("chr")
            .unwrap()
            .as_series()
            .unwrap()
            .unique()
            .unwrap();
        assert_eq!(batch_chroms.len(), 1, "Chr should be unique");
    }

    fn check_positions_sorted(batch: &BsxBatch) {
        assert!(
            batch
                .data()
                .column("position")
                .unwrap()
                .as_series()
                .unwrap()
                .is_sorted(SortOptions::default().with_order_descending(false))
                .unwrap(),
            "Got unsorted batch {}",
            {
                let mut prev = 0;
                let mut idx = 0;
                for (i, val) in batch
                    .data()
                    .column("position")
                    .unwrap()
                    .u64()
                    .unwrap()
                    .into_iter()
                    .enumerate()
                {
                    if let Some(v) = val {
                        if v >= prev {
                            prev = v;
                        } else {
                            idx = i;
                            break;
                        }
                    }
                }

                batch.data().slice((idx - 2) as i64, 5)
            }
        );
    }

    fn check_positions_unique(batch: &BsxBatch) {
        assert_eq!(
            batch
                .data()
                .column("position")
                .unwrap()
                .as_series()
                .unwrap()
                .unique()
                .unwrap()
                .len(),
            batch.data().height(),
            "Got non-unique positions"
        );
    }

    fn write_report<W>(&self, sink: W, report_type: ReportTypeSchema) -> Result<(), Box<dyn Error>>
    where
        W: Write,
    {
        let mut writer = CsvWriter::new(sink)
            .with_separator(b'\t')
            .include_header(matches!(report_type, ReportTypeSchema::BedGraph))
            .include_bom(false)
            .with_float_scientific(Some(false))
            .batched(&report_type.schema())
            .unwrap();

        for chr_name in self.reference_metadata.order() {
            let bsx_data = self.reference_metadata.get(chr_name).unwrap();

            let mut report_data = report_type.report_mutate_from_bsx(bsx_data.clone())?;
            report_data.align_chunks_par();

            writer.write_batch(&report_data)?;
        }

        writer.finish()?;
        Ok(())
    }

    fn create_fai<W>(fasta_path: &Path, mut fai_sink: W) -> io::Result<()>
    where
        W: Write,
    {
        rust_htslib::faidx::build(fasta_path).expect("Failed to build .fai");

        let index_path = fasta_path.with_added_extension("fai");
        let mut fai_handle = File::open(index_path.clone())?;
        let mut index: Vec<u8> = Vec::new();
        fai_handle.read_to_end(&mut index)?;

        std::fs::remove_file(index_path.clone())?;

        fai_sink.write_all(&index)?;
        Ok(())
    }

    fn generate_in_range<N>(rng: &mut ChaCha8Rng, range: Range<N>) -> N
    where
        N: SampleUniform + PartialOrd,
    {
        rng.gen_range(range)
    }

    /// Nucleotide probabilities - in the same order as A, T, G, C
    fn generate_sequence(
        rng: &mut ChaCha8Rng,
        length: usize,
        nuc_probabilities: [f64; 4],
    ) -> Vec<u8> {
        let prob = nuc_probabilities;
        let a = prob[0];
        let t = prob[1] + a;
        let g = prob[2] + t;
        let c = prob[3] + g;
        assert_eq!(c, 1f64);

        rand_distr::Uniform::new_inclusive(0f64, 1f64)
            .sample_iter(rng)
            .take(length)
            .map(|n| {
                if n <= a {
                    b'A'
                } else if n <= t {
                    b'T'
                } else if n <= g {
                    b'G'
                } else {
                    b'C'
                }
            })
            .collect_vec()
    }

    fn generate_counts_total<N>(rng: &mut ChaCha8Rng, mean: f64, std: f64, length: usize) -> Vec<N>
    where
        N: SampleUniform + PrimInt,
    {
        let counts_dist = rand_distr::Normal::new(mean, std).unwrap();
        counts_dist
            .sample_iter(rng)
            .take(length)
            .map(|v| N::from(v.abs().floor()).unwrap())
            .collect_vec()
    }

    fn generate_counts_m<N>(
        rng: &mut ChaCha8Rng,
        counts_total: &[N],
        encoded_contexts: &[Option<bool>],
        contexts_prob: &[(Context, f64)],
    ) -> Vec<N>
    where
        N: PrimInt + Sync + Send,
    {
        let rng_mutex = Mutex::new(rng);
        let encoded_context_prob: HashMap<Option<bool>, f64> = HashMap::from_iter(
            contexts_prob
                .iter()
                .map(|(context, prob)| (context.to_bool(), *prob)),
        );
        counts_total
            .par_iter()
            .zip(encoded_contexts.par_iter())
            .map(|(count_total, context)| {
                let context_prob = encoded_context_prob[context];
                let dist =
                    rand_distr::Binomial::new(count_total.clone().to_u64().unwrap(), context_prob)
                        .unwrap();
                N::from(dist.sample(*rng_mutex.lock().unwrap())).unwrap()
            })
            .collect()
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct CytosinesHashmap {
    c_positions: HashMap<String, HashSet<u64>>,
}

impl CytosinesHashmap {
    fn new() -> Self {
        Self {
            c_positions: HashMap::new(),
        }
    }

    fn get_mut(&mut self, name: &str) -> Option<&mut HashSet<u64>> {
        self.c_positions.get_mut(name)
    }

    fn append_chr(&mut self, name: String, cytosines: &[u64]) {
        self.c_positions
            .insert(name, HashSet::from_iter(cytosines.iter().cloned()));
    }

    fn extend_positions(&mut self, name: String, cytosines: &[u64]) {
        self.c_positions.get_mut(&name).unwrap().extend(cytosines);
    }

    fn cytoines_total(&self) -> usize {
        self.c_positions.values().map(HashSet::len).sum()
    }
}

fn context_data_to_bsx<N>(
    context_data: ContextData,
    counts_m: Vec<N>,
    counts_total: Vec<N>,
    chr_name: &str,
) -> PolarsResult<DataFrame>
where
    N: PrimInt,
{
    let mut result_df = context_data.into_dataframe()?;
    result_df.with_column(Series::new(
        "count_m".into(),
        counts_m.iter().map(|v| v.to_u64().unwrap()).collect_vec(),
    ))?;
    result_df.with_column(Series::new(
        "count_total".into(),
        counts_total
            .iter()
            .map(|v| v.to_u64().unwrap())
            .collect_vec(),
    ))?;

    let mut result_lazy = result_df.lazy();
    result_lazy = decode_context(result_lazy, "context", "context");
    result_lazy = decode_strand(result_lazy, "strand", "strand");

    result_lazy = result_lazy
        .with_columns([
            lit(chr_name).alias("chr"),
            (col("count_m")
                .cast(DataType::Float64)
                .div(col("count_total")))
            .alias("density"),
        ])
        .cast(BsxBatch::hashmap(), true);
    result_lazy.collect()
}

#[test]
fn report_reading() {
    pretty_env_logger::init();
    // let chunk_size = 10_000;
    let report_schemas = [
        ReportTypeSchema::BedGraph,
        ReportTypeSchema::Bismark,
        ReportTypeSchema::Coverage,
        ReportTypeSchema::CgMap,
    ];
    let config = TestReportConfig {
        seed: Some(1234),
        ..Default::default()
    };
    let env = TestReportEnv::new(config);

    for report_type in report_schemas {
        match env.test_report_type(report_type) {
            Ok(()) => {}
            Err(err) => {
                panic!("{}", err)
            }
        }
    }
}
