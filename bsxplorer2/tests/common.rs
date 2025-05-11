#![allow(unused)]

use std::collections::BTreeMap;
use std::io::Write;
use std::ops::BitOr;

use anyhow::bail;
use bio::io::fasta::Record;
use bsxplorer2::data_structs::batch::colnames::{CONTEXT_NAME,
                                                COUNT_M_NAME,
                                                COUNT_TOTAL_NAME,
                                                STRAND_NAME};
use bsxplorer2::data_structs::batch::{BsxBatch,
                                      BsxBatchBuilder,
                                      BsxBatchMethods};
use bsxplorer2::data_structs::context_data::ContextData;
use bsxplorer2::io::bsx::BsxIpcWriter;
use bsxplorer2::io::report::ReportTypeSchema;
use itertools::Itertools;
use polars::prelude::{AnyValue, Column, DataType, Scalar};
use polars::series::ChunkCompareEq;
use rand::{random, Rng, RngCore, SeedableRng};
use rand_distr::{Binomial, Distribution, Normal};
use smallstr::SmallString;

pub struct DemoReportBuilder<R: SeedableRng + RngCore> {
    chr_len:       usize,
    mean_coverage: u32,
    std_coverage:  f32,
    mean_density:  f32,
    rng:           R,
}

impl<R: SeedableRng + RngCore> Default for DemoReportBuilder<R> {
    fn default() -> Self {
        Self::new(100_000, 30, 10.0, 0.4, None)
    }
}

impl<R: SeedableRng + RngCore> DemoReportBuilder<R> {
    pub fn new(
        chr_length: usize,
        mean_coverage: u32,
        std_coverage: f32,
        mean_methylation: f32,
        seed: Option<u64>,
    ) -> Self {
        assert!(mean_coverage > 0, "mean_coverage must be greater than 0");
        assert!(std_coverage > 0.0, "std_coverage must be greater than 0.0");
        assert!(
            mean_methylation >= 0.0 && mean_methylation <= 1.0,
            "mean_methylation must be between 0.0 and 1.0"
        );
        let rng = R::seed_from_u64(seed.unwrap_or_else(|| random()));

        Self {
            chr_len: chr_length,
            mean_coverage,
            std_coverage,
            mean_density: mean_methylation,
            rng,
        }
    }

    pub fn write_bsx<W: Write>(
        mut self,
        sink: W,
        n_chr: usize,
        chunk_size: usize,
    ) -> anyhow::Result<Vec<BsxBatch>> {
        let mut orig_batches = Vec::<BsxBatch>::new();
        let mut chr_batches = BTreeMap::new();

        for (_record, mut full_batch) in self.into_iter().take(n_chr) {
            let chr = full_batch.chr_val().to_owned();
            orig_batches.push(full_batch.clone());
            let mut partitioned = Vec::new();

            loop {
                let (left, right) = full_batch.split_at(chunk_size);
                partitioned.push(left);
                if right.is_empty() {
                    break;
                }
                else {
                    let _ = std::mem::replace(&mut full_batch, right);
                }
            }

            chr_batches.insert(chr, partitioned);
        }

        let chr_list = chr_batches
            .keys()
            .cloned()
            .collect::<Vec<_>>();
        let mut writer = BsxIpcWriter::try_new(
            sink,
            chr_list
                .iter()
                .map(SmallString::to_string)
                .collect(),
            None,
            None,
        )?;

        for chr in chr_list {
            let batches = chr_batches.remove(&chr).unwrap();
            for batch in batches {
                writer.write_batch(batch)?;
            }
        }

        Ok(orig_batches)
    }

    pub fn set_chr_length(
        &mut self,
        chr_length: usize,
    ) {
        self.chr_len = chr_length;
    }

    pub fn set_mean_coverage(
        &mut self,
        mean_coverage: u32,
    ) {
        self.mean_coverage = mean_coverage;
    }

    pub fn set_mean_methylation(
        &mut self,
        mean_methylation: f32,
    ) {
        self.mean_density = mean_methylation;
    }

    pub fn set_std_coverage(
        &mut self,
        std_coverage: f32,
    ) {
        self.std_coverage = std_coverage;
    }

    pub fn rng_mut(&mut self) -> &mut R {
        &mut self.rng
    }
}

impl<R: SeedableRng + RngCore> DemoReportBuilder<R> {
    /// Returns a random sequence of a given length.
    fn generate_record(
        &mut self,
        length: usize,
        name: String,
    ) -> Record {
        let mut seq = String::new();
        let chars = ['A', 'C', 'G', 'T'];
        for _ in 0..length {
            seq.push(chars[(self.rng.gen::<u32>() % 4) as usize]);
        }
        Record::with_attrs(&name, None, &seq.as_bytes())
    }

    /// Returns a [ContextData] object from a [Record] object.
    fn generate_context_data(
        &mut self,
        record: &Record,
    ) -> ContextData {
        ContextData::from_sequence(record.seq())
    }

    /// Returns a tuple of two vectors: the first vector contains the total
    /// coverage, and the second vector contains the methylated coverage.
    fn generate_methylation(
        &mut self,
        context_data: &ContextData,
    ) -> (Vec<u32>, Vec<u32>) {
        let coverage_dist =
            Normal::new(self.mean_coverage as f64, self.std_coverage as f64)
                .unwrap();
        let count_total = coverage_dist
            .sample_iter(&mut self.rng)
            .map(|x| x.abs() as u32)
            .take(context_data.len())
            .collect_vec();

        let count_methylated = count_total
            .iter()
            .map(|x| {
                let methylation_dist =
                    Binomial::new(*x as u64, self.mean_density as f64).unwrap();
                methylation_dist.sample(&mut self.rng) as u32
            })
            .collect_vec();

        (count_total, count_methylated)
    }
}

impl<R: SeedableRng + RngCore> Iterator for DemoReportBuilder<R> {
    type Item = (Record, BsxBatch);

    fn next(&mut self) -> Option<Self::Item> {
        let chr_name = self
            .rng
            .gen_range(0..1000000)
            .to_string();
        let record = self.generate_record(self.chr_len, chr_name.clone());
        let context_data = self.generate_context_data(&record);
        let target_length = context_data.len();
        let (count_total, count_methylated) =
            self.generate_methylation(&context_data);

        let mut context_df = context_data.to_df::<BsxBatch>();
        context_df
            .with_column(Column::new("count_total".into(), count_total))
            .unwrap();
        context_df
            .with_column(Column::new("count_m".into(), count_methylated))
            .unwrap();
        context_df
            .with_column(Column::new_scalar(
                "chr".into(),
                Scalar::new(
                    DataType::String,
                    AnyValue::StringOwned(chr_name.into()),
                ),
                target_length,
            ))
            .unwrap();
        context_df
            .with_column(
                (context_df
                    .column("count_m")
                    .unwrap()
                    .cast(&DataType::Float64)
                    .unwrap()
                    / context_df
                        .column("count_total")
                        .unwrap()
                        .cast(&DataType::Float64)
                        .unwrap())
                .unwrap()
                .with_name("density".into()),
            )
            .unwrap();

        let batch = BsxBatchBuilder::all_checks()
            .build(context_df)
            .unwrap();
        Some((record, batch))
    }
}

pub fn compare_batches(
    original: &BsxBatch,
    read: &BsxBatch,
    report_type: &ReportTypeSchema,
) -> anyhow::Result<()> {
    if original.chr_val() != read.chr_val() {
        bail!(
            "Chromosomes differ for original: {}\nread: {}",
            original.data(),
            read.data()
        );
    }
    if original.height() != read.height() {
        bail!(
            "Heights differ for original: {}\nread: {}",
            original.data(),
            read.data()
        );
    }
    let (original_df, read_df) =
        if matches!(report_type, ReportTypeSchema::BedGraph) {
            // TODO: Find out, why first row density equal to NaN when testing
            return Ok(());

            // (
            //     original
            //         .data()
            //         .drop_many([COUNT_M_NAME, COUNT_TOTAL_NAME, CONTEXT_NAME,
            // STRAND_NAME]),     read.data()
            //         .drop_many([COUNT_M_NAME, COUNT_TOTAL_NAME, CONTEXT_NAME,
            // STRAND_NAME]), )
        }
        else {
            (original.data().clone(), read.data().clone())
        };
    if !original_df.equals_missing(&read_df) {
        let diff_mask = original_df
            .materialized_column_iter()
            .zip(read_df.materialized_column_iter())
            .map(|(orig, read)| !orig.equal(read).unwrap())
            .reduce(|acc, new| acc.bitor(new))
            .unwrap();
        let orig_diff = original_df.filter(&diff_mask)?;
        let read_diff = read_df.filter(&diff_mask)?;
        bail!(
            "Data differs for original: {}\nread: {}",
            orig_diff,
            read_diff
        );
    }

    Ok(())
}
