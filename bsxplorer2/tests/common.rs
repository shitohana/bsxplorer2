#![allow(unused)]

use std::collections::BTreeMap;
use std::io::Write;
use std::ops::BitOr;

use anyhow::bail;
use bio::io::fasta::Record;
use bsxplorer2::data_structs::batch::{
    create_caregorical_dtype,
    BsxBatch,
    BsxBatchBuilder,
    BsxColumns,
};
use bsxplorer2::data_structs::ContextData;
use bsxplorer2::io::bsx::BsxFileWriter;
use bsxplorer2::io::report::ReportType;
use itertools::Itertools;
use polars::prelude::{
    AnyValue,
    Column,
    DataType,
    Scalar,
};
use polars::series::ChunkCompareEq;
use rand::{
    random,
    Rng,
    RngCore,
    SeedableRng,
};
use rand_distr::{
    Binomial,
    Distribution,
    Normal,
};

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
            Normal::new(self.mean_coverage as f64, self.std_coverage as f64).unwrap();
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
        let chr_name = self.rng.gen_range(0..1000000).to_string();
        let record = self.generate_record(self.chr_len, chr_name.clone());
        let context_data = self.generate_context_data(&record);
        let target_length = context_data.len();
        let (count_total, count_methylated) = self.generate_methylation(&context_data);

        let mut context_df = context_data.to_df();
        context_df
            .with_column(Column::new("count_total".into(), count_total))
            .unwrap();
        context_df
            .with_column(Column::new("count_m".into(), count_methylated))
            .unwrap();
        context_df
            .with_column(Column::new_scalar(
                "chr".into(),
                Scalar::new(DataType::String, AnyValue::StringOwned(chr_name.into())),
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

        let batch = BsxBatchBuilder::all_checks().cast_only(context_df).unwrap();
        Some((record, unsafe { BsxBatchBuilder::build_unchecked(batch) }))
    }
}

pub fn compare_batches(
    original: &BsxBatch,
    read: &BsxBatch,
    report_type: &ReportType,
) -> anyhow::Result<()> {
    if original.seqname() != read.seqname() {
        bail!(
            "Chromosomes differ for original: {}\nread: {}",
            original.data(),
            read.data()
        );
    }
    if original.len() != read.len() {
        bail!(
            "Heights differ for original: {}\nread: {}",
            original.data(),
            read.data()
        );
    }
    let (original_df, read_df) = {
        // Cast the "chr" column to string in both dataframes
        let mut original_df = original.data().clone();
        let mut read_df = read.data().clone();

        if let Ok(chr_col) = original_df.column("chr") {
            original_df
                .with_column(chr_col.cast(&DataType::String).unwrap())
                .unwrap();
        }

        if let Ok(chr_col) = read_df.column("chr") {
            read_df
                .with_column(chr_col.cast(&DataType::String).unwrap())
                .unwrap();
        }

        if matches!(report_type, ReportType::BedGraph) {
            original_df = original_df.drop(BsxColumns::CountM.as_str())?;
            original_df = original_df.drop(BsxColumns::CountTotal.as_str())?;
            read_df = read_df.drop(BsxColumns::CountM.as_str())?;
            read_df = read_df.drop(BsxColumns::CountTotal.as_str())?;
        }

        (original_df, read_df)
    };

    if !original_df.equals_missing(&read_df) {
        let mut colnames = original_df
            .schema()
            .iter_names_cloned()
            .filter(|v| v != BsxColumns::Density.as_str())
            .collect_vec();

        let diff_mask = colnames
            .iter()
            .map(|name| {
                (
                    original_df.column(name).unwrap().as_materialized_series(),
                    read_df.column(name).unwrap().as_materialized_series(),
                )
            })
            .map(|(orig, read)| {
                !orig.equal(read).expect(
                    format!("Equality check failed for {}, {}", orig, read).as_str(),
                )
            })
            .reduce(|acc, new| acc.bitor(new))
            .unwrap();

        if !diff_mask.any() {
            return Ok(());
        }

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
