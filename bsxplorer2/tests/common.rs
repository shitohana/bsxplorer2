#![allow(unused)]

use bio::io::fasta::Record;
use bsxplorer2::data_structs::{
    batch::{BsxBatch, BsxBatchBuilder},
    context_data::ContextData,
};
use itertools::Itertools;
use polars::prelude::{AnyValue, Column, DataType, Scalar};
use rand::{Rng, RngCore, SeedableRng};
use rand_distr::{Binomial, Distribution, Normal};

pub struct DemoReportBuilder<R: SeedableRng + RngCore> {
    chr_length: usize,
    mean_coverage: u32,
    std_coverage: f32,
    mean_methylation: f32,
    rng: R,
}

impl<R: SeedableRng + RngCore> Iterator for DemoReportBuilder<R> {
    type Item = (Record, BsxBatch);

    fn next(&mut self) -> Option<Self::Item> {
        let chr_name = self.rng.gen_range(0..1000000).to_string();
        let record = self.generate_record(self.chr_length, chr_name.clone());
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
                (context_df.column("count_m").unwrap().cast(&DataType::Float64).unwrap()
                    / context_df
                        .column("count_total")
                        .unwrap().cast(&DataType::Float64).unwrap())
                .unwrap().with_name("density".into()),
            )
            .unwrap();

        let batch = BsxBatchBuilder::all_checks()
            .build(context_df)
            .unwrap();
        Some((record, batch))
    }
}

impl<R: SeedableRng + RngCore> Default for DemoReportBuilder<R> {
    fn default() -> Self {
        Self::new(100_000, 30, 10.0, 0.4, None)
    }
}

impl<R: SeedableRng + RngCore> DemoReportBuilder<R> {
    pub(crate) fn new(
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
        let rng =
            R::seed_from_u64(seed.unwrap_or_else(|| rand::thread_rng().gen()));

        Self {
            chr_length,
            mean_coverage,
            std_coverage,
            mean_methylation,
            rng,
        }
    }

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

    /// Returns a tuple of two vectors: the first vector contains the total coverage, and the second vector contains the methylated coverage.
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
                    Binomial::new(*x as u64, self.mean_methylation as f64)
                        .unwrap();
                methylation_dist.sample(&mut self.rng) as u32
            })
            .collect_vec();

        (count_total, count_methylated)
    }

    fn set_chr_length(
        &mut self,
        chr_length: usize,
    ) {
        self.chr_length = chr_length;
    }

    fn set_mean_coverage(
        &mut self,
        mean_coverage: u32,
    ) {
        self.mean_coverage = mean_coverage;
    }

    fn set_mean_methylation(
        &mut self,
        mean_methylation: f32,
    ) {
        self.mean_methylation = mean_methylation;
    }

    fn set_std_coverage(
        &mut self,
        std_coverage: f32,
    ) {
        self.std_coverage = std_coverage;
    }
}
