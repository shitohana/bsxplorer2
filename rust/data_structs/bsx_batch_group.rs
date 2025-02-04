use crate::data_structs::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use crate::utils::types::{Context, IPCEncodedEnum, Strand};
use anyhow::anyhow;
use itertools::{izip, Itertools};
use log::warn;
use polars::datatypes::BooleanChunked;
use polars::error::PolarsError;
use polars::frame::DataFrame;
use polars::prelude::NamedFrom;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use statrs::statistics::Statistics;
use std::collections::HashMap;
use std::fmt::Display;
use std::hash::Hash;

macro_rules! check_eq_batches {
    ($expression: expr, $err_string: expr) => {
        'eval: {
            let mut evaled = $expression;
            if evaled.all_equal() {
                break 'eval Ok(());
            } else {
                break 'eval Err(anyhow!($err_string).context(format!("{:?}", evaled)));
            }
        }
    };
}

pub struct EncodedBsxBatchGroup<R: Display + Eq> {
    batches: Vec<EncodedBsxBatch>,
    labels: Option<Vec<R>>,
}

impl<R: Display + Eq> EncodedBsxBatchGroup<R> {
    pub fn labels(&self) -> &Option<Vec<R>> {
        &self.labels
    }
}

impl<R: Display + Eq + Hash + Clone + Default> EncodedBsxBatchGroup<R> {
    pub fn try_new(batches: Vec<EncodedBsxBatch>, labels: Option<Vec<R>>) -> anyhow::Result<Self> {
        if batches.is_empty() {
            return Err(anyhow!("Empty batches group"));
        }
        check_eq_batches!(
            batches.iter().map(EncodedBsxBatch::height),
            "Batches lengths differ"
        )?;
        check_eq_batches!(
            batches
                .iter()
                .map(|b| EncodedBsxBatch::as_contig(b, &mut None))
                .map(|c| c.unwrap()),
            "Batches differ"
        )?;
        if labels
            .as_ref()
            .map(|l| l.len() != batches.len())
            .unwrap_or(false)
        {
            return Err(anyhow!("Batches and labels lengths differ"));
        }

        Ok(EncodedBsxBatchGroup { batches, labels })
    }

    unsafe fn new_unchecked(batches: Vec<EncodedBsxBatch>) -> Self {
        Self {
            batches,
            labels: None,
        }
    }

    /// Check full positions, strands and contexts consistency between batches
    fn check_expensive(&self) -> bool {
        check_eq_batches!(
            self.batches
                .iter()
                .map(|b| { b.data().select(["position", "strand", "context"]).unwrap() }),
            "Batches data differ"
        )
        .is_ok()
    }

    /// Apply [BooleanChunked] mask to the batches data
    fn filter_mask(self, mask: &BooleanChunked) -> anyhow::Result<Self> {
        let filtered_batches: anyhow::Result<Vec<EncodedBsxBatch>, PolarsError> = self
            .batches
            .into_par_iter()
            .map(|batch| batch.filter_mask(&mask))
            .collect();
        Ok(Self {
            batches: filtered_batches?,
            labels: self.labels,
        })
    }

    /// Filter batches data by [Context]. Filter mask is calculated once and
    /// applied to all batches
    pub fn filter_context(self, context: Context) -> anyhow::Result<Self> {
        let context_bool = context.to_bool();
        let df0 = &self.batches[0];

        let mask = {
            let context_bool = context_bool;
            let mask_iter = df0
                .data()
                .column("context")?
                .bool()?
                .iter()
                .map(|v| v == context_bool)
                .collect_vec();
            BooleanChunked::new("mask".into(), mask_iter)
        };

        self.filter_mask(&mask)
    }

    /// Filter batches data by [Strand]. Filter mask is calculated once and
    /// applied to all batches
    pub fn filter_strand(self, strand: Strand) -> anyhow::Result<Self> {
        let strand_bool = strand.to_bool();
        let df0 = &self.batches[0];

        let mask = {
            let mask_iter = df0
                .data()
                .column("strand")?
                .bool()?
                .iter()
                .map(|v| v == strand_bool)
                .collect_vec();
            BooleanChunked::new("mask".into(), mask_iter)
        };

        self.filter_mask(&mask)
    }

    /// Replaces low coverage methylation sites with count_total = 0, density = NaN
    pub fn mark_low_counts(self, threshold: i16) -> anyhow::Result<Self> {
        let marked: anyhow::Result<Vec<EncodedBsxBatch>, PolarsError> = self
            .batches
            .into_par_iter()
            .map(|batch| batch.mark_low_counts(threshold))
            .collect();
        Ok(Self {
            batches: marked?,
            labels: self.labels,
        })
    }

    /// Filter rows where number of missing (NaN) methylation site densities <= `n_missing`
    pub fn filter_n_missing(self, n_missing: usize) -> anyhow::Result<Self> {
        let mut threshold = n_missing;
        if threshold > self.n_samples() {
            warn!(
                "Value for filter_n_missing {}, is bigger, than sample count {}. Setting it to {}",
                threshold,
                self.n_samples(),
                self.n_samples()
            );
            threshold = self.n_samples();
        }

        let mask = {
            let nan_vecs = self
                .batches
                .par_iter()
                .map(|batch| -> anyhow::Result<Vec<bool>, PolarsError> {
                    Ok(batch
                        .data()
                        .column("density")?
                        .f32()?
                        .iter()
                        .map(|v| v.map(|x| x.is_nan()).unwrap_or(true))
                        .collect())
                })
                .collect::<anyhow::Result<Vec<Vec<bool>>, PolarsError>>()?;
            let (width, height) = (self.n_samples(), self.height());

            let mask_vec: Vec<bool> = (0..height)
                .into_par_iter()
                .map(|j| {
                    let mut nan_count = 0usize;
                    for i in 0..width {
                        if nan_vecs[i][j] {
                            nan_count += 1;
                        }
                    }
                    nan_count <= threshold
                })
                .collect();
            BooleanChunked::new("mask".into(), mask_vec)
        };

        self.filter_mask(&mask)
    }

    pub fn take_data(self) -> Vec<DataFrame> {
        self.batches.into_iter().map(DataFrame::from).collect()
    }

    #[inline]
    pub fn n_samples(&self) -> usize {
        self.batches.len()
    }

    #[inline]
    pub fn height(&self) -> usize {
        self.batches[0].height()
    }

    pub fn split_groups(self) -> HashMap<R, EncodedBsxBatchGroup<R>> {
        let mut out = HashMap::<R, Vec<EncodedBsxBatch>>::new();
        let batches_len = self.batches.len();
        for (batch, label) in izip!(
            self.batches,
            self.labels.unwrap_or(vec![R::default(); batches_len])
        ) {
            out.entry(label).or_insert_with(Vec::new).push(batch);
        }
        HashMap::from_iter(
            out.into_iter()
                .map(|(k, v)| (k, unsafe { EncodedBsxBatchGroup::new_unchecked(v) })),
        )
    }

    #[inline]
    pub fn batches(&self) -> &Vec<EncodedBsxBatch> {
        &self.batches
    }

    pub fn get_average_density(&self, na_rm: bool) -> anyhow::Result<Vec<f64>> {
        let density_cols = self
            .batches
            .iter()
            .map(EncodedBsxBatch::get_density_vals)
            .collect::<anyhow::Result<Vec<_>, _>>()?;

        let (height, width) = (self.height(), self.n_samples());

        let res = (0..height)
            .into_par_iter()
            .map(|j| {
                (0..width)
                    .into_iter()
                    .map(|i| density_cols[i][j])
                    .filter(|&x| if na_rm { !x.is_nan() } else { true })
                    .map(|x| x as f64)
                    .mean()
            })
            .collect::<Vec<_>>();

        Ok(res)
    }

    pub fn get_positions(&self) -> anyhow::Result<Vec<u32>> {
        Ok(self.batches[0].get_position_vals()?)
    }

    pub fn get_chr(&self) -> anyhow::Result<String> {
        Ok(self.batches[0].chr()?)
    }
}
