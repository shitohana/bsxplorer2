//! This module provides functionalities for grouping and analyzing multiple `EncodedBsxBatch` instances.
//!
//! It includes the `EncodedBsxBatchGroup` struct, which represents a collection of `EncodedBsxBatch`
//! instances, optionally associated with labels. This allows for operations on groups of batches,
//! such as filtering, aggregating statistics, and splitting into subgroups. The module also provides
//! utility functions for checking data consistency across batches and performing various analyses.
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

/// Represents a group of encoded BSX batches.
///
/// This struct holds multiple `EncodedBsxBatch` instances, potentially associated with labels.
/// It provides methods for filtering, aggregating, and analyzing the batches as a group.
pub struct EncodedBsxBatchGroup<R: Display + Eq> {
    /// The collection of encoded BSX batches.
    batches: Vec<EncodedBsxBatch>,
    /// Optional labels associated with each batch.
    labels: Option<Vec<R>>,
}

impl<R: Display + Eq> EncodedBsxBatchGroup<R> {
    /// Returns a reference to the optional labels associated with the batches.
    pub fn labels(&self) -> &Option<Vec<R>> {
        &self.labels
    }
}
#[allow(dead_code)]
impl<R: Display + Eq + Hash + Clone + Default> EncodedBsxBatchGroup<R>
where
    R: Display + Eq + Hash + Clone + Default,
{
    /// Attempts to create a new `EncodedBsxBatchGroup`.
    ///
    /// This function performs several checks to ensure data consistency across the batches:
    /// - It verifies that the input `batches` vector is not empty.
    /// - It checks if all batches have the same height (number of rows).
    /// - It ensures that all non-empty batches have the same contig.
    /// - If labels are provided, it verifies that the number of labels matches the number of batches.
    ///
    /// # Arguments
    ///
    /// * `batches` - A vector of `EncodedBsxBatch` instances.
    /// * `labels` - An optional vector of labels, one for each batch.
    ///
    /// # Returns
    ///
    /// Returns a Result containing the new `EncodedBsxBatchGroup` if successful, or an error if any of the checks fail.
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
                .filter(|b| b.height() > 0)
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

    /// Creates a new `EncodedBsxBatchGroup` without performing any consistency checks.
    ///
    /// # Safety
    ///
    /// This function is unsafe because it bypasses the checks performed by `try_new`.
    /// It should only be used when the caller can guarantee that the input batches are consistent.
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

    /// Extracts and returns the data from each batch as a vector of `DataFrame`s.
    pub fn take_data(self) -> Vec<DataFrame> {
        self.batches.into_iter().map(DataFrame::from).collect()
    }

    /// Returns the number of samples (batches) in the group.
    #[inline]
    pub fn n_samples(&self) -> usize {
        self.batches.len()
    }

    /// Returns the height (number of rows) of the batches in the group.
    /// Assumes that all batches have the same height.
    #[inline]
    pub fn height(&self) -> usize {
        self.batches[0].height()
    }

    /// Splits the `EncodedBsxBatchGroup` into multiple groups based on the associated labels.
    ///
    /// This function iterates through the batches and their corresponding labels (if present),
    /// and groups the batches based on the labels. Batches with the same label are placed into the same group.
    ///
    /// # Returns
    ///
    /// Returns a HashMap where the keys are the unique labels and the values are `EncodedBsxBatchGroup` instances
    /// containing the batches associated with that label.
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

    /// Returns a reference to the underlying vector of `EncodedBsxBatch` instances.
    #[inline]
    pub fn batches(&self) -> &Vec<EncodedBsxBatch> {
        &self.batches
    }

    /// Calculates the average methylation density for each methylation site across all batches.
    ///
    /// # Arguments
    ///
    /// * `na_rm` - If true, missing values (NaN) are ignored in the calculation. If false, NaN values propagate.
    ///
    /// # Returns
    ///
    /// Returns a Result containing a vector of average densities (one for each methylation site), or an error if the underlying data is invalid.
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

    /// Calculates the sum of methylated counts (`count_m`) for each methylation site across all batches.
    ///
    /// # Returns
    ///
    /// Returns a Result containing a vector of summed methylated counts (one for each methylation site), or an error if the underlying data is invalid.

    pub fn get_sum_counts_m(&self) -> anyhow::Result<Vec<u32>> {
        let count_m_cols = self
            .batches
            .iter()
            .map(EncodedBsxBatch::get_counts_m)
            .collect::<anyhow::Result<Vec<_>, _>>()?;

        let (height, width) = (self.height(), self.n_samples());

        let res = (0..height)
            .into_par_iter()
            .map(|j| {
                (0..width)
                    .into_iter()
                    .map(|i| count_m_cols[i][j])
                    .map(|x| x as u32)
                    .sum()
            })
            .collect::<Vec<_>>();

        Ok(res)
    }

    /// Calculates the sum of total counts (`count_total`) for each methylation site across all batches.
    ///
    /// # Returns
    ///
    /// Returns a Result containing a vector of summed total counts (one for each methylation site), or an error if the underlying data is invalid.
    pub fn get_sum_counts_total(&self) -> anyhow::Result<Vec<u32>> {
        let count_m_cols = self
            .batches
            .iter()
            .map(EncodedBsxBatch::get_counts_total)
            .collect::<anyhow::Result<Vec<_>, _>>()?;

        let (height, width) = (self.height(), self.n_samples());

        let res = (0..height)
            .into_par_iter()
            .map(|j| {
                (0..width)
                    .into_iter()
                    .map(|i| count_m_cols[i][j])
                    .map(|x| x as u32)
                    .sum()
            })
            .collect::<Vec<_>>();

        Ok(res)
    }

    /// Retrieves the genomic positions of the methylation sites.
    ///
    /// Assumes that all batches have the same positions.
    ///
    /// # Returns
    ///
    /// Returns a Result containing a vector of genomic positions (one for each methylation site), or an error if the underlying data is invalid.
    pub fn get_positions(&self) -> anyhow::Result<Vec<u32>> {
        Ok(self.batches[0].get_position_vals()?)
    }

    /// Retrieves the chromosome (contig) name.
    ///
    /// Assumes that all batches have the same chromosome.
    ///
    /// # Returns
    ///
    /// Returns a Result containing the chromosome name, or an error if the underlying data is invalid.
    pub fn get_chr(&self) -> anyhow::Result<String> {
        Ok(self.batches[0].chr()?)
    }
}
