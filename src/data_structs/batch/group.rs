/*******************************************************************************
 Copyright (c) 2025
 The Prosperity Public License 3.0.0

 Contributor: [shitohana](https://github.com/shitohana)

 Source Code: https://github.com/shitohana/BSXplorer
 ******************************************************************************/

//! This module provides functionalities for grouping and analyzing multiple
//! `EncodedBsxBatch` instances.
//!
//! It includes the `EncodedBsxBatchGroup` struct, which represents a collection
//! of `EncodedBsxBatch` instances, optionally associated with labels. This
//! allows for operations on groups of batches, such as filtering, aggregating
//! statistics, and splitting into subgroups. The module also provides
//! utility functions for checking data_structs consistency across batches and
//! performing various analyses.
//!
//! The module supports operations like:
//! - Filtering batches by context or strand
//! - Marking low-coverage methylation sites
//! - Calculating average methylation densities
//! - Summing counts across batches
//! - Splitting batches into labeled groups
use std::collections::HashMap;
use std::fmt::Display;
use std::hash::Hash;

use anyhow::{anyhow, Context as AnyhowContext};
use itertools::{izip, Itertools};
use log::{debug, info, warn};
use ndarray::{Array2, Axis};
use polars::datatypes::BooleanChunked;
use polars::error::PolarsError;
use polars::frame::DataFrame;
use polars::prelude::{ChunkedArray, NamedFrom, PolarsNumericType};
use rayon::iter::{IntoParallelIterator,
                  ParallelIterator};

use crate::data_structs::batch::encoded::EncodedBsxBatch;
use crate::data_structs::batch::traits::BsxBatchMethods;
use crate::utils::types::{Context, IPCEncodedEnum, Strand};

/// Macro to check if all items in a collection are equal
///
/// This macro evaluates an expression that produces a collection, checks if all
/// elements are equal, and returns an appropriate Result.
///
/// # Arguments
///
/// * `$expression` - Expression that produces a collection to check
/// * `$err_string` - Error message to return if elements are not equal
///
/// # Returns
///
/// * `Ok(())` if all elements are equal
/// * `Err(anyhow::Error)` with context if elements are not equal
macro_rules! check_eq_batches {
    ($expression: expr, $err_string: expr) => {
        'eval: {
            let mut evaled = $expression;
            if evaled.all_equal() {
                break 'eval Ok(());
            }
            else {
                break 'eval Err(anyhow!($err_string)
                    .context(format!("Values: {:?}", evaled)));
            }
        }
    };
}

pub(crate) fn apply_over_column<P, R, T, F>(
    columns: Vec<&ChunkedArray<P>>, 
    func: F
) -> anyhow::Result<T>
where
    P: PolarsNumericType<Native = R>,
    F: FnOnce(Array2<Option<R>>) -> T
{
    assert!(columns.len() > 0);
    
    let matrix: Array2<Option<R>> = Array2::from_shape_vec(
        (columns.len(), columns[0].len()),
        columns.into_iter().map(|a| a.to_vec()).concat()
    )?;
    
    Ok(func(matrix))
}

/// Represents a group of encoded BSX batches with optional labels.
///
/// This struct holds multiple `EncodedBsxBatch` instances, potentially
/// associated with labels. It provides methods for filtering, aggregating, and
/// analyzing the batches as a group.
///
/// All batches in a group are expected to have:
/// - The same number of rows (height)
/// - The same chromosome (contig)
/// - Identical genomic positions, strands, and contexts
pub struct EncodedBsxBatchGroup<R: Display + Eq> {
    /// The collection of encoded BSX batches, each representing a sample
    batches: Vec<EncodedBsxBatch>,
    /// Optional labels associated with each batch (e.g., condition, tissue
    /// type, etc.)
    labels:  Option<Vec<R>>,
}

impl<R: Display + Eq> EncodedBsxBatchGroup<R> {
    /// Returns a reference to the optional labels associated with the batches.
    ///
    /// # Returns
    ///
    /// A reference to the optional vector of labels, one for each batch.
    pub fn labels(&self) -> &Option<Vec<R>> { &self.labels }
}

#[allow(dead_code)]
impl<R: Display + Eq + Hash + Clone + Default> EncodedBsxBatchGroup<R>
where
    R: Display + Eq + Hash + Clone + Default,
{
    /// Attempts to create a new `EncodedBsxBatchGroup` with validation.
    ///
    /// This function performs several checks to ensure data_structs consistency
    /// across the batches:
    /// - It verifies that the input `batches` vector is not empty.
    /// - It checks if all batches have the same height (number of rows).
    /// - It ensures that all non-empty batches have the same contig.
    /// - If labels are provided, it verifies that the number of labels matches
    ///   the number of batches.
    ///
    /// # Arguments
    ///
    /// * `batches` - A vector of `EncodedBsxBatch` instances representing
    ///   samples.
    /// * `labels` - An optional vector of labels, one for each batch.
    ///
    /// # Returns
    ///
    /// Returns a Result containing the new `EncodedBsxBatchGroup` if
    /// successful, or an error if any of the checks fail.
    pub fn try_new(
        batches: Vec<EncodedBsxBatch>,
        labels: Option<Vec<R>>,
    ) -> anyhow::Result<Self> {
        if batches.is_empty() {
            return Err(anyhow!(
                "Cannot create batch group: Empty batches provided"
            ));
        }

        check_eq_batches!(
            batches
                .iter()
                .map(EncodedBsxBatch::height),
            "Batch heights do not match"
        )
        .with_context(|| "Batches must have the same number of rows")?;

        check_eq_batches!(
            batches
                .iter()
                .filter(|b| b.height() > 0)
                .map(|b| EncodedBsxBatch::as_contig(b))
                .map(|c| c.unwrap()),
            "Contigs differ between batches"
        )
        .with_context(|| {
            "All batches must reference the same chromosome/contig"
        })?;

        if let Some(ref labels_vec) = labels {
            if labels_vec.len() != batches.len() {
                return Err(anyhow!(
                    "Number of labels ({}) does not match number of batches \
                     ({})",
                    labels_vec.len(),
                    batches.len()
                ));
            }
        }

        debug!(
            "Successfully created batch group with {} batches",
            batches.len()
        );
        Ok(EncodedBsxBatchGroup { batches, labels })
    }

    /// Creates a new `EncodedBsxBatchGroup` without performing any consistency
    /// checks.
    ///
    /// # Safety
    ///
    /// This function is unsafe because it bypasses the checks performed by
    /// `try_new`. It should only be used when the caller can guarantee
    /// that:
    /// - The input batches are not empty
    /// - All batches have the same height (number of rows)
    /// - All batches have the same contig (chromosome)
    ///
    /// # Arguments
    ///
    /// * `batches` - A vector of `EncodedBsxBatch` instances that are known to
    ///   be consistent.
    ///
    /// # Returns
    ///
    /// A new `EncodedBsxBatchGroup` with no labels.
    unsafe fn new_unchecked(batches: Vec<EncodedBsxBatch>) -> Self {
        debug!(
            "Creating unchecked batch group with {} batches",
            batches.len()
        );
        Self {
            batches,
            labels: None,
        }
    }

    /// Performs an expensive validation check of full positions, strands and
    /// contexts consistency between batches.
    ///
    /// This is a more thorough check than what's done in `try_new`, and
    /// verifies that all batches have identical positions, strands, and
    /// contexts in their data_structs frames.
    ///
    /// # Returns
    ///
    /// Returns `true` if all batches have consistent data_structs, `false`
    /// otherwise.
    fn check_expensive(&self) -> bool {
        debug!(
            "Performing expensive consistency check on {} batches",
            self.batches.len()
        );
        let result = check_eq_batches!(
            self.batches.iter().map(|b| {
                b.data()
                    .select(["position", "strand", "context"])
                    .unwrap()
            }),
            "Data columns (position, strand, context) differ between batches"
        )
        .is_ok();

        if result {
            debug!("Expensive consistency check passed");
        }
        else {
            warn!(
                "Expensive consistency check failed - batches have \
                 inconsistent data_structs"
            );
        }

        result
    }

    /// Applies a Boolean mask to filter all batches in the group.
    ///
    /// This internal method applies the same mask to all batches, maintaining
    /// consistency.
    ///
    /// # Arguments
    ///
    /// * `mask` - A Boolean chunked array specifying which rows to keep (true)
    ///   or filter out (false).
    ///
    /// # Returns
    ///
    /// A new `EncodedBsxBatchGroup` containing only the rows that match the
    /// mask.
    fn filter_mask(
        self,
        mask: &BooleanChunked,
    ) -> anyhow::Result<Self> {
        debug!("Applying filter mask to {} batches", self.batches.len());

        let filtered_batches: Result<Vec<EncodedBsxBatch>, PolarsError> = self
            .batches
            .into_par_iter()
            .map(|batch| batch.filter_mask(&mask))
            .collect();

        let filtered_batches = filtered_batches
            .with_context(|| "Failed to apply filter mask to batches")?;

        info!(
            "Successfully filtered batch group to {} rows",
            if filtered_batches.is_empty() {
                0
            }
            else {
                filtered_batches[0].height()
            }
        );

        Ok(Self {
            batches: filtered_batches,
            labels:  self.labels,
        })
    }

    /// Filters all batches in the group to keep only rows with the specified
    /// genomic context.
    ///
    /// This method computes a Boolean mask once and applies it to all batches,
    /// which is more efficient than filtering each batch individually.
    ///
    /// # Arguments
    ///
    /// * `context` - The genomic context to keep (e.g., CpG, CHG, CHH).
    ///
    /// # Returns
    ///
    /// A new `EncodedBsxBatchGroup` containing only the rows with the specified
    /// context.
    pub fn filter_context(
        self,
        context: Context,
    ) -> anyhow::Result<Self> {
        let context_bool = context.to_bool();
        let df0 = &self.batches[0];

        info!("Filtering batch group by context: {:?}", context);

        let mask = {
            let context_bool = context_bool;
            let mask_iter = df0
                .data()
                .column("context")
                .with_context(|| {
                    "Failed to get 'context' column from batch data_structs"
                })?
                .bool()
                .with_context(|| {
                    "Failed to convert 'context' column to boolean"
                })?
                .iter()
                .map(|v| v == context_bool)
                .collect_vec();
            BooleanChunked::new("mask".into(), mask_iter)
        };

        self.filter_mask(&mask)
    }

    /// Filters all batches in the group to keep only rows with the specified
    /// DNA strand.
    ///
    /// This method computes a Boolean mask once and applies it to all batches,
    /// which is more efficient than filtering each batch individually.
    ///
    /// # Arguments
    ///
    /// * `strand` - The DNA strand to keep (either positive or negative).
    ///
    /// # Returns
    ///
    /// A new `EncodedBsxBatchGroup` containing only the rows with the specified
    /// strand.
    pub fn filter_strand(
        self,
        strand: Strand,
    ) -> anyhow::Result<Self> {
        let strand_bool = strand.to_bool();
        let df0 = &self.batches[0];

        info!("Filtering batch group by strand: {:?}", strand);

        let mask = {
            let mask_iter = df0
                .data()
                .column("strand")
                .with_context(|| {
                    "Failed to get 'strand' column from batch data_structs"
                })?
                .bool()
                .with_context(|| {
                    "Failed to convert 'strand' column to boolean"
                })?
                .iter()
                .map(|v| v == strand_bool)
                .collect_vec();
            BooleanChunked::new("mask".into(), mask_iter)
        };

        self.filter_mask(&mask)
    }

    /// Replaces low coverage methylation sites with count_total = 0 and density
    /// = NaN.
    ///
    /// This method is useful for filtering out unreliable methylation
    /// measurements that don't have sufficient read coverage.
    ///
    /// # Arguments
    ///
    /// * `threshold` - The minimum coverage (read count) required. Sites with
    ///   coverage below this value will be marked as having zero total count
    ///   and NaN density.
    ///
    /// # Returns
    ///
    /// A new `EncodedBsxBatchGroup` with low-coverage sites marked.

    pub fn mark_low_counts(
        self,
        threshold: i16,
    ) -> anyhow::Result<Self> {
        info!(
            "Marking methylation sites with coverage below {} as low-count",
            threshold
        );

        let marked: Result<Vec<EncodedBsxBatch>, PolarsError> = self
            .batches
            .into_par_iter()
            .map(|batch| batch.mark_low_counts(threshold))
            .collect();

        let marked = marked.with_context(|| {
            format!(
                "Failed to mark low count sites with threshold {}",
                threshold
            )
        })?;

        debug!(
            "Successfully marked low-count sites in {} batches",
            marked.len()
        );

        Ok(Self {
            batches: marked,
            labels:  self.labels,
        })
    }

    /// Filters rows where the number of missing (NaN) methylation site
    /// densities is less than or equal to `n_missing`.
    ///
    /// This method is useful for removing sites that have too many missing
    /// values across samples.
    ///
    /// # Arguments
    ///
    /// * `n_missing` - The maximum number of missing values allowed for a site
    ///   to be kept. If a site has more than `n_missing` NaN values across all
    ///   samples, it will be removed.
    ///
    /// # Returns
    ///
    /// A new `EncodedBsxBatchGroup` with sites having too many missing values
    /// removed.
    pub fn filter_n_missing(
        self,
        n_missing: usize,
    ) -> anyhow::Result<Self> {
        let mut threshold = n_missing;
        if threshold > self.n_samples() {
            warn!(
                "Value for filter_n_missing ({}) exceeds sample count ({}). \
                 Adjusting to {}",
                threshold,
                self.n_samples(),
                self.n_samples()
            );
            threshold = self.n_samples();
        }

        let mask = {
            let mask_vec = apply_over_column(
                self.batches.iter().map(EncodedBsxBatch::density).collect_vec(),
                |arr| {
                    arr
                        .map(|x| if x.is_none() {1usize} else {0})
                        .sum_axis(Axis(0))
                        .map(|value| value <= &threshold)
                        .to_vec()
                }
            )?;
            
            BooleanChunked::new("mask".into(), mask_vec)
        };

        self.filter_mask(&mask)
    }

    /// Extracts and returns the data_structs from each batch as a vector of
    /// `DataFrame`s.
    ///
    /// This method consumes the batch group and returns the underlying
    /// data_structs.
    ///
    /// # Returns
    ///
    /// A vector of Polars DataFrames, one for each batch.
    pub fn take_data(self) -> Vec<DataFrame> {
        info!(
            "Extracting {} DataFrames from batch group",
            self.batches.len()
        );
        self.batches
            .into_iter()
            .map(DataFrame::from)
            .collect()
    }

    /// Returns the number of samples (batches) in the group.
    ///
    /// # Returns
    ///
    /// The number of batches in this group.
    #[inline]
    pub fn n_samples(&self) -> usize { self.batches.len() }

    /// Returns the height (number of rows) of the batches in the group.
    ///
    /// Assumes that all batches have the same height, which is enforced by
    /// `try_new`.
    ///
    /// # Returns
    ///
    /// The number of rows in each batch.
    #[inline]
    pub fn height(&self) -> usize { self.batches[0].height() }

    /// Splits the `EncodedBsxBatchGroup` into multiple groups based on the
    /// associated labels.
    ///
    /// This function iterates through the batches and their corresponding
    /// labels (if present), and groups the batches based on the labels.
    /// Batches with the same label are placed into the same group.
    ///
    /// If no labels were provided when creating this group, all batches will be
    /// assigned the default label.
    ///
    /// # Returns
    ///
    /// A HashMap where the keys are the unique labels and the values are
    /// `EncodedBsxBatchGroup` instances containing the batches associated
    /// with that label.
    pub fn split_groups(self) -> HashMap<R, EncodedBsxBatchGroup<R>> {
        info!("Splitting batch group by labels into separate groups");

        let mut out = HashMap::<R, Vec<EncodedBsxBatch>>::new();
        let batches_len = self.batches.len();

        let labels = self.labels.unwrap_or_else(|| {
            debug!("No labels found, using default label for all batches");
            vec![R::default(); batches_len]
        });

        for (batch, label) in izip!(self.batches, labels) {
            out.entry(label.clone())
                .or_insert_with(Vec::new)
                .push(batch);
        }

        let result = HashMap::from_iter(out.into_iter().map(|(k, v)| {
            debug!("Created group for label '{}' with {} batches", k, v.len());
            (k, unsafe { EncodedBsxBatchGroup::new_unchecked(v) })
        }));

        info!("Split into {} distinct groups", result.len());
        result
    }

    /// Returns a reference to the underlying vector of `EncodedBsxBatch`
    /// instances.
    ///
    /// # Returns
    ///
    /// A reference to the vector of batches in this group.
    #[inline]
    pub fn batches(&self) -> &Vec<EncodedBsxBatch> { &self.batches }

    /// Calculates the average methylation density for each methylation site
    /// across all batches.
    ///
    /// # Arguments
    ///
    /// * `na_rm` - If true, missing values (NaN) are ignored in the
    ///   calculation. If false, NaN values result in NaN averages.
    ///
    /// # Returns
    ///
    /// A Result containing a vector of average densities (one for each
    /// methylation site).
    pub fn get_average_density(
        &self,
        na_rm: bool,
    ) -> anyhow::Result<Vec<f32>> {
        apply_over_column::<_, f32, _, _>(
            self.batches.iter().map(EncodedBsxBatch::density).collect_vec(),
            |arr| {
                arr
                    .map(|x| x.unwrap_or(if na_rm { 0.0 } else { f32::NAN }))
                    .mean_axis(Axis(0))
                    .unwrap().to_vec()
            }
        )
    }

    /// Calculates the sum of methylated counts (`count_m`) for each methylation
    /// site across all batches.
    ///
    /// # Returns
    ///
    /// A Result containing a vector of summed methylated counts (one for each
    /// methylation site).
    pub fn get_sum_counts_m(&self) -> anyhow::Result<Vec<i16>> {
        apply_over_column(
            self.batches.iter().map(EncodedBsxBatch::count_m).collect_vec(),
            |arr| {
                arr
                    .map(|x| x.unwrap_or(0))
                    .sum_axis(Axis(0))
                    .to_vec()
            }
        )
    }

    /// Calculates the sum of total counts (`count_total`) for each methylation
    /// site across all batches.
    ///
    /// # Returns
    ///
    /// A Result containing a vector of summed total counts (one for each
    /// methylation site).
    pub fn get_sum_counts_total(&self) -> anyhow::Result<Vec<i16>> {
        apply_over_column(
            self.batches.iter().map(EncodedBsxBatch::count_total).collect_vec(),
            |arr| {
                arr
                    .map(|x| x.unwrap_or(0))
                    .sum_axis(Axis(0))
                    .to_vec()
            }
        )
    }

    /// Retrieves the genomic positions of the methylation sites.
    ///
    /// Assumes that all batches have the same positions, which is enforced by
    /// batch group creation.
    ///
    /// # Returns
    ///
    /// A Result containing a vector of genomic positions (one for each
    /// methylation site).
    pub fn get_positions(&self) -> anyhow::Result<Vec<u32>> {
        debug!("Retrieving genomic positions from batch group");
        self.batches[0]
            .get_position_vals()
            .with_context(|| {
                "Failed to retrieve position values from first batch"
            })
    }

    /// Retrieves the chromosome (contig) name.
    ///
    /// Assumes that all batches have the same chromosome, which is enforced by
    /// batch group creation.
    ///
    /// # Returns
    ///
    /// A Result containing the chromosome name.
    pub fn get_chr(&self) -> anyhow::Result<String> {
        debug!("Retrieving chromosome name from batch group");
        // todo change signature to &str
        self.batches[0].chr_val().map(|v| v.to_string()).with_context(|| {
            "Failed to retrieve chromosome name from first batch"
        })
    }
}
