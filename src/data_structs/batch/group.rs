use std::collections::HashMap;
use std::fmt::Display;
use std::hash::Hash;

use anyhow::{anyhow, Context as AnyhowContext};
use itertools::{izip, Itertools};
use log::warn;
use ndarray::{Array2, Axis};
use polars::datatypes::BooleanChunked;
use polars::error::PolarsError;
use polars::frame::DataFrame;
use polars::prelude::{ChunkedArray, NamedFrom, PolarsNumericType};
use rayon::iter::{IntoParallelIterator,
                  ParallelIterator};

use crate::data_structs::batch::encoded::EncodedBsxBatch;
use crate::data_structs::batch::lazy::LazyBsxBatch;
use crate::data_structs::batch::traits::BsxBatchMethods;
use crate::utils::types::{Context, IPCEncodedEnum, Strand};

/// Checks if all items in a collection are equal
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
    columns: Vec<&ChunkedArray<P>>,    func: F
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

/// Group of encoded BSX batches with optional labels
pub struct EncodedBsxBatchGroup<R: Display + Eq> {
    /// The collection of encoded BSX batches, each representing a sample
    batches: Vec<EncodedBsxBatch>,
    /// Optional labels associated with each batch
    labels:  Option<Vec<R>>,
}

impl<R: Display + Eq> EncodedBsxBatchGroup<R> {
    /// Returns a reference to the optional labels
    pub fn labels(&self) -> &Option<Vec<R>> { &self.labels }
}

#[allow(dead_code)]
impl<R: Display + Eq + Hash + Clone + Default> EncodedBsxBatchGroup<R>
where
    R: Display + Eq + Hash + Clone + Default,
{
    /// Creates a new batch group with validation
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

        Ok(EncodedBsxBatchGroup { batches, labels })
    }

    /// Creates a batch group without consistency checks
    unsafe fn new_unchecked(batches: Vec<EncodedBsxBatch>) -> Self {
        Self {
            batches,
            labels: None,
        }
    }

    /// Checks if all batches have consistent data
    fn check_expensive(&self) -> bool {
        let result = check_eq_batches!(
            self.batches.iter().map(|b| {
                b.data()
                    .select(["position", "strand", "context"])
                    .unwrap()
            }),
            "Data columns (position, strand, context) differ between batches"
        )
        .is_ok();

        result
    }

    /// Applies a Boolean mask to filter all batches
    fn filter_mask(
        self,
        mask: &BooleanChunked,
    ) -> anyhow::Result<Self> {
        let filtered_batches: Result<Vec<EncodedBsxBatch>, PolarsError> = self
            .batches
            .into_par_iter()
            .map(|batch| batch.filter_mask(&mask))
            .collect();

        let filtered_batches = filtered_batches
            .with_context(|| "Failed to apply filter mask to batches")?;

        Ok(Self {
            batches: filtered_batches,
            labels:  self.labels,
        })
    }

    /// Filters batches to keep only rows with specified genomic context
    pub fn filter_context(
        self,
        context: Context,
    ) -> anyhow::Result<Self> {
        let context_bool = context.to_bool();
        let df0 = &self.batches[0];

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

    /// Replaces low coverage methylation sites
    pub fn mark_low_counts(
        self,
        threshold: i16,
    ) -> anyhow::Result<Self> {
        let marked: anyhow::Result<Vec<EncodedBsxBatch>> = self
            .batches
            .into_par_iter()
            .map(
                |batch| LazyBsxBatch::from(batch)
                    .mark_low_coverage(threshold as u32)
                    .try_into()
            )
            .collect();

        let marked = marked.with_context(|| {
            format!(
                "Failed to mark low count sites with threshold {}",
                threshold
            )
        })?;
        Ok(Self {
            batches: marked,
            labels:  self.labels,
        })
    }

    /// Filters rows with too many missing values
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

    /// Extracts the underlying DataFrames
    pub fn take_data(self) -> Vec<DataFrame> {
        self.batches
            .into_iter()
            .map(DataFrame::from)
            .collect()
    }

    /// Returns the number of samples
    #[inline]
    pub fn n_samples(&self) -> usize { self.batches.len() }

    /// Returns the number of rows in each batch
    #[inline]
    pub fn height(&self) -> usize { self.batches[0].height() }

    /// Splits into multiple groups based on labels
    pub fn split_groups(self) -> HashMap<R, EncodedBsxBatchGroup<R>> {
        let mut out = HashMap::<R, Vec<EncodedBsxBatch>>::new();
        let batches_len = self.batches.len();

        let labels = self.labels.unwrap_or_else(|| {
            vec![R::default(); batches_len]
        });

        for (batch, label) in izip!(self.batches, labels) {
            out.entry(label.clone())
                .or_insert_with(Vec::new)
                .push(batch);
        }

        let result = HashMap::from_iter(out.into_iter().map(|(k, v)| {
            (k, unsafe { EncodedBsxBatchGroup::new_unchecked(v) })
        }));

        result
    }

    /// Returns the batches in this group
    #[inline]
    pub fn batches(&self) -> &Vec<EncodedBsxBatch> { &self.batches }

    /// Calculates average methylation density
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

    /// Calculates sum of methylated counts
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

    /// Calculates sum of total counts
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

    /// Gets genomic positions
    pub fn get_positions(&self) -> anyhow::Result<Vec<u32>> {
        self.batches[0]
            .position()
            .to_vec_null_aware()
            .left()
            .ok_or(anyhow!("Nulls in position column"))
            .with_context(|| {
                "Failed to retrieve position values from first batch"
            })
    }

    /// Gets chromosome name
    pub fn get_chr(&self) -> anyhow::Result<String> {
        // todo change signature to &str
        self.batches[0].chr_val().map(|v| v.to_string()).with_context(|| {
            "Failed to retrieve chromosome name from first batch"
        })
    }
}
