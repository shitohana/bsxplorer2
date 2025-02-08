use crate::data_structs::bsx_batch_group::EncodedBsxBatchGroup;
use crate::io::bsx::multiple_reader::MultiBsxFileReader;
use crate::io::bsx::region_read::BsxFileReader;
use crate::utils::types::Context;
use anyhow::{anyhow, Result};
use bio_types::annot::contig::Contig;
use bio_types::annot::refids::RefIDSet;
use bio_types::strand::NoStrand;
use itertools::Itertools;
use log::debug;
use serde::Serialize;
use statrs::distribution::{ContinuousCDF, StudentsT};
use statrs::statistics::Statistics;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufReader, BufWriter, Read, Seek, Write};
use std::sync::Arc;
// ---------------------------------------------------------------------------
// MethyLasso segmentation (fused–lasso) for a single condition.
//
// This function accepts a slice of replicates (each as a Polars DataFrame),
// filters the rows by the specified methylation context (if provided),
// then averages the methylation densities across replicates. The resulting
// signal is denoised using TV denoising (fused–lasso). Finally, the denoised
// signal is scanned to extract contiguous segments (regions with constant
// methylation level).
//
// ---------------------------------------------------------------------------
// The function takes vector of [EncodedBsxBatch] from replicates samples.
// Batches data preferably should not be filtered by context, because (optional)
// filtering is done applying this method. `n_missing` parameter denotes the
// number of possible missing cytosine densities per methylation site. `lambda`
// parameter takes λ for Total Variation Denoising via Condat’s Fast Algorithm.
// This function solves the following problem for a 1D signal y:
//
// minimize_{x} ½ * Σᵢ (yᵢ – xᵢ)² + λ * Σ |xᵢ₊₁ – xᵢ|
//
// The solution is a piecewise–constant approximation of y.

#[macro_export]
macro_rules! generate_stat_method {
    ($name: ident, $field:ident, $statistic:ident) => {
        pub fn $name(&self) -> f64 {
            use statrs::statistics::Statistics;
            self.$field.iter().$statistic()
        }
    };
}

fn segment_tv1d(y: &[f64], positions: &[u32], lambda: f64) -> Result<Vec<u32>> {
    // Apply total variation denoising to the averaged signal.
    let denoised = tv1d::condat(y, lambda);

    // Scan the denoised signal to extract segments (each segment is a contiguous run
    // of (approximately) constant methylation level).
    let mut boundaries = Vec::new();
    if denoised.is_empty() {
        return Ok(boundaries);
    }
    let mut current_value = denoised[0];
    for i in 1..denoised.len() {
        if (denoised[i] - current_value).abs() > f64::EPSILON {
            boundaries.push(positions[i - 1]);
            current_value = denoised[i];
        }
    }
    // Final segment.
    boundaries.push(*positions.last().unwrap());

    Ok(boundaries)
}

#[derive(Debug, Clone, Serialize)]
pub struct MethyLassoConfig {
    pub context: Context,
    pub n_missing: usize,
    pub min_coverage: i16,
    pub diff_threshold: f64,
    pub p_value: f64,
    pub type_density: HashMap<RegionType, f64>,
    pub min_cpgs: usize,
    pub segment_model: SegmentModel,
}

impl Default for MethyLassoConfig {
    fn default() -> Self {
        Self {
            context: Context::CG,
            n_missing: 0,
            min_coverage: 5,
            diff_threshold: 0.05,
            p_value: 0.01,
            type_density: HashMap::from_iter([
                (RegionType::UMR, 0.1),
                (RegionType::LMR, 0.3),
                (RegionType::PMD, 1.0),
            ]),
            min_cpgs: 10,
            segment_model: SegmentModel::default(),
        }
    }
}

impl MethyLassoConfig {
    pub fn finish<F, R>(&self, readers: Vec<(R, F)>) -> Result<MethyLassoDmrIterator<F, R>>
    where
        F: Read + Seek + Send + Sync + 'static,
        R: Display + Eq + Hash + Clone + Default + std::fmt::Debug,
    {
        let sample_mapping: HashMap<uuid::Uuid, (R, F)> = HashMap::from_iter(
            readers
                .into_iter()
                .map(|(group, handle)| (uuid::Uuid::new_v4(), (group, handle))),
        );

        let group_mapping: HashMap<uuid::Uuid, R> = HashMap::from_iter(
            sample_mapping
                .iter()
                .map(|(id, (group, _))| (id.clone(), group.clone())),
        );
        let group_order = {
            let group_set: HashSet<R> =
                HashSet::from_iter(group_mapping.values().map(|x| x.clone()));
            if group_set.len() != 2 {
                return Err(anyhow!(
                    "There should be only two groups of samples! ({:?})",
                    group_set
                ));
            }
            let groups_vec = group_set.into_iter().collect_vec();
            (groups_vec[0].clone(), groups_vec[1].clone())
        };
        let readers_mapping: HashMap<uuid::Uuid, BsxFileReader<F>> = HashMap::from_iter(
            sample_mapping
                .into_iter()
                .map(|(id, (_, handle))| (id, BsxFileReader::new(handle))),
        );
        let multi_reader =
            MultiBsxFileReader::try_new(readers_mapping).map_err(|e| anyhow!("{:?}", e))?;
        let config_copy = self.clone();

        let out = MethyLassoDmrIterator {
            group_mapping,
            multi_reader,
            config: config_copy,
            ref_idset: RefIDSet::new(),
            boundaries: BTreeSet::new(),
            unprocessed_sites: PairwiseSites::new(
                MethylationSites::default(),
                MethylationSites::default(),
            ),
            group_pair: group_order,
            current_chr: Arc::new(String::default()),
            segment_model: self.segment_model.clone(),
        };

        Ok(out)
    }
}

#[derive(Debug, Default, Hash, Copy, Clone, Eq, PartialEq, Serialize)]
pub enum RegionType {
    DMR,
    UMR, // Unmethylated region (or DMV)
    LMR, // Low-methylated region
    PMD, // Partially methylated domain
    #[default]
    Unknown,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MethylationSites {
    pub chr: Arc<String>,
    pub positions: Vec<u32>,
    density: Vec<f64>,
}

impl MethylationSites {
    pub fn drain_until(&mut self, pos: u32) -> Option<MethylationSites> {
        let ref_bound = self.positions.partition_point(|&x| x <= pos);
        if ref_bound == 0 {
            return None;
        }
        Some(Self {
            chr: Arc::clone(&self.chr),
            positions: self.positions.drain(..ref_bound).collect(),
            density: self.density.drain(..ref_bound).collect(),
        })
    }

    pub fn get_chr(&self) -> &Arc<String> {
        &self.chr
    }

    pub fn get_positions(&self) -> &Vec<u32> {
        &self.positions
    }

    pub fn get_density(&self) -> &Vec<f64> {
        &self.density
    }

    pub fn new(chr: Arc<String>, positions: Vec<u32>, density: Vec<f64>) -> Self {
        assert_eq!(positions.len(), density.len());
        assert_ne!(positions.len(), 0);
        Self {
            chr,
            positions,
            density,
        }
    }

    pub fn as_contig(&self) -> Contig<Arc<String>, NoStrand> {
        let start = *self.positions.first().unwrap();
        let end = *self.positions.last().unwrap();
        let len = if start == end { 1 } else { end + 1 - start };
        Contig::new(
            self.chr.clone(),
            start as isize,
            len as usize,
            NoStrand::Unknown,
        )
    }

    pub(crate) fn append(&mut self, other: MethylationSites) {
        let mut other = other;
        self.chr = other.clone().chr;
        self.positions.append(&mut other.positions);
        self.density.append(&mut other.density);
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub(crate) fn drain(&mut self) -> Option<Self> {
        if self.positions.is_empty() {
            return None;
        }
        Some(Self {
            chr: Arc::clone(&self.chr),
            positions: self.positions.drain(..).collect(),
            density: self.density.drain(..).collect(),
        })
    }

    generate_stat_method!(density_mean, density, mean);
    generate_stat_method!(density_var, density, variance);
    generate_stat_method!(density_std, density, std_dev);
}

impl Default for MethylationSites {
    fn default() -> Self {
        Self {
            chr: Arc::new(String::default()),
            positions: Vec::new(),
            density: Vec::new(),
        }
    }
}

impl Iterator for MethylationSites {
    type Item = (u32, f64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.positions.is_empty() {
            return None;
        }
        Some((self.positions.remove(0), self.density.remove(0)))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PairwiseSites {
    pub left: MethylationSites,
    right: MethylationSites,
}

impl PairwiseSites {
    fn new(left: MethylationSites, right: MethylationSites) -> Self {
        Self { left, right }
    }

    fn drain_until(&mut self, pos: u32) -> Option<PairwiseSites> {
        let left_drain = self.left.drain_until(pos)?;
        let right_drain = self.right.drain_until(pos)?;
        Some(PairwiseSites::new(left_drain, right_drain))
    }

    pub fn meth_diff(&self) -> f64 {
        (self.mean_left() - self.mean_right()).abs()
    }

    pub fn mean_left(&self) -> f64 {
        self.left.density.iter().mean()
    }

    pub fn mean_right(&self) -> f64 {
        self.right.density.iter().mean()
    }

    fn compute_t_test_pvalue(&self) -> f64 {
        let n_left = self.left.len();
        let n_right = self.right.len();
        if n_left < 2 || n_right < 2 {
            return 1.0; // not enough data for a t-test
        }
        let (mean_left, var_left) = (self.left.density_mean(), self.left.density_var());
        let (mean_right, var_right) = (self.right.density_mean(), self.right.density_var());

        // If both variances are 0
        if var_left == 0.0 && var_right == 0.0 {
            return if mean_left == mean_right { 1.0 } else { 0.0 };
        }

        let standard_error = (var_left / n_left as f64 + var_right / n_right as f64).sqrt();
        if standard_error == 0.0 {
            return 1.0;
        }
        let t_stat = (mean_left - mean_right) / standard_error;
        // Welch's degrees of freedom:
        let df_num = (var_left / n_left as f64 + var_right / n_right as f64).powi(2);
        let df_den = (var_left.powi(2)) / ((n_left as f64).powi(2) * (n_left as f64 - 1.0))
            + (var_right.powi(2)) / ((n_right as f64).powi(2) * (n_right as f64 - 1.0));
        let df = if df_den != 0.0 { df_num / df_den } else { 1.0 };
        let t_dist = StudentsT::new(0.0, 1.0, df).unwrap();
        let p_value = 2.0 * (1.0 - t_dist.cdf(t_stat.abs()));
        p_value
    }

    fn append(&mut self, other: PairwiseSites) {
        self.left.append(other.left);
        self.right.append(other.right);
    }

    fn drain(&mut self) -> Option<PairwiseSites> {
        let left = self.left.drain()?;
        let right = self.right.drain()?;
        (Self { left, right }).into()
    }

    pub fn into_methylation_region(self, config: &MethyLassoConfig) -> MethylationRegion {
        let mean_left = self.mean_left();
        let mean_right = self.mean_right();
        let p_value = self.compute_t_test_pvalue();
        let meth_diff = (mean_left - mean_right).abs();
        let overall = (mean_left + mean_right) / 2.0;

        let region_type = if p_value <= config.p_value && meth_diff >= config.diff_threshold {
            RegionType::DMR
        } else if overall <= *config.type_density.get(&RegionType::UMR).unwrap_or(&1.0) {
            RegionType::UMR
        } else if overall <= *config.type_density.get(&RegionType::LMR).unwrap_or(&1.0) {
            RegionType::LMR
        } else {
            RegionType::PMD
        };
        MethylationRegion {
            pairwise_sites: self,
            region_type,
            p_value,
            mean_left,
            mean_right,
        }
    }

    pub fn len(&self) -> usize {
        self.left.len()
    }

    pub fn left(&self) -> &MethylationSites {
        &self.left
    }

    pub fn right(&self) -> &MethylationSites {
        &self.right
    }
}

/// Perform Welch’s t–test on two samples and return the two–tailed p–value.
/// (If there isn’t enough data, we return p = 1.0 so that no merge occurs.)
fn welch_t_test(sample1: &[f64], sample2: &[f64]) -> f64 {
    let n1 = sample1.len();
    let n2 = sample2.len();
    if n1 < 2 || n2 < 2 {
        return 1.0;
    }
    let mean1 = sample1.mean();
    let mean2 = sample2.mean();
    let var1 = sample1.variance();
    let var2 = sample2.variance();
    let se = (var1 / (n1 as f64) + var2 / (n2 as f64)).sqrt();
    if se == 0.0 {
        return 1.0;
    }
    let t_stat = (mean1 - mean2).abs() / se;
    let df = (var1 / (n1 as f64) + var2 / (n2 as f64)).powi(2)
        / ((var1 / (n1 as f64)).powi(2) / ((n1 - 1) as f64)
            + (var2 / (n2 as f64)).powi(2) / ((n2 - 1) as f64));
    let t_dist = StudentsT::new(0.0, 1.0, df).unwrap();
    let p_value = 2.0 * (1.0 - t_dist.cdf(t_stat));
    p_value
}

#[derive(Debug, Clone, Serialize)]
pub struct SegmentModel {
    pub min_seg_length: usize,
    pub num_lambdas: usize,
    pub lambda_low: f64,
    pub penalty_weight: f64,
    pub seg_count_weight: f64,
    pub merge_pvalue: f64,
    pub segmentation_tol: f64,
}

impl Default for SegmentModel {
    fn default() -> Self {
        Self {
            min_seg_length: 10,
            num_lambdas: 100,
            lambda_low: 1e-3,
            penalty_weight: 1e5,
            seg_count_weight: 1e4,
            merge_pvalue: 1e-2,
            segmentation_tol: 1e-6,
        }
    }
}

impl SegmentModel {
    /// Optimize λ by grid search (in log–space) over a specified range.
    /// Returns the optimal λ and the corresponding denoised signal.
    ///
    /// You can tune:
    ///   - num_lambdas: grid resolution.
    ///   - lambda_low: lower bound (should be > 0).
    ///   - lambda_high: an upper bound (for example, the overall data range).
    pub fn denoise(&self, signal: &[f64]) -> (f64, Vec<f64>) {
        let mut best_lambda = self.lambda_low;
        let mut best_score = f64::INFINITY;
        let mut best_x = Vec::new();

        let data_min = signal.iter().cloned().fold(f64::INFINITY, f64::min);
        let data_max = signal.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let lambda_high = data_max - data_min; // roughly the overall data range.

        // Grid search in log-space.
        for i in 0..self.num_lambdas {
            let t = i as f64 / (self.num_lambdas - 1) as f64;
            let lambda = self.lambda_low * (lambda_high / self.lambda_low).powf(t);
            let denoised = tv1d::condat(signal, lambda);
            let segments = self.get_segments(&denoised);
            let score = self.objective_function(signal, &denoised, &segments);

            // Uncomment the next line to see the grid search details.
            // println!("λ = {:.6} => {} segments, score = {:.3}", lambda, segments.len(), score);

            if score < best_score {
                best_score = score;
                best_lambda = lambda;
                best_x = denoised;
            }
        }
        (best_lambda, best_x)
    }

    pub fn map_segments(positions: &[u32], segments: Vec<(usize, usize, f64)>) -> Vec<u32> {
        Vec::from_iter(segments.into_iter().map(|(_, end, _)| positions[end]))
    }

    /// Extract segments from a (piecewise–constant) signal.
    /// Two adjacent values are considered equal if their difference is below tol.
    pub fn get_segments(&self, x: &[f64]) -> Vec<(usize, usize, f64)> {
        let mut segments = Vec::new();
        if x.is_empty() {
            return segments;
        }
        let mut start = 0;
        let mut current_val = x[0];
        for (i, &val) in x.iter().enumerate() {
            if (val - current_val).abs() > self.segmentation_tol {
                segments.push((start, i - 1, current_val));
                start = i;
                current_val = val;
            }
        }
        segments.push((start, x.len() - 1, current_val));
        segments
    }

    /// A modified objective function that penalizes both the residual error
    /// and (heavily) any segments shorter than a desired minimum length, as well as
    /// the overall number of segments. (Lower is better.)
    ///
    /// - rss: sum of squared residuals.
    /// - seg_penalty: for each segment shorter than min_seg_length, we add
    ///       penalty_weight * (min_seg_length - seg_length)².
    /// - count_penalty: seg_count_weight * (# segments)
    fn objective_function(&self, y: &[f64], x: &[f64], segments: &[(usize, usize, f64)]) -> f64 {
        let rss: f64 = y
            .iter()
            .zip(x.iter())
            .map(|(yi, xi)| (yi - xi).powi(2))
            .sum();
        let mut seg_penalty = 0.0;
        for &(start, end, _) in segments {
            let seg_len = end - start + 1;
            if seg_len < self.min_seg_length {
                let diff = self.min_seg_length as f64 - seg_len as f64;
                seg_penalty += self.penalty_weight * diff * diff;
            }
        }
        let count_penalty = self.seg_count_weight * segments.len() as f64;
        rss + seg_penalty + count_penalty
    }

    /// Merge adjacent segments if the t–test on their underlying data does not show a significant difference.
    /// The t–test is performed on the original noisy data from each segment.
    /// If the p–value exceeds p_threshold, the segments are merged.
    fn merge_adjacent_segments(
        &self,
        y: &[f64],
        segments: &[(usize, usize, f64)],
    ) -> Vec<(usize, usize, f64)> {
        if segments.is_empty() {
            return vec![];
        }
        let mut merged_segments = Vec::new();
        let mut current_seg = segments[0];

        for seg in segments.iter().skip(1) {
            let sample1 = &y[current_seg.0..=current_seg.1];
            let sample2 = &y[seg.0..=seg.1];
            let p_value = welch_t_test(sample1, sample2);

            // If the difference is not significant, merge the segments.
            if p_value > self.merge_pvalue {
                let new_start = current_seg.0;
                let new_end = seg.1;
                let merged_data = &y[new_start..=new_end];
                let new_mean = merged_data.mean();
                current_seg = (new_start, new_end, new_mean);
            } else {
                merged_segments.push(current_seg);
                current_seg = *seg;
            }
        }
        merged_segments.push(current_seg);
        merged_segments
    }

    /// Iteratively merge segments until no further merging occurs.
    pub fn merge_segments(
        &self,
        y: &[f64],
        segments: &[(usize, usize, f64)],
    ) -> Vec<(usize, usize, f64)> {
        let mut current_segments = segments.to_vec();
        loop {
            let merged = self.merge_adjacent_segments(y, &current_segments);
            if merged.len() == current_segments.len() {
                break;
            }
            current_segments = merged;
        }
        current_segments
    }
}

pub struct MethyLassoDmrIterator<F, R>
where
    F: Read + Seek + Send + Sync + 'static,
    R: Display + Eq + Hash + Clone + Default,
{
    group_mapping: HashMap<uuid::Uuid, R>,
    multi_reader: MultiBsxFileReader<uuid::Uuid, F>,
    segment_model: SegmentModel,
    config: MethyLassoConfig,
    ref_idset: RefIDSet<Arc<String>>,
    boundaries: BTreeSet<u32>,
    current_chr: Arc<String>,
    unprocessed_sites: PairwiseSites,
    group_pair: (R, R),
}

impl<F, R> MethyLassoDmrIterator<F, R>
where
    F: Read + Seek + Send + Sync + 'static,
    R: Display + Eq + Hash + Clone + Default + std::fmt::Debug,
{
    fn read_batch_group(&mut self) -> Option<Result<EncodedBsxBatchGroup<R>>> {
        if let Some(batches) = self.multi_reader.next() {
            let mut data = Vec::new();
            let mut labels = Vec::new();

            for (id, batch) in batches {
                labels.push(self.group_mapping.get(&id).unwrap().clone());
                data.push(batch)
            }

            let batch_group = EncodedBsxBatchGroup::try_new(data, Some(labels));
            Some(batch_group)
        } else {
            None
        }
    }

    pub fn blocks_total(&self) -> usize {
        self.multi_reader.blocks_total()
    }

    pub fn current_block(&self) -> usize {
        self.multi_reader.current_batch_idx()
    }

    fn process_group(&mut self, group: EncodedBsxBatchGroup<R>) -> Result<()> {
        let segment_model = SegmentModel::default();

        // 1. Check correct groups
        if !(group
            .labels()
            .as_ref()
            .map(|groups| groups.iter().unique().count() == 2)
            .unwrap_or(false))
        {
            return Err(anyhow!(
                "There should be EXACTLY two sample groups! {:?}",
                group.labels()
            ));
        }
        // 2. Apply filters
        let group = group
            .filter_context(self.config.context)?
            .mark_low_counts(self.config.min_coverage)?
            .filter_n_missing(self.config.n_missing)?;

        if group.height() == 0 {
            return Ok(());
        }

        // 3. Divide groups
        let individual_groups = group.split_groups();
        let group_left = individual_groups.get(&self.group_pair.0).unwrap();
        let group_right = individual_groups.get(&self.group_pair.1).unwrap();
        // 6. Extract constant vars
        let chr = self.ref_idset.intern(group_left.get_chr()?.as_str());
        let positions = {
            self.unprocessed_sites
                .left
                .positions
                .clone()
                .into_iter()
                .chain(group_left.get_positions()?.into_iter())
                .collect_vec()
        };
        // 4. Calculate avg densities
        let avg_density_left = {
            let new_vals = group_left.get_average_density(true)?;
            self.unprocessed_sites
                .left
                .drain()
                .map(|m| m.density)
                .unwrap_or(Vec::new())
                .into_iter()
                .chain(new_vals.into_iter())
                .collect_vec()
        };
        let avg_density_right = {
            let new_vals = group_right.get_average_density(true)?;
            self.unprocessed_sites
                .right
                .drain()
                .map(|m| m.density)
                .unwrap_or(Vec::new())
                .into_iter()
                .chain(new_vals.into_iter())
                .collect_vec()
        };

        let mut boundaries_union = {
            let (_, denoised_left) = segment_model.denoise(&avg_density_left);
            let (_, denoised_right) = segment_model.denoise(&avg_density_right);

            let initial_seg_left = segment_model.get_segments(&denoised_left);
            let initial_seg_right = segment_model.get_segments(&denoised_right);

            let merged_left = segment_model.merge_segments(&avg_density_left, &initial_seg_left);
            let merged_right = segment_model.merge_segments(&avg_density_right, &initial_seg_right);

            let boundaries_left =
                BTreeSet::from_iter(SegmentModel::map_segments(&positions, merged_left));
            let boundaries_right =
                BTreeSet::from_iter(SegmentModel::map_segments(&positions, merged_right));
            let boundaries_union =
                BTreeSet::from_iter(boundaries_left.union(&boundaries_right).cloned());

            boundaries_union
        };

        // 11. Create pairwise sites
        let pairwise = PairwiseSites::new(
            MethylationSites::new(chr.clone(), positions.clone(), avg_density_left),
            MethylationSites::new(chr.clone(), positions.clone(), avg_density_right),
        );

        self.boundaries.append(&mut boundaries_union);
        self.unprocessed_sites.append(pairwise);
        self.current_chr = chr;

        debug!("Processed new methylation batch");
        Ok(())
    }
}

#[derive(Debug)]
pub struct MethylationRegion {
    pub pairwise_sites: PairwiseSites,
    pub region_type: RegionType,
    pub p_value: f64,
    pub mean_left: f64,
    pub mean_right: f64,
}

impl Display for MethylationRegion {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let meth_diff = (self.mean_right - self.mean_left).abs();
        let overall = (self.mean_left + self.mean_right) / 2.0;
        let n_cytosines = self.pairwise_sites.len();
        let ref_id = self.pairwise_sites.left.chr.as_str();
        let start = self.pairwise_sites.left.positions.first().unwrap();
        let end = self.pairwise_sites.left.positions.last().unwrap();
        write!(
            f,
            "{:?}\t| {ref_id}:{start}-{end}
\t| N cytosines: {n_cytosines}
\t| Overall density: {overall:.3}
\t| Methylation difference: {meth_diff:.5}
\t| Methylation mean: A: {:.3}  B: {:.3} 
\t| DMR p-value: {:.3e}\n",
            self.region_type, self.mean_left, self.mean_right, self.p_value
        )?;
        Ok(())
    }
}

impl MethylationRegion {
    pub fn pairwise_sites(&self) -> &PairwiseSites {
        &self.pairwise_sites
    }

    pub fn region_type(&self) -> RegionType {
        self.region_type
    }
}

impl<F, R> Iterator for MethyLassoDmrIterator<F, R>
where
    F: Read + Seek + Send + Sync + 'static,
    R: Display + Eq + Hash + Clone + Default + std::fmt::Debug,
{
    type Item = MethylationRegion;

    fn next(&mut self) -> Option<Self::Item> {
        // If enough boundaries, return data
        let pairwise_data = if self.boundaries.len() > 1 {
            let boundary = self.boundaries.pop_first().unwrap();
            self.unprocessed_sites
                .drain_until(boundary)
                .expect("Bound to drained methylation data")
        // If it is last boundary,
        } else {
            // Try to fill the buffer
            if let Some(new_group) = self.read_batch_group() {
                let new_group = new_group.expect("Failed to read new group");

                let chr_differ = new_group
                    .get_chr()
                    .map(|val| self.current_chr != self.ref_idset.intern(val.as_str()))
                    .unwrap_or(true);

                // Return the rest of unprocessed
                if chr_differ {
                    let pairwise_data = if let Some(data) = self.unprocessed_sites.drain() {
                        data
                    } else {
                        self.boundaries.clear();
                        self.process_group(new_group)
                            .expect("Failed to process new group");

                        return self.next();
                    };

                    pairwise_data
                // Process new data and try again
                } else {
                    self.process_group(new_group)
                        .expect("Failed to process new group");
                    return self.next();
                }
            // Filling buffer failed
            } else {
                // If some unprocessed left
                if !self.boundaries.is_empty() {
                    self.boundaries.clear();
                    if let Some(data) = self.unprocessed_sites.drain() {
                        data
                    } else {
                        return self.next();
                    }
                // Fully released
                } else {
                    return None;
                }
            }
        };
        if pairwise_data.len() < self.config.min_cpgs {
            return self.next();
        }

        Some(pairwise_data.into_methylation_region(&self.config))
    }
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct MethylLassoRunConfig {
    pub analysis_config: MethyLassoConfig,
    pub sample_paths: Vec<String>,
    pub sample_labels: Vec<String>,
    pub selected_regions: Vec<RegionType>,
    pub output: String,
}

impl MethylLassoRunConfig {
    pub fn run(&self) -> Result<()> {
        let files = self
            .sample_paths
            .iter()
            .map(|path| File::open(path))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .map(BufReader::new)
            .collect_vec();

        let iterator = self.analysis_config.clone().finish(
            self.sample_labels
                .iter()
                .cloned()
                .zip(files.into_iter())
                .collect_vec(),
        )?;

        let mut sink = BufWriter::new(File::create(&self.output)?);

        for region in iterator {
            if self.selected_regions.contains(&region.region_type) {
                let meth_diff = region.mean_right - region.mean_left;
                let overall = (region.mean_left + region.mean_right) / 2.0;
                let n_cytosines = region.pairwise_sites.len();
                let ref_id = region.pairwise_sites.left.chr.as_str();
                let start = region.pairwise_sites.left.positions.first().unwrap();
                let end = region.pairwise_sites.left.positions.last().unwrap();
                let region_type = region.region_type;
                let mean_a = region.mean_left;
                let mean_b = region.mean_right;
                let pvalue = region.p_value;

                let line = format!("{ref_id}\t{start}\t{end}\t{n_cytosines}\t{region_type:?}\t{pvalue:.3e}\t{mean_a:.5}\t{mean_b:.5}\t{meth_diff:.5}\t{overall}");
                writeln!(sink, "{line}")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio_types::annot::loc::Loc;
    use std::sync::Arc;

    #[test]
    fn test_creation() {
        let chr = Arc::new("chr1".to_string());
        let positions = vec![1, 2, 3, 4, 5];
        let density = vec![0.1, 0.2, 0.3, 0.4, 0.5];
        let sites = MethylationSites::new(chr.clone(), positions.clone(), density.clone());

        assert_eq!(sites.get_chr(), &chr);
        assert_eq!(sites.get_positions(), &positions);
        assert_eq!(sites.get_density(), &density);
    }

    #[test]
    #[should_panic]
    fn test_creation_with_mismatched_lengths() {
        let chr = Arc::new("chr1".to_string());
        let positions = vec![1, 2, 3];
        let density = vec![0.1, 0.2]; // Mismatched length
        MethylationSites::new(chr, positions, density);
    }

    #[test]
    #[should_panic]
    fn test_creation_with_empty_data() {
        let chr = Arc::new("chr1".to_string());
        let positions = vec![];
        let density = vec![];
        MethylationSites::new(chr, positions, density);
    }

    #[test]
    fn test_drain_until() {
        let chr = Arc::new("chr1".to_string());
        let positions = vec![1, 2, 3, 4, 5, 10, 15];
        let density = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
        let mut sites = MethylationSites::new(chr.clone(), positions, density);

        let drained = sites.drain_until(5).unwrap();
        assert_eq!(drained.get_positions(), &vec![1, 2, 3, 4, 5]);
        assert_eq!(drained.get_density(), &vec![0.1, 0.2, 0.3, 0.4, 0.5]);
        assert_eq!(sites.get_positions(), &vec![10, 15]);
        assert_eq!(sites.get_density(), &vec![0.6, 0.7]);
    }

    #[test]
    fn test_drain_until_no_drain() {
        let chr = Arc::new("chr1".to_string());
        let positions = vec![10, 15, 20];
        let density = vec![0.6, 0.7, 0.8];
        let mut sites = MethylationSites::new(chr.clone(), positions, density);

        let drained = sites.drain_until(5);
        assert!(drained.is_none());
        assert_eq!(sites.get_positions(), &vec![10, 15, 20]);
        assert_eq!(sites.get_density(), &vec![0.6, 0.7, 0.8]);
    }

    #[test]
    fn test_append() {
        let chr = Arc::new("chr1".to_string());
        let mut sites1 = MethylationSites::new(chr.clone(), vec![1, 2, 3], vec![0.1, 0.2, 0.3]);
        let sites2 = MethylationSites::new(chr.clone(), vec![4, 5], vec![0.4, 0.5]);

        sites1.append(sites2);
        assert_eq!(sites1.get_positions(), &vec![1, 2, 3, 4, 5]);
        assert_eq!(sites1.get_density(), &vec![0.1, 0.2, 0.3, 0.4, 0.5]);
    }

    #[test]
    fn test_drain() {
        let chr = Arc::new("chr1".to_string());
        let mut sites = MethylationSites::new(chr.clone(), vec![1, 2, 3], vec![0.1, 0.2, 0.3]);

        let drained = sites.drain().unwrap();
        assert_eq!(drained.get_positions(), &vec![1, 2, 3]);
        assert_eq!(drained.get_density(), &vec![0.1, 0.2, 0.3]);
        assert!(sites.get_positions().is_empty());
        assert!(sites.get_density().is_empty());
    }

    #[test]
    fn test_len() {
        let chr = Arc::new("chr1".to_string());
        let sites = MethylationSites::new(
            chr.clone(),
            vec![1, 2, 3, 4, 5],
            vec![0.1, 0.2, 0.3, 0.4, 0.5],
        );

        assert_eq!(sites.len(), 5);
    }

    #[test]
    fn test_as_contig() {
        let chr = Arc::new("chr1".to_string());
        let sites = MethylationSites::new(chr.clone(), vec![10, 20, 30], vec![0.1, 0.2, 0.3]);

        let contig = sites.as_contig();
        assert_eq!(contig.refid(), &chr);
        assert_eq!(contig.start(), 10);
        assert_eq!(contig.length(), 21); // 30 - 10 + 1
        assert_eq!(contig.strand(), NoStrand::Unknown);
    }

    #[test]
    fn test_density_statistics() {
        let chr = Arc::new("chr1".to_string());
        let sites = MethylationSites::new(
            chr.clone(),
            vec![1, 2, 3, 4, 5],
            vec![0.1, 0.2, 0.3, 0.4, 0.5],
        );

        let mean = sites.density_mean();
        let variance = sites.density_var();
        let std_dev = sites.density_std();

        assert!(
            (mean - 0.3).abs() < 1e-6,
            "Expected mean ≈ 0.3, got {}",
            mean
        );
        assert!(
            (variance - 0.025).abs() < 1e-6,
            "Expected variance ≈ 0.025, got {}",
            variance
        );
        assert!(
            (std_dev - 0.1581).abs() < 1e-4,
            "Expected std dev ≈ 0.1581, got {}",
            std_dev
        );
    }

    use assert_approx_eq::assert_approx_eq;

    /// Helper function to create a dummy `MethylationSites`
    fn create_methylation_sites(
        chr: &str,
        positions: Vec<u32>,
        density: Vec<f64>,
    ) -> MethylationSites {
        MethylationSites::new(Arc::new(chr.to_string()), positions, density)
    }

    /// Helper function to create a dummy `PairwiseSites`
    fn create_pairwise_sites() -> PairwiseSites {
        let left =
            create_methylation_sites("chr1", vec![1, 2, 3, 4, 5], vec![0.1, 0.2, 0.3, 0.4, 0.5]);
        let right =
            create_methylation_sites("chr1", vec![1, 2, 3, 4, 5], vec![0.2, 0.3, 0.4, 0.5, 0.6]);
        PairwiseSites::new(left, right)
    }

    #[test]
    fn test_pairwise_sites_creation() {
        let pairwise = create_pairwise_sites();
        assert_eq!(pairwise.left.get_positions(), &vec![1, 2, 3, 4, 5]);
        assert_eq!(pairwise.right.get_positions(), &vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_meth_diff() {
        let pairwise = create_pairwise_sites();
        let diff = pairwise.meth_diff();
        assert_approx_eq!(diff, 0.1, 1e-6);
    }

    #[test]
    fn test_mean_calculations() {
        let pairwise = create_pairwise_sites();
        assert_approx_eq!(pairwise.mean_left(), 0.3, 1e-6);
        assert_approx_eq!(pairwise.mean_right(), 0.4, 1e-6);
    }

    #[test]
    fn test_compute_t_test_pvalue() {
        let pairwise = create_pairwise_sites();
        let p_value = pairwise.compute_t_test_pvalue();
        assert!(
            p_value >= 0.0 && p_value <= 1.0,
            "P-value should be between 0 and 1."
        );
    }

    #[test]
    fn test_append_pairwise() {
        let mut pairwise1 = create_pairwise_sites();
        let pairwise2 = create_pairwise_sites();
        pairwise1.append(pairwise2);

        assert_eq!(pairwise1.left.get_positions().len(), 10);
        assert_eq!(pairwise1.right.get_positions().len(), 10);
    }

    #[test]
    fn test_drain_until_pairwise() {
        let mut pairwise = create_pairwise_sites();
        let drained = pairwise.drain_until(3).unwrap();
        assert_eq!(drained.left.get_positions(), &vec![1, 2, 3]);
        assert_eq!(drained.right.get_positions(), &vec![1, 2, 3]);
        assert_eq!(pairwise.left.get_positions(), &vec![4, 5]);
        assert_eq!(pairwise.right.get_positions(), &vec![4, 5]);
    }

    #[test]
    fn test_drain_pairwise() {
        let mut pairwise = create_pairwise_sites();
        let _ = pairwise.drain().unwrap();
        assert!(pairwise.left.get_positions().is_empty());
        assert!(pairwise.right.get_positions().is_empty());
    }

    #[test]
    fn test_into_methylation_region() {
        let config = MethyLassoConfig::default();

        let pairwise = create_pairwise_sites();
        let methylation_region = pairwise.into_methylation_region(&config);

        assert!(matches!(methylation_region.region_type, RegionType::PMD));
    }

    fn default_config() -> MethyLassoConfig {
        MethyLassoConfig::default()
    }

    /// Helper function to create a `PairwiseSites` with specified left and right densities
    fn create_pairwise_sites_custom(
        left_density: Vec<f64>,
        right_density: Vec<f64>,
    ) -> PairwiseSites {
        let positions = vec![1, 2, 3, 4, 5]; // Shared genomic positions
        let left = create_methylation_sites("chr1", positions.clone(), left_density);
        let right = create_methylation_sites("chr1", positions, right_density);
        PairwiseSites::new(left, right)
    }

    #[test]
    fn test_dmr_classification() {
        let config = default_config();
        // DMR: Significant p-value and methylation difference above threshold
        let pairwise = create_pairwise_sites_custom(
            vec![0.9, 0.9, 0.9, 0.9, 0.9],
            vec![0.1, 0.1, 0.1, 0.1, 0.1],
        );
        let region = pairwise.into_methylation_region(&config);

        assert_eq!(region.region_type, RegionType::DMR);
    }

    #[test]
    fn test_umr_classification() {
        let config = default_config();
        // UMR: Both left and right methylation averages ≤ 0.1
        let pairwise = create_pairwise_sites_custom(
            vec![0.05, 0.08, 0.09, 0.06, 0.07],
            vec![0.08, 0.07, 0.06, 0.09, 0.05],
        );
        let region = pairwise.into_methylation_region(&config);

        assert_eq!(region.region_type, RegionType::UMR);
    }

    #[test]
    fn test_lmr_classification() {
        let config = default_config();
        // LMR: Overall average methylation between 0.1 and 0.3
        let pairwise = create_pairwise_sites_custom(
            vec![0.2, 0.25, 0.27, 0.22, 0.24],
            vec![0.2, 0.21, 0.23, 0.25, 0.22],
        );
        let region = pairwise.into_methylation_region(&config);

        assert_eq!(region.region_type, RegionType::LMR);
    }

    #[test]
    fn test_pmd_classification() {
        let config = default_config();
        // PMD: Overall average methylation above 0.3
        let pairwise = create_pairwise_sites_custom(
            vec![0.5, 0.6, 0.7, 0.65, 0.55],
            vec![0.4, 0.45, 0.5, 0.48, 0.52],
        );
        let region = pairwise.into_methylation_region(&config);

        assert_eq!(region.region_type, RegionType::PMD);
    }
}
