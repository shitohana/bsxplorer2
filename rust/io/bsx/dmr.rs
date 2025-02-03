use crate::bsx_batch_group::EncodedBsxBatchGroup;
use crate::io::bsx::multiple_reader::MultiBsxFileReader;
use crate::io::bsx::region_read::BsxFileReader;
use crate::utils::types::Context;
use anyhow::{anyhow, Result};
use bio_types::annot::contig::Contig;
use bio_types::annot::refids::RefIDSet;
use bio_types::strand::NoStrand;
use itertools::Itertools;
use log::debug;
use statrs::distribution::{ContinuousCDF, StudentsT};
use statrs::statistics::Statistics;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::hash::Hash;
use std::io::{Read, Seek};
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

#[derive(Debug, Clone)]
pub struct MethyLassoConfig {
    pub context: Context,
    pub n_missing: usize,
    pub min_coverage: i16,
    pub lambda: f64,
    pub diff_threshold: f64,
    pub p_value: f64,
    pub type_density: HashMap<RegionType, f64>,
    pub min_cpgs: usize,
}

impl Default for MethyLassoConfig {
    fn default() -> Self {
        Self {
            context: Context::CG,
            n_missing: 0,
            min_coverage: 5,
            lambda: 3.0,
            diff_threshold: 0.05,
            p_value: 0.01,
            type_density: HashMap::from_iter([
                (RegionType::UMR, 0.1),
                (RegionType::LMR, 0.3),
                (RegionType::PMD, 1.0),
            ]),
            min_cpgs: 10,
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
        };

        Ok(out)
    }
}

#[derive(Debug, Default, Hash, Copy, Clone, Eq, PartialEq)]
pub enum RegionType {
    DMR,
    UMR, // Unmethylated region (or DMV)
    LMR, // Low-methylated region
    PMD, // Partially methylated domain
    #[default]
    Unknown,
}

#[derive(Debug, Clone, PartialEq)]
struct MethylationSites {
    chr: Arc<String>,
    positions: Vec<u32>,
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

    pub fn append(&mut self, other: MethylationSites) {
        let mut other = other;
        self.chr = other.clone().chr;
        self.positions.append(&mut other.positions);
        self.density.append(&mut other.density);
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn drain(&mut self) -> Option<Self> {
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
    left: MethylationSites,
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

    fn meth_diff(&self) -> f64 {
        (self.mean_left() - self.mean_right()).abs()
    }

    fn mean_left(&self) -> f64 {
        self.left.density.iter().mean()
    }

    fn mean_right(&self) -> f64 {
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

    fn into_methylation_region(self, config: &MethyLassoConfig) -> MethylationRegion {
        let mean_left = self.mean_left();
        let mean_right = self.mean_right();
        let p_value = self.compute_t_test_pvalue();
        let meth_diff = (mean_left - mean_right).abs();
        let overall = (mean_left + mean_right) / 2.0;

        let region_type = if p_value <= config.p_value && meth_diff >= config.diff_threshold {
            RegionType::DMR
        } else if overall <= *config.type_density.get(&RegionType::UMR).unwrap() {
            RegionType::UMR
        } else if overall <= *config.type_density.get(&RegionType::LMR).unwrap() {
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
}

pub struct MethyLassoDmrIterator<F, R>
where
    F: Read + Seek + Send + Sync + 'static,
    R: Display + Eq + Hash + Clone + Default,
{
    group_mapping: HashMap<uuid::Uuid, R>,
    multi_reader: MultiBsxFileReader<uuid::Uuid, F>,
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

    fn process_group(&mut self, group: EncodedBsxBatchGroup<R>) -> Result<()> {
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
        let mut positions = self.unprocessed_sites.left.positions.clone();
        positions.append(&mut group_left.get_positions()?);
        // 4. Calculate avg densities
        let mut avg_density_left = self
            .unprocessed_sites
            .left
            .drain()
            .map(|m| m.density)
            .unwrap_or(Vec::new());
        avg_density_left.append(&mut group_left.get_average_density(true)?);
        let mut avg_density_right = self
            .unprocessed_sites
            .right
            .drain()
            .map(|m| m.density)
            .unwrap_or(Vec::new());
        avg_density_right.append(&mut group_right.get_average_density(true)?);
        // 7. Get segment boundaries
        let boundaries_left = BTreeSet::from_iter(segment_tv1d(
            &avg_density_left,
            &positions,
            self.config.lambda,
        )?);
        let boundaries_right = BTreeSet::from_iter(segment_tv1d(
            &avg_density_left,
            &positions,
            self.config.lambda,
        )?);
        // 8. Union boundaries
        let mut boundaries_union =
            BTreeSet::from_iter(boundaries_left.union(&boundaries_right).cloned());
        // 9. Cache boundaries
        self.boundaries.append(&mut boundaries_union);
        // 10. Create methylation stats objects
        let sites_left = MethylationSites::new(chr.clone(), positions.clone(), avg_density_left);
        let sites_right = MethylationSites::new(chr.clone(), positions.clone(), avg_density_right);
        // 11. Create pairwise sites
        let pairwise = PairwiseSites::new(sites_left, sites_right);
        self.unprocessed_sites.append(pairwise);

        self.current_chr = chr;
        debug!("Processed new methylation batch");
        Ok(())
    }
}

#[derive(Debug)]
pub struct MethylationRegion {
    pairwise_sites: PairwiseSites,
    region_type: RegionType,
    p_value: f64,
    mean_left: f64,
    mean_right: f64,
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
