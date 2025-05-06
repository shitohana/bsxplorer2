use std::collections::BTreeMap;
use std::fmt::Write;

use anyhow::{Context as AnyhowContext, Result};
use hashbrown::HashMap;
use itertools::Itertools;
use log::{debug, info, trace};
use serde::{Deserialize, Serialize, Serializer};

use crate::data_structs::enums::{Context, Strand};

/// Serializes a HashMap in a deterministic order based on keys.
/// This is useful for consistent JSON outputs and testing.
fn serialize_sorted_map<S, K: Ord + Serialize, V: Serialize>(
    map: &HashMap<K, V>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer, {
    let sorted_map: BTreeMap<_, _> = map.iter().collect();
    sorted_map.serialize(serializer)
}

/// Represents comprehensive methylation statistics for DNA sequencing
/// data_structs.
///
/// Stores information about methylation levels, variance, coverage
/// distribution, and context-specific (CG, CHG, CHH) and strand-specific
/// methylation data_structs.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MethylationStats {
    /// Overall mean methylation level across all analyzed positions
    mean_methylation: f64,

    /// Variance of methylation levels across all analyzed positions
    methylation_var: f64,

    /// Maps coverage levels to their frequency
    /// Key: coverage depth, Value: number of positions with that coverage
    #[serde(serialize_with = "serialize_sorted_map")]
    coverage_distribution: HashMap<u16, u32>,

    /// Methylation statistics per sequence context
    /// Key: context (CG, CHG, CHH), Value: (sum of methylation levels, count
    /// of positions)
    #[serde(serialize_with = "serialize_sorted_map")]
    context_methylation: HashMap<Context, (f64, u32)>,

    /// Methylation statistics per DNA strand
    /// Key: strand (Forward, Reverse), Value: (sum of methylation levels,
    /// count of positions)
    #[serde(serialize_with = "serialize_sorted_map")]
    strand_methylation: HashMap<Strand, (f64, u32)>,
}

/// A flattened representation of methylation statistics.
///
/// This structure provides direct access to key statistics without the need to
/// query maps, making it suitable for summary outputs and reporting.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MethylationStatFlat {
    /// Overall mean methylation level
    pub mean_methylation: f64,

    /// Overall variance in methylation levels
    pub methylation_var: f64,

    /// Average coverage across all positions
    pub mean_coverage: f64,

    /// Mean methylation level in CG context
    pub cg_mean_methylation: f64,

    /// Number of positions in CG context
    pub cg_coverage: u32,

    /// Mean methylation level in CHG context
    pub chg_mean_methylation: f64,

    /// Number of positions in CHG context
    pub chg_coverage: u32,

    /// Mean methylation level in CHH context
    pub chh_mean_methylation: f64,

    /// Number of positions in CHH context
    pub chh_coverage: u32,

    /// Mean methylation level on forward strand
    pub fwd_mean_methylation: f64,

    /// Number of positions on forward strand
    pub fwd_coverage: u32,

    /// Mean methylation level on reverse strand
    pub rev_mean_methylation: f64,

    /// Number of positions on reverse strand
    pub rev_coverage: u32,
}

impl From<MethylationStats> for MethylationStatFlat {
    /// Converts detailed methylation statistics into a flattened format
    /// for easier reporting and analysis.
    fn from(value: MethylationStats) -> Self {
        trace!("Converting MethylationStats to MethylationStatFlat");

        let mean_coverage = if !value.coverage_distribution.is_empty() {
            value.total_coverage() as f64
                / value.coverage_distribution.len() as f64
        }
        else {
            0.0
        };

        MethylationStatFlat {
            mean_methylation: value.mean_methylation,
            methylation_var: value.methylation_var,
            mean_coverage,
            cg_mean_methylation: value
                .context_methylation
                .get(&Context::CG)
                .map_or(0.0, |(sum_m, count)| {
                    if *count > 0 {
                        sum_m / *count as f64
                    }
                    else {
                        0.0
                    }
                }),
            cg_coverage: value
                .context_methylation
                .get(&Context::CG)
                .map_or(0, |(_, count)| *count),
            chg_mean_methylation: value
                .context_methylation
                .get(&Context::CHG)
                .map_or(0.0, |(sum_m, count)| {
                    if *count > 0 {
                        sum_m / *count as f64
                    }
                    else {
                        0.0
                    }
                }),
            chg_coverage: value
                .context_methylation
                .get(&Context::CHG)
                .map_or(0, |(_, count)| *count),
            chh_mean_methylation: value
                .context_methylation
                .get(&Context::CHH)
                .map_or(0.0, |(sum_m, count)| {
                    if *count > 0 {
                        sum_m / *count as f64
                    }
                    else {
                        0.0
                    }
                }),
            chh_coverage: value
                .context_methylation
                .get(&Context::CHH)
                .map_or(0, |(_, count)| *count),
            fwd_mean_methylation: value
                .strand_methylation
                .get(&Strand::Forward)
                .map_or(0.0, |(sum_m, count)| {
                    if *count > 0 {
                        sum_m / *count as f64
                    }
                    else {
                        0.0
                    }
                }),
            fwd_coverage: value
                .strand_methylation
                .get(&Strand::Forward)
                .map_or(0, |(_, count)| *count),
            rev_mean_methylation: value
                .strand_methylation
                .get(&Strand::Reverse) // Fixed: was using Forward before
                .map_or(0.0, |(sum_m, count)| {
                    if *count > 0 {
                        sum_m / *count as f64
                    } else {
                        0.0
                    }
                }),
            rev_coverage: value
                .strand_methylation
                .get(&Strand::Reverse)
                .map_or(0, |(_, count)| *count),
        }
    }
}

impl MethylationStats {
    /// Creates an empty methylation statistics instance with all values
    /// initialized to zero.
    ///
    /// This is useful as a starting point before aggregating data_structs from
    /// multiple sources.
    pub fn new() -> Self {
        debug!("Creating new empty MethylationStats");
        Self {
            mean_methylation:      0.0,
            methylation_var:       0.0,
            coverage_distribution: HashMap::new(),
            context_methylation:   HashMap::new(),
            strand_methylation:    HashMap::new(),
        }
    }

    /// Creates a MethylationStats instance from pre-calculated statistics.
    ///
    /// # Arguments
    ///
    /// * `mean_methylation` - The overall mean methylation level
    /// * `variance_methylation` - The variance in methylation levels
    /// * `coverage_distribution` - HashMap mapping coverage depths to their
    ///   frequencies
    /// * `context_methylation` - HashMap mapping sequence contexts to their
    ///   methylation statistics
    /// * `strand_methylation` - HashMap mapping DNA strands to their
    ///   methylation statistics
    ///
    /// # Returns
    ///
    /// A new MethylationStats instance with the provided data_structs, handling
    /// NaN values gracefully.
    pub fn from_data(
        mean_methylation: f64,
        variance_methylation: f64,
        coverage_distribution: HashMap<u16, u32>,
        context_methylation: HashMap<Context, (f64, u32)>,
        strand_methylation: HashMap<Strand, (f64, u32)>,
    ) -> MethylationStats {
        debug!(
            "Creating MethylationStats from data_structs with mean: {}, \
             variance: {}",
            mean_methylation, variance_methylation
        );

        // Handle NaN values by replacing them with 0.0
        let mean_methylation = if mean_methylation.is_nan() {
            debug!("Detected NaN in mean_methylation, replacing with 0.0");
            0.0
        }
        else {
            mean_methylation
        };

        let variance_methylation = if variance_methylation.is_nan() {
            debug!("Detected NaN in variance_methylation, replacing with 0.0");
            0.0
        }
        else {
            variance_methylation
        };

        Self {
            mean_methylation,
            methylation_var: variance_methylation,
            context_methylation,
            coverage_distribution,
            strand_methylation,
        }
    }

    /// Merges another MethylationStats instance into this one.
    ///
    /// This performs a weighted merge of statistics, with weights based on
    /// coverage. The method updates mean methylation, variance, coverage
    /// distribution, and context/strand-specific statistics.
    ///
    /// # Arguments
    ///
    /// * `other` - The MethylationStats instance to merge into this one
    pub fn merge(
        &mut self,
        other: &MethylationStats,
    ) {
        let weight_self = self.total_coverage() as f64;
        let weight_other = other.total_coverage() as f64;
        let total_weight = weight_self + weight_other;

        debug!(
            "Merging MethylationStats (self weight: {}, other weight: {})",
            weight_self, weight_other
        );

        // Weighted mean methylation calculation
        if total_weight > 0.0 {
            let delta = self.mean_methylation - other.mean_methylation;
            self.mean_methylation = (weight_self * self.mean_methylation
                + weight_other * other.mean_methylation)
                / total_weight;

            // Correct variance computation considering mean shift
            self.methylation_var = ((weight_self * self.methylation_var
                + weight_other * other.methylation_var)
                + (weight_self * weight_other / total_weight)
                    * (delta * delta))
                / total_weight;

            trace!(
                "Updated mean methylation: {}, variance: {}",
                self.mean_methylation,
                self.methylation_var
            );
        }
        else {
            debug!(
                "Skipping mean/variance calculations due to zero total weight"
            );
        }

        // Merge coverage distribution
        for (&count_total, &frequency) in &other.coverage_distribution {
            *self
                .coverage_distribution
                .entry(count_total)
                .or_insert(0) += frequency;
            trace!(
                "Updated coverage {} to frequency {}",
                count_total,
                self.coverage_distribution[&count_total]
            );
        }

        // Merge context methylation
        for (context, &(sum_methylation, count)) in &other.context_methylation {
            let entry = self
                .context_methylation
                .entry(*context)
                .or_insert((0.0, 0));
            entry.0 += sum_methylation;
            entry.1 += count;
            trace!(
                "Updated context {:?} to (sum: {}, count: {})",
                context,
                entry.0,
                entry.1
            );
        }

        // Merge strand methylation
        for (strand, &(sum_methylation, count)) in &other.strand_methylation {
            let entry = self
                .strand_methylation
                .entry(*strand)
                .or_insert((0.0, 0));
            entry.0 += sum_methylation;
            entry.1 += count;
            trace!(
                "Updated strand {:?} to (sum: {}, count: {})",
                strand,
                entry.0,
                entry.1
            );
        }
    }

    /// Finalizes methylation statistics by converting accumulated sums to mean
    /// values.
    ///
    /// This should be called after all merging operations are complete and
    /// before retrieving final statistics.
    pub fn finalize_methylation(&mut self) {
        debug!("Finalizing methylation statistics");

        for (context, (sum_methylation, count)) in &mut self.context_methylation
        {
            if *count > 0 {
                *sum_methylation /= *count as f64;
                trace!(
                    "Finalized {:?} context methylation to {}",
                    context,
                    *sum_methylation
                );
            }
            else {
                debug!(
                    "No data_structs for context {:?}, keeping sum at {}",
                    context, *sum_methylation
                );
            }
        }

        for (strand, (sum_methylation, count)) in &mut self.strand_methylation {
            if *count > 0 {
                *sum_methylation /= *count as f64;
                trace!(
                    "Finalized {:?} strand methylation to {}",
                    strand,
                    *sum_methylation
                );
            }
            else {
                debug!(
                    "No data_structs for strand {:?}, keeping sum at {}",
                    strand, *sum_methylation
                );
            }
        }
    }

    /// Computes the total sequencing coverage across all positions.
    ///
    /// # Returns
    ///
    /// The sum of all frequency values in the coverage distribution.
    pub fn total_coverage(&self) -> u32 {
        let total = self
            .coverage_distribution
            .iter()
            .map(|(_, &freq)| freq)
            .filter(|&freq| freq > 0)
            .sum();

        trace!("Calculated total coverage: {}", total);
        total
    }

    /// Computes the genome-wide mean methylation level.
    ///
    /// # Returns
    ///
    /// The mean methylation level, or 0.0 if no coverage data_structs is
    /// available.
    pub fn mean_methylation(&self) -> f64 {
        if self.total_coverage() == 0 {
            debug!(
                "No coverage data_structs, returning 0.0 for mean methylation"
            );
            0.0
        }
        else {
            trace!(
                "Returning calculated mean methylation: {}",
                self.mean_methylation
            );
            self.mean_methylation
        }
    }

    /// Merges multiple MethylationStats instances to create a genome-wide
    /// summary.
    ///
    /// # Arguments
    ///
    /// * `stats_list` - A slice of MethylationStats instances to merge
    ///
    /// # Returns
    ///
    /// A new MethylationStats instance containing the merged data_structs.
    pub fn merge_multiple(stats_list: &[MethylationStats]) -> Self {
        info!("Merging {} MethylationStats instances", stats_list.len());

        let mut genome_stats = MethylationStats::new();
        for (idx, stats) in stats_list.iter().enumerate() {
            debug!(
                "Merging MethylationStats instance {}/{}",
                idx + 1,
                stats_list.len()
            );
            genome_stats.merge(stats);
        }

        debug!("Finalizing merged genome-wide methylation statistics");
        genome_stats.finalize_methylation();

        info!(
            "Genome-wide methylation: mean={:.6}, var={:.6}, total_coverage={}",
            genome_stats.mean_methylation,
            genome_stats.methylation_var,
            genome_stats.total_coverage()
        );

        genome_stats
    }

    /// Generates a detailed text representation of the methylation statistics.
    ///
    /// This method finalizes methylation calculations before generating the
    /// output.
    ///
    /// # Returns
    ///
    /// A formatted string containing detailed methylation statistics.
    pub fn display_long(&mut self) -> Result<String> {
        debug!(
            "Generating detailed text representation of methylation statistics"
        );

        let mut buf = String::new();
        self.finalize_methylation();

        writeln!(&mut buf, "Methylation mean: {:.6}", self.mean_methylation)
            .context("Failed to write mean methylation")?;

        writeln!(
            &mut buf,
            "Methylation variance: {:.6}",
            self.methylation_var
        )
        .context("Failed to write methylation variance")?;

        writeln!(&mut buf).context("Failed to write newline")?;

        writeln!(&mut buf, "Cytosine coverage distribution: ")
            .context("Failed to write coverage distribution header")?;

        writeln!(&mut buf, "coverage\tcount")
            .context("Failed to write coverage distribution column headers")?;

        for (key, value) in self
            .coverage_distribution
            .iter()
            .sorted_by_key(|(k, _)| **k)
        {
            writeln!(&mut buf, "{}\t{}", key, value).context(format!(
                "Failed to write coverage entry for key {}",
                key
            ))?;
        }

        writeln!(&mut buf)
            .context("Failed to write newline after coverage distribution")?;

        writeln!(&mut buf, "Methylation per context: ")
            .context("Failed to write context methylation header")?;

        writeln!(&mut buf, "coverage\tmean\tcount")
            .context("Failed to write context methylation column headers")?;

        for (key, (mean, count)) in self.context_methylation.iter() {
            writeln!(&mut buf, "{}\t{:.6}\t{}", key, mean, count).context(
                format!("Failed to write context entry for {:?}", key),
            )?;
        }

        writeln!(&mut buf)
            .context("Failed to write newline after context methylation")?;

        writeln!(&mut buf, "Methylation per strand")
            .context("Failed to write strand methylation header")?;

        writeln!(&mut buf, "strand\tmean\tcount")
            .context("Failed to write strand methylation column headers")?;

        for (key, (mean, count)) in self.strand_methylation.iter() {
            writeln!(&mut buf, "{}\t{:.6}\t{}", key, mean, count).context(
                format!("Failed to write strand entry for {:?}", key),
            )?;
        }

        trace!("Successfully generated detailed methylation statistics text");
        Ok(buf)
    }

    pub fn coverage_distribution(&self) -> &HashMap<u16, u32> {
        &self.coverage_distribution
    }

    pub fn methylation_var(&self) -> f64 { self.methylation_var }

    pub fn context_methylation(&self) -> &HashMap<Context, (f64, u32)> {
        &self.context_methylation
    }

    pub fn strand_methylation(&self) -> &HashMap<Strand, (f64, u32)> {
        &self.strand_methylation
    }
}

impl Default for MethylationStats {
    fn default() -> Self { Self::new() }
}

#[cfg(test)]
mod tests {
    use hashbrown::HashMap;
    use num::pow::Pow;

    use super::*;

    /// Helper function to create a dummy `MethylationStats`
    fn sample_stats(
        mean: f64,
        variance: f64,
        coverage: Vec<(u16, u32)>,
        context: Vec<(Context, f64, u32)>,
    ) -> MethylationStats {
        let coverage_distribution = coverage
            .into_iter()
            .collect::<HashMap<u16, u32>>();
        let context_methylation = context
            .into_iter()
            .map(|(ctx, sum_m, count)| (ctx, (sum_m, count)))
            .collect::<HashMap<Context, (f64, u32)>>();

        MethylationStats {
            mean_methylation: mean,
            methylation_var: variance,
            coverage_distribution,
            context_methylation,
            strand_methylation: HashMap::new(),
        }
    }

    /// Test creation of an empty `MethylationStats`
    #[test]
    fn test_methylation_stats_new() {
        let stats = MethylationStats::new();
        assert_eq!(stats.mean_methylation, 0.0);
        assert_eq!(stats.methylation_var, 0.0);
        assert!(stats.coverage_distribution.is_empty());
        assert!(stats.context_methylation.is_empty());
    }

    /// Test merging of two `MethylationStats` structures
    #[test]
    fn test_merge_methylation_stats() {
        let mut stats1 = sample_stats(
            0.3,
            0.01,
            vec![(10, 100), (20, 50)],
            vec![(Context::CG, 0.3, 150)],
        );
        let stats2 = sample_stats(0.5, 0.02, vec![(10, 200), (30, 75)], vec![
            (Context::CG, 0.5, 275),
            (Context::CHG, 0.6, 100),
        ]);

        let weight1 = 150.0;
        let weight2 = 275.0;
        let total_weight = weight1 + weight2;
        let mean1 = 0.3;
        let mean2 = 0.5;
        let var1 = 0.01;
        let var2 = 0.02;
        let delta = mean1 - mean2;

        let expected_mean = (weight1 * mean1 + weight2 * mean2) / total_weight;
        let expected_variance = ((weight1 * var1 + weight2 * var2)
            + (weight1 * weight2 / total_weight) * (delta * delta))
            / total_weight;

        stats1.merge(&stats2);
        stats1.finalize_methylation();

        // Check if mean methylation is correctly weighted
        assert!((stats1.mean_methylation - expected_mean).abs() < 1e-6);

        // Check if variance is correctly calculated
        assert!((stats1.methylation_var - expected_variance).abs() < 1e-6);

        // Check merged coverage distribution
        assert_eq!(stats1.coverage_distribution.get(&10), Some(&300));
        assert_eq!(stats1.coverage_distribution.get(&20), Some(&50));
        assert_eq!(stats1.coverage_distribution.get(&30), Some(&75));

        // Check finalized context methylation (ensuring it's averaged)
        assert!(
            (stats1
                .context_methylation
                .get(&Context::CG)
                .unwrap()
                .0
                - (0.3 + 0.5) / 425.0)
                .abs()
                < 1e-6
        );
        assert!(
            (stats1
                .context_methylation
                .get(&Context::CHG)
                .unwrap()
                .0
                - (0.6 / 100.0))
                .abs()
                < 1e-6
        );
    }

    /// Test genome-wide mean methylation calculation
    #[test]
    fn test_genome_wide_mean_methylation() {
        let stats =
            sample_stats(0.4, 0.02, vec![(10, 100), (20, 200)], vec![(
                Context::CG,
                0.4,
                300,
            )]);
        assert!((stats.mean_methylation() - 0.4).abs() < 1e-6);
    }

    /// Test merging multiple methylation statistics
    #[test]
    fn test_merge_multiple_methylation_stats() {
        let stats1 = sample_stats(0.2, 0.01, vec![(10, 100)], vec![(
            Context::CG,
            0.2,
            100,
        )]);
        let stats2 = sample_stats(0.4, 0.02, vec![(20, 200)], vec![(
            Context::CG,
            0.4,
            200,
        )]);
        let stats3 = sample_stats(0.6, 0.03, vec![(30, 300)], vec![(
            Context::CHG,
            0.6,
            300,
        )]);

        let merged =
            MethylationStats::merge_multiple(&vec![stats1, stats2, stats3]);

        let weight1 = 100.0;
        let weight2 = 200.0;
        let weight3 = 300.0;
        let total_weight = weight1 + weight2 + weight3;
        let mean1 = 0.2;
        let mean2 = 0.4;
        let mean3 = 0.6;
        let var1 = 0.01;
        let var2 = 0.02;
        let var3 = 0.03;

        let expected_mean =
            (weight1 * mean1 + weight2 * mean2 + weight3 * mean3)
                / total_weight;
        let expected_variance =
            ((weight1 * var1 + weight2 * var2 + weight3 * var3)
                + (weight1 * weight2 / total_weight) * (mean1 - mean2).pow(2)
                + (weight1 * weight3 / total_weight) * (mean1 - mean3).pow(2)
                + (weight2 * weight3 / total_weight) * (mean2 - mean3).pow(2))
                / total_weight;

        assert!((merged.mean_methylation - expected_mean).abs() < 1e-6);
        assert!(
            (merged.methylation_var as f64 - expected_variance as f64).abs()
                < 1e-6
        );
        assert_eq!(merged.coverage_distribution.get(&10), Some(&100));
        assert_eq!(merged.coverage_distribution.get(&20), Some(&200));
        assert_eq!(merged.coverage_distribution.get(&30), Some(&300));

        // Finalize and check context-specific methylation
        assert!(
            (merged
                .context_methylation
                .get(&Context::CG)
                .unwrap()
                .0
                - (0.2 + 0.4) / 300.0)
                .abs()
                < 1e-6
        );
        assert!(
            (merged
                .context_methylation
                .get(&Context::CHG)
                .unwrap()
                .0
                - (0.6 / 300.0))
                .abs()
                < 1e-6
        );
    }

    /// Test finalization of context methylation levels
    #[test]
    fn test_finalize_context_methylation() {
        let mut stats = sample_stats(0.5, 0.01, vec![], vec![
            (Context::CG, 1.5, 3),
            (Context::CHH, 0.9, 3),
        ]);
        stats.finalize_methylation();

        assert!(
            (stats
                .context_methylation
                .get(&Context::CG)
                .unwrap()
                .0
                - 0.5)
                .abs()
                < 1e-6
        );
        assert!(
            (stats
                .context_methylation
                .get(&Context::CHH)
                .unwrap()
                .0
                - 0.3)
                .abs()
                < 1e-6
        );
    }

    /// Test serialization and deserialization
    #[test]
    fn test_serialization_deserialization() {
        let stats = sample_stats(
            0.35,
            0.015,
            vec![(10, 150), (20, 250)],
            vec![(Context::CG, 0.35, 400)],
        );
        let json = serde_json::to_string(&stats).expect("Serialization failed");
        let deserialized: MethylationStats =
            serde_json::from_str(&json).expect("Deserialization failed");

        assert_eq!(stats.mean_methylation, deserialized.mean_methylation);
        assert_eq!(stats.methylation_var, deserialized.methylation_var);
        assert_eq!(
            stats.coverage_distribution,
            deserialized.coverage_distribution
        );
        assert_eq!(stats.context_methylation, deserialized.context_methylation);
    }
}
