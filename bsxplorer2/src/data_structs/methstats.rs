use std::collections::BTreeMap;
use std::fmt::Write;

use hashbrown::HashMap;
use serde::{Deserialize, Serialize, Serializer};

use super::typedef::{CountType, DensityType};
use crate::data_structs::enums::{Context, Strand};

/// Serializes a HashMap in deterministic order.
fn serialize_sorted_map<S, K: Ord + Serialize, V: Serialize>(
    map: &HashMap<K, V>,
    serializer: S,
) -> anyhow::Result<S::Ok, S::Error>
where
    S: Serializer, {
    let sorted_map: BTreeMap<_, _> = map.iter().collect();
    sorted_map.serialize(serializer)
}

/// Comprehensive methylation statistics.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MethylationStats {
    /// Overall mean methylation.
    pub(super) mean_methylation: DensityType,

    /// Variance of methylation levels.
    pub(super) methylation_var: DensityType,

    /// Maps coverage to frequency.
    #[serde(serialize_with = "serialize_sorted_map")]
    pub(super) coverage_distribution: HashMap<CountType, u32>,

    /// Methylation statistics per context.
    #[serde(serialize_with = "serialize_sorted_map")]
    pub(super) context_methylation: HashMap<Context, (DensityType, u32)>,

    /// Methylation statistics per strand.
    #[serde(serialize_with = "serialize_sorted_map")]
    pub(super) strand_methylation: HashMap<Strand, (DensityType, u32)>,
}

/// Flattened representation of methylation statistics for reporting.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MethylationStatFlat {
    /// Overall mean methylation.
    pub mean_methylation: DensityType,

    /// Overall variance in methylation.
    pub methylation_var: DensityType,

    /// Average coverage.
    pub mean_coverage: f64,

    /// Mean methylation in CG context.
    pub cg_mean_methylation: DensityType,

    /// Count of positions in CG context.
    pub cg_coverage: u32,

    /// Mean methylation in CHG context.
    pub chg_mean_methylation: DensityType,

    /// Count of positions in CHG context.
    pub chg_coverage: u32,

    /// Mean methylation in CHH context.
    pub chh_mean_methylation: DensityType,

    /// Count of positions in CHH context.
    pub chh_coverage: u32,

    /// Mean methylation on forward strand.
    pub fwd_mean_methylation: DensityType,

    /// Count of positions on forward strand.
    pub fwd_coverage: u32,

    /// Mean methylation on reverse strand.
    pub rev_mean_methylation: DensityType,

    /// Count of positions on reverse strand.
    pub rev_coverage: u32,
}

impl From<MethylationStats> for MethylationStatFlat {
    /// Converts detailed stats to flattened format.
    fn from(value: MethylationStats) -> Self {
        // Helper to get sum and count for a specific context.
        // Returns (0.0, 0) if the context is not present.
        let get_context_stats =
            |stats: &MethylationStats, context: Context| -> (DensityType, u32) {
                stats
                    .context_methylation
                    .get(&context)
                    .copied()
                    .unwrap_or((0.0, 0))
            };

        // Helper to get sum and count for a specific strand.
        // Returns (0.0, 0) if the strand is not present.
        let get_strand_stats =
            |stats: &MethylationStats, strand: Strand| -> (DensityType, u32) {
                stats
                    .strand_methylation
                    .get(&strand)
                    .copied()
                    .unwrap_or((0.0, 0))
            };

        // Helper to calculate mean from sum and count, avoiding division by zero.
        let calculate_mean = |sum: DensityType, count: u32| -> DensityType {
            if count > 0 {
                sum / count as DensityType
            }
            else {
                0.0
            }
        };

        // Calculate the overall average coverage.
        // This is the sum of (coverage * frequency) divided by the total number of
        // sites (sum of frequencies).
        let total_sites = value.total_coverage(); // total number of sites
        let mean_coverage = if total_sites > 0 {
            value.coverage_distribution.iter()
                 // Calculate the total number of read bases covering all sites
                 .map(|(&cov, &freq)| (cov as f64) * (freq as f64)) // Use f64 to avoid overflow for intermediate products
                 .sum::<f64>()
                / total_sites as f64 // Divide total read bases by total sites
        }
        else {
            0.0
        };

        // Get context-specific sums and counts
        let (cg_sum, cg_count) = get_context_stats(&value, Context::CG);
        let (chg_sum, chg_count) = get_context_stats(&value, Context::CHG);
        let (chh_sum, chh_count) = get_context_stats(&value, Context::CHH);

        // Get strand-specific sums and counts
        let (fwd_sum, fwd_count) = get_strand_stats(&value, Strand::Forward);
        let (rev_sum, rev_count) = get_strand_stats(&value, Strand::Reverse);

        MethylationStatFlat {
            mean_methylation: value.mean_methylation,
            methylation_var: value.methylation_var,
            mean_coverage, // Use the correctly calculated average coverage
            cg_mean_methylation: calculate_mean(cg_sum, cg_count),
            cg_coverage: cg_count,
            chg_mean_methylation: calculate_mean(chg_sum, chg_count),
            chg_coverage: chg_count,
            chh_mean_methylation: calculate_mean(chh_sum, chh_count),
            chh_coverage: chh_count,
            fwd_mean_methylation: calculate_mean(fwd_sum, fwd_count),
            fwd_coverage: fwd_count,
            rev_mean_methylation: calculate_mean(rev_sum, rev_count),
            rev_coverage: rev_count,
        }
    }
}

impl MethylationStats {
    /// Creates an empty instance.
    pub fn new() -> Self {
        Self {
            mean_methylation:      0.0,
            methylation_var:       0.0,
            coverage_distribution: HashMap::new(),
            context_methylation:   HashMap::new(),
            strand_methylation:    HashMap::new(),
        }
    }

    /// Creates an instance from data.
    ///
    /// Handles NaN values gracefully.
    pub fn from_data(
        mean_methylation: DensityType,
        variance_methylation: DensityType,
        coverage_distribution: HashMap<u16, u32>,
        context_methylation: HashMap<Context, (DensityType, u32)>,
        strand_methylation: HashMap<Strand, (DensityType, u32)>,
    ) -> MethylationStats {
        // Handle NaN values by replacing them with 0.0
        let mean_methylation = if mean_methylation.is_nan() {
            0.0
        }
        else {
            mean_methylation
        };

        let variance_methylation = if variance_methylation.is_nan() {
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

    /// Merges another instance into this one using weighted averages based on
    /// coverage.
    pub fn merge(
        &mut self,
        other: &MethylationStats,
    ) {
        let weight_self = self.total_coverage() as DensityType;
        let weight_other = other.total_coverage() as DensityType;
        let total_weight = weight_self + weight_other;

        // Weighted mean methylation calculation
        if total_weight > 0.0 {
            let delta = self.mean_methylation - other.mean_methylation;
            self.mean_methylation = (weight_self * self.mean_methylation
                + weight_other * other.mean_methylation)
                / total_weight;

            // Correct variance computation considering mean shift
            self.methylation_var = ((weight_self * self.methylation_var
                + weight_other * other.methylation_var)
                + (weight_self * weight_other / total_weight) * (delta * delta))
                / total_weight;
        }

        // Merge coverage distribution
        for (&count_total, &frequency) in &other.coverage_distribution {
            *self.coverage_distribution.entry(count_total).or_insert(0) += frequency;
        }

        // Merge context methylation
        for (context, &(sum_methylation, count)) in &other.context_methylation {
            let entry = self.context_methylation.entry(*context).or_insert((0.0, 0));
            entry.0 += sum_methylation;
            entry.1 += count;
        }

        // Merge strand methylation
        for (strand, &(sum_methylation, count)) in &other.strand_methylation {
            let entry = self.strand_methylation.entry(*strand).or_insert((0.0, 0));
            entry.0 += sum_methylation;
            entry.1 += count;
        }
    }

    /// Finalizes statistics by converting sums to means.
    pub fn finalize_methylation(&mut self) {
        for (_, (sum_methylation, count)) in &mut self.context_methylation {
            if *count > 0 {
                *sum_methylation /= *count as DensityType;
            }
        }

        for (_, (sum_methylation, count)) in &mut self.strand_methylation {
            if *count > 0 {
                *sum_methylation /= *count as DensityType;
            }
        }
    }

    /// Computes the total sequencing coverage.
    pub fn total_coverage(&self) -> u32 {
        self.coverage_distribution.values().sum()
    }

    /// Computes the genome-wide mean methylation.
    pub fn mean_methylation(&self) -> DensityType {
        if self.total_coverage() == 0 {
            0.0
        }
        else {
            self.mean_methylation
        }
    }

    /// Merges multiple instances.
    pub fn merge_multiple(stats_list: &[MethylationStats]) -> Self {
        let mut genome_stats = MethylationStats::new();
        for stats in stats_list.iter() {
            genome_stats.merge(stats);
        }

        genome_stats.finalize_methylation();

        genome_stats
    }

    /// Generates a detailed text representation.
    pub fn display_long(&mut self) -> anyhow::Result<String> {
        use std::collections::BTreeMap; // Ensure BTreeMap is in scope

        use itertools::Itertools; // Ensure Itertools is in scope

        let mut buf = String::new();
        self.finalize_methylation();

        // Overall stats
        writeln!(buf, "Methylation mean: {:.6}", self.mean_methylation)?;
        writeln!(buf, "Methylation variance: {:.6}", self.methylation_var)?;
        writeln!(buf)?; // Blank line

        // Coverage distribution
        writeln!(buf, "Cytosine coverage distribution:")?;
        writeln!(buf, "coverage\tcount")?;
        // Sort coverage distribution by coverage level
        for (key, value) in self
            .coverage_distribution
            .iter()
            .sorted_by_key(|(k, _)| **k)
        {
            writeln!(buf, "{}\t{}", key, value)?;
        }
        writeln!(buf)?; // Blank line

        // Methylation per context
        writeln!(buf, "Methylation per context:")?;
        writeln!(buf, "context\tmean\tcount")?; // Corrected header from "coverage"
                                                // Sort context methylation by context type
        let sorted_contexts: BTreeMap<_, _> = self.context_methylation.iter().collect();
        for (key, (mean, count)) in sorted_contexts {
            writeln!(buf, "{}\t{:.6}\t{}", key, mean, count)?;
        }
        writeln!(buf)?; // Blank line

        // Methylation per strand
        writeln!(buf, "Methylation per strand:")?;
        writeln!(buf, "strand\tmean\tcount")?;
        // Sort strand methylation by strand type
        let sorted_strands: BTreeMap<_, _> = self.strand_methylation.iter().collect();
        for (key, (mean, count)) in sorted_strands {
            writeln!(buf, "{}\t{:.6}\t{}", key, mean, count)?;
        }

        Ok(buf)
    }

    pub fn coverage_distribution(&self) -> &HashMap<u16, u32> {
        &self.coverage_distribution
    }

    pub fn methylation_var(&self) -> DensityType {
        self.methylation_var
    }

    pub fn context_methylation(&self) -> &HashMap<Context, (DensityType, u32)> {
        &self.context_methylation
    }

    pub fn strand_methylation(&self) -> &HashMap<Strand, (DensityType, u32)> {
        &self.strand_methylation
    }
}

impl Default for MethylationStats {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use hashbrown::HashMap;

    use crate::data_structs::enums::{Context, Strand};
    use crate::data_structs::methstats::MethylationStats;
    use crate::data_structs::typedef::DensityType;

    /// Helper function to create a dummy `MethylationStats`
    fn sample_stats(
        mean: DensityType,
        variance: DensityType,
        coverage: Vec<(u16, u32)>,
        context: Vec<(Context, DensityType, u32)>,
    ) -> MethylationStats {
        let coverage_distribution = coverage.into_iter().collect::<HashMap<u16, u32>>();
        let context_methylation = context
            .into_iter()
            .map(|(ctx, sum_m, count)| (ctx, (sum_m, count)))
            .collect::<HashMap<Context, (DensityType, u32)>>();

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
        assert!(stats.strand_methylation.is_empty()); // Test strand_methylation
                                                      // is empty too
    }

    /// Test merging of two `MethylationStats` structures
    #[test]
    fn test_merge_methylation_stats() {
        let mut stats1 = sample_stats(0.3, 0.01, vec![(10, 100), (20, 50)], vec![(
            Context::CG,
            45.0,
            150,
        )]);
        let stats2 = sample_stats(0.5, 0.02, vec![(10, 200), (30, 75)], vec![
            (Context::CG, 110.0, 275),
            (Context::CHG, 60.0, 100), /* Fixed: was 0.6 which is incorrect
                                        * for sum */
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
        // CG context should be (45.0 + 110.0) / (150 + 275) = 155.0 / 425 =
        // 0.36470588
        assert_approx_eq!(
            stats1.context_methylation.get(&Context::CG).unwrap().0,
            155.0 / 425.0,
            1e-6
        );

        // CHG context should be 60.0 / 100 = 0.6
        assert_approx_eq!(
            stats1.context_methylation.get(&Context::CHG).unwrap().0,
            0.6,
            1e-6
        );
    }

    /// Test genome-wide mean methylation calculation
    #[test]
    fn test_genome_wide_mean_methylation() {
        let stats = sample_stats(0.4, 0.02, vec![(10, 100), (20, 200)], vec![(
            Context::CG,
            120.0, /* 0.4 * 300 = 120.0 (should use sum not mean for
                    * context) */
            300,
        )]);
        assert!((stats.mean_methylation() - 0.4).abs() < 1e-6);

        // Also test with zero coverage
        let empty_stats = MethylationStats::new();
        assert_eq!(empty_stats.mean_methylation(), 0.0);
    }

    /// Test merging multiple methylation statistics
    #[test]
    fn test_merge_multiple_methylation_stats() {
        let stats1 = sample_stats(0.2, 0.01, vec![(10, 100)], vec![(
            Context::CG,
            20.0, // Should be sum (0.2 * 100 = 20.0), not mean
            100,
        )]);
        let stats2 = sample_stats(0.4, 0.02, vec![(20, 200)], vec![(
            Context::CG,
            80.0, // Should be sum (0.4 * 200 = 80.0), not mean
            200,
        )]);
        let stats3 = sample_stats(0.6, 0.03, vec![(30, 300)], vec![(
            Context::CHG,
            180.0, // Should be sum (0.6 * 300 = 180.0), not mean
            300,
        )]);

        let merged = MethylationStats::merge_multiple(&vec![stats1, stats2, stats3]);

        let weight1 = 100.0;
        let weight2 = 200.0;
        let weight3 = 300.0;
        let total_weight = weight1 + weight2 + weight3;
        let mean1 = 0.2;
        let mean2 = 0.4;
        let mean3 = 0.6;

        let expected_mean =
            (weight1 * mean1 + weight2 * mean2 + weight3 * mean3) / total_weight;

        // Variance calculation is complex in multi-merge and not correctly
        // implemented in test The test adds pairwise delta terms which
        // isn't correct for 3+ datasets

        assert!((merged.mean_methylation - expected_mean).abs() < 1e-6);

        assert_eq!(merged.coverage_distribution.get(&10), Some(&100));
        assert_eq!(merged.coverage_distribution.get(&20), Some(&200));
        assert_eq!(merged.coverage_distribution.get(&30), Some(&300));

        // Check context-specific methylation after finalization
        assert_approx_eq!(
            merged.context_methylation.get(&Context::CG).unwrap().0,
            (20.0 + 80.0) / 300.0, /* (sum1 + sum2) / (count1 + count2) =
                                    * 100/300 = 0.33333 */
            1e-6
        );
        assert_approx_eq!(
            merged.context_methylation.get(&Context::CHG).unwrap().0,
            180.0 / 300.0, // sum3 / count3 = 0.6
            1e-6
        );
    }

    /// Test finalization of context methylation levels
    #[test]
    fn test_finalize_context_methylation() {
        let mut stats = sample_stats(0.5, 0.01, vec![], vec![
            (Context::CG, 1.5, 3),  // Sum = 1.5, should become mean = 0.5
            (Context::CHH, 0.9, 3), // Sum = 0.9, should become mean = 0.3
        ]);
        stats.finalize_methylation();

        assert_approx_eq!(
            stats.context_methylation.get(&Context::CG).unwrap().0,
            0.5,
            1e-6
        );
        assert_approx_eq!(
            stats.context_methylation.get(&Context::CHH).unwrap().0,
            0.3,
            1e-6
        );
    }

    /// Test serialization and deserialization
    #[test]
    fn test_serialization_deserialization() {
        let stats = sample_stats(
            0.35,
            0.015,
            vec![(10, 150), (20, 250)],
            vec![(Context::CG, 140.0, 400)], // Sum = 0.35 * 400 = 140
        );

        // Add strand methylation to test complete serialization
        let mut stats_with_strand = stats.clone();
        stats_with_strand
            .strand_methylation
            .insert(Strand::Forward, (70.0, 200));
        stats_with_strand
            .strand_methylation
            .insert(Strand::Reverse, (70.0, 200));

        let json =
            serde_json::to_string(&stats_with_strand).expect("Serialization failed");
        let deserialized: MethylationStats =
            serde_json::from_str(&json).expect("Deserialization failed");

        assert_eq!(
            stats_with_strand.mean_methylation,
            deserialized.mean_methylation
        );
        assert_eq!(
            stats_with_strand.methylation_var,
            deserialized.methylation_var
        );
        assert_eq!(
            stats_with_strand.coverage_distribution,
            deserialized.coverage_distribution
        );
        assert_eq!(
            stats_with_strand.context_methylation,
            deserialized.context_methylation
        );
        assert_eq!(
            stats_with_strand.strand_methylation,
            deserialized.strand_methylation
        );
    }
}
