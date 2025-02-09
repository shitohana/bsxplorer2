use crate::utils::types::{Context, IPCEncodedEnum, Strand};
use itertools::Itertools;
use serde::{Deserialize, Serialize, Serializer};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write;

fn serialize_sorted_map<S, K: Ord + Serialize, V: Serialize>(
    map: &HashMap<K, V>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let sorted_map: BTreeMap<_, _> = map.iter().collect();
    sorted_map.serialize(serializer)
}
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MethylationStats {
    mean_methylation: f64,
    methylation_var: f64,
    #[serde(serialize_with = "serialize_sorted_map")]
    coverage_distribution: HashMap<u16, u32>, // (count_total -> frequency)
    #[serde(serialize_with = "serialize_sorted_map")]
    context_methylation: HashMap<Context, (f64, u32)>, // (context -> (sum_methylation, count))
    #[serde(serialize_with = "serialize_sorted_map")]
    strand_methylation: HashMap<Strand, (f64, u32)>, // (strand -> (sum_methylation, count))
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MethylationStatFlat {
    pub mean_methylation: f64,
    pub methylation_var: f64,
    pub mean_coverage: f64,
    pub cg_mean_methylation: f64,
    pub cg_coverage: u32,
    pub chg_mean_methylation: f64,
    pub chg_coverage: u32,
    pub chh_mean_methylation: f64,
    pub chh_coverage: u32,
    pub fwd_mean_methylation: f64,
    pub fwd_coverage: u32,
    pub rev_mean_methylation: f64,
    pub rev_coverage: u32,
}

impl From<MethylationStats> for MethylationStatFlat {
    fn from(value: MethylationStats) -> Self {
        MethylationStatFlat {
            mean_methylation: value.mean_methylation,
            methylation_var: value.methylation_var,
            mean_coverage: value.total_coverage() as f64 / value.coverage_distribution.len() as f64,
            cg_mean_methylation: value
                .context_methylation
                .get(&Context::CG)
                .map_or(0.0, |(sum_m, count)| sum_m / *count as f64),
            cg_coverage: value
                .context_methylation
                .get(&Context::CG)
                .map_or(0, |(_, count)| *count),
            chg_mean_methylation: value
                .context_methylation
                .get(&Context::CHG)
                .map_or(0.0, |(sum_m, count)| sum_m / *count as f64),
            chg_coverage: value
                .context_methylation
                .get(&Context::CHG)
                .map_or(0, |(_, count)| *count),
            chh_mean_methylation: value
                .context_methylation
                .get(&Context::CHH)
                .map_or(0.0, |(sum_m, count)| sum_m / *count as f64),
            chh_coverage: value
                .context_methylation
                .get(&Context::CHH)
                .map_or(0, |(_, count)| *count),
            fwd_mean_methylation: value
                .strand_methylation
                .get(&Strand::Forward)
                .map_or(0.0, |(sum_m, count)| sum_m / *count as f64),
            fwd_coverage: value
                .strand_methylation
                .get(&Strand::Forward)
                .map_or(0, |(_, count)| *count),
            rev_mean_methylation: value
                .strand_methylation
                .get(&Strand::Forward)
                .map_or(0.0, |(sum_m, count)| sum_m / *count as f64),
            rev_coverage: value
                .strand_methylation
                .get(&Strand::Reverse)
                .map_or(0, |(_, count)| *count),
        }
    }
}

impl MethylationStats {
    /// Create an empty methylation statistics instance
    pub fn new() -> Self {
        Self {
            mean_methylation: 0.0,
            methylation_var: 0.0,
            coverage_distribution: HashMap::new(),
            context_methylation: HashMap::new(),
            strand_methylation: HashMap::new(),
        }
    }

    pub fn from_data(
        mean_methylation: f64,
        variance_methylation: f64,
        coverage_distribution: HashMap<u16, u32>,
        context_methylation: HashMap<Context, (f64, u32)>,
        strand_methylation: HashMap<Strand, (f64, u32)>,
    ) -> MethylationStats {
        let mean_methylation = if mean_methylation.is_nan() {
            0.0
        } else {
            mean_methylation
        };
        let variance_methylation = if variance_methylation.is_nan() {
            0.0
        } else {
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

    pub fn merge(&mut self, other: &MethylationStats) {
        let weight_self = self.total_coverage() as f64;
        let weight_other = other.total_coverage() as f64;
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
            let entry = self
                .context_methylation
                .entry(context.clone())
                .or_insert((0.0, 0));
            entry.0 += sum_methylation;
            entry.1 += count;
        }

        // Merge strand methylation
        for (strand, &(sum_methylation, count)) in &other.strand_methylation {
            let entry = self
                .strand_methylation
                .entry(strand.clone())
                .or_insert((0.0, 0));
            entry.0 += sum_methylation;
            entry.1 += count;
        }
    }

    /// Compute final per-context methylation levels (mean methylation per context)
    pub fn finalize_methylation(&mut self) {
        for (_, (sum_methylation, count)) in &mut self.context_methylation {
            if *count > 0 {
                *sum_methylation /= *count as f64; // Convert sum to actual mean methylation per context
            }
        }
        for (_, (sum_methylation, count)) in &mut self.strand_methylation {
            if *count > 0 {
                *sum_methylation /= *count as f64; // Convert sum to actual mean methylation per strand
            }
        }
    }

    /// Compute total sequencing coverage
    pub fn total_coverage(&self) -> u32 {
        self.coverage_distribution
            .iter()
            .map(|(_, &freq)| freq)
            .filter(|&freq| freq > 0)
            .sum()
    }

    /// Compute genome-wide mean methylation
    pub fn mean_methylation(&self) -> f64 {
        if self.total_coverage() == 0 {
            0.0
        } else {
            self.mean_methylation
        }
    }

    /// Merge multiple batches to create a genome-wide summary
    pub fn merge_multiple(stats_list: &[MethylationStats]) -> Self {
        let mut genome_stats = MethylationStats::new();
        for stats in stats_list {
            genome_stats.merge(stats);
        }
        genome_stats.finalize_methylation();
        genome_stats
    }

    pub fn display_long(&mut self) -> String {
        let mut buf = String::new();
        self.finalize_methylation();

        writeln!(&mut buf, "Methylation mean: {:.6}", self.mean_methylation).unwrap();
        writeln!(
            &mut buf,
            "Methylation variance: {:.6}",
            self.methylation_var
        )
        .unwrap();
        writeln!(&mut buf).unwrap();
        writeln!(&mut buf, "Cytosine coverage distribution: ").unwrap();
        writeln!(&mut buf, "coverage\tcount").unwrap();
        for (key, value) in self
            .coverage_distribution
            .iter()
            .sorted_by_key(|(k, _)| **k)
        {
            writeln!(&mut buf, "{}\t{}", key, value).unwrap();
        }
        writeln!(&mut buf).unwrap();
        writeln!(&mut buf, "Methylation per context: ").unwrap();
        writeln!(&mut buf, "coverage\tmean\tcount").unwrap();
        for (key, (mean, count)) in self.context_methylation.iter() {
            writeln!(&mut buf, "{}\t{:.6}\t{}", key.to_string(), mean, count).unwrap();
        }
        writeln!(&mut buf).unwrap();
        writeln!(&mut buf, "Methylation per strand").unwrap();
        writeln!(&mut buf, "strand\tmean\tcount").unwrap();
        for (key, (mean, count)) in self.strand_methylation.iter() {
            writeln!(&mut buf, "{}\t{:.6}\t{}", key.to_string(), mean, count).unwrap();
        }
        buf
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::pow::Pow;
    use std::collections::HashMap;

    /// Helper function to create a dummy `MethylationStats`
    fn sample_stats(
        mean: f64,
        variance: f64,
        coverage: Vec<(u16, u32)>,
        context: Vec<(Context, f64, u32)>,
    ) -> MethylationStats {
        let coverage_distribution = coverage.into_iter().collect::<HashMap<u16, u32>>();
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
        let stats2 = sample_stats(
            0.5,
            0.02,
            vec![(10, 200), (30, 75)],
            vec![(Context::CG, 0.5, 275), (Context::CHG, 0.6, 100)],
        );

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
            (stats1.context_methylation.get(&Context::CG).unwrap().0 - (0.3 + 0.5) / 425.0).abs()
                < 1e-6
        );
        assert!(
            (stats1.context_methylation.get(&Context::CHG).unwrap().0 - (0.6 / 100.0)).abs() < 1e-6
        );
    }

    /// Test genome-wide mean methylation calculation
    #[test]
    fn test_genome_wide_mean_methylation() {
        let stats = sample_stats(
            0.4,
            0.02,
            vec![(10, 100), (20, 200)],
            vec![(Context::CG, 0.4, 300)],
        );
        assert!((stats.mean_methylation() - 0.4).abs() < 1e-6);
    }

    /// Test merging multiple methylation statistics
    #[test]
    fn test_merge_multiple_methylation_stats() {
        let stats1 = sample_stats(0.2, 0.01, vec![(10, 100)], vec![(Context::CG, 0.2, 100)]);
        let stats2 = sample_stats(0.4, 0.02, vec![(20, 200)], vec![(Context::CG, 0.4, 200)]);
        let stats3 = sample_stats(0.6, 0.03, vec![(30, 300)], vec![(Context::CHG, 0.6, 300)]);

        let merged = MethylationStats::merge_multiple(&vec![stats1, stats2, stats3]);

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

        let expected_mean = (weight1 * mean1 + weight2 * mean2 + weight3 * mean3) / total_weight;
        let expected_variance = ((weight1 * var1 + weight2 * var2 + weight3 * var3)
            + (weight1 * weight2 / total_weight) * (mean1 - mean2).pow(2)
            + (weight1 * weight3 / total_weight) * (mean1 - mean3).pow(2)
            + (weight2 * weight3 / total_weight) * (mean2 - mean3).pow(2))
            / total_weight;

        assert!((merged.mean_methylation - expected_mean).abs() < 1e-6);
        assert!(((merged.methylation_var - expected_variance) as f64).abs() < 1e-6);
        assert_eq!(merged.coverage_distribution.get(&10), Some(&100));
        assert_eq!(merged.coverage_distribution.get(&20), Some(&200));
        assert_eq!(merged.coverage_distribution.get(&30), Some(&300));

        // Finalize and check context-specific methylation
        assert!(
            (merged.context_methylation.get(&Context::CG).unwrap().0 - (0.2 + 0.4) / 300.0).abs()
                < 1e-6
        );
        assert!(
            (merged.context_methylation.get(&Context::CHG).unwrap().0 - (0.6 / 300.0)).abs() < 1e-6
        );
    }

    /// Test finalization of context methylation levels
    #[test]
    fn test_finalize_context_methylation() {
        let mut stats = sample_stats(
            0.5,
            0.01,
            vec![],
            vec![(Context::CG, 1.5, 3), (Context::CHH, 0.9, 3)],
        );
        stats.finalize_methylation();

        assert!((stats.context_methylation.get(&Context::CG).unwrap().0 - 0.5).abs() < 1e-6);
        assert!((stats.context_methylation.get(&Context::CHH).unwrap().0 - 0.3).abs() < 1e-6);
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
