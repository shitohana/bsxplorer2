use std::collections::BTreeMap;
use std::fmt::Write;

use hashbrown::HashMap;
use serde::{
    Deserialize,
    Serialize,
    Serializer,
};

use super::typedef::{
    CountType,
    DensityType,
};
use crate::data_structs::enums::{
    Context,
    Strand,
};

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
