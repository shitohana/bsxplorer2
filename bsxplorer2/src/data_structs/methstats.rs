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

impl MethylationStats {
    /// Creates an empty instance.
    pub fn new() -> Self {
        Self {
            mean_methylation:      0.0,
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

        Self {
            mean_methylation,
            context_methylation,
            coverage_distribution,
            strand_methylation,
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

    /// Generates a detailed text representation.
    pub fn display_long(&mut self) -> anyhow::Result<String> {
        use std::collections::BTreeMap; // Ensure BTreeMap is in scope

        use itertools::Itertools; // Ensure Itertools is in scope

        let mut buf = String::new();
        self.finalize_methylation();

        // Overall stats
        writeln!(buf, "Methylation mean: {:.6}", self.mean_methylation)?;
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
