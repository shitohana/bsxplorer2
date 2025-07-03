#![allow(unused)]
use std::ops::{Add, AddAssign};

use hashbrown::HashMap;
use serde::{Deserialize, Serialize};

use crate::data_structs::typedef::{CountType, DensityType};
use crate::data_structs::MethylationStats;
use crate::getter_fn;
use crate::prelude::{Context, Strand};

#[derive(Clone, Debug, Copy, Serialize, Deserialize)]
struct MethAgg {
    sum: DensityType,
    count: DensityType
}

impl AddAssign for MethAgg {
    fn add_assign(&mut self, rhs: Self) {
        self.sum += rhs.sum;
        self.count += rhs.count;
    }
}

impl Add for MethAgg {
    type Output = MethAgg;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.sum += rhs.sum;
        self.count += rhs.count;
        self
    }
}

impl Default for MethAgg {
    fn default() -> Self {
        Self::new()
    }
}

impl MethAgg {
    getter_fn!(*sum, DensityType);
    getter_fn!(*count, DensityType);

    pub fn new() -> Self {
        Self { sum: 0.0, count: 0.0 }
    }
    
    pub fn add_density(&mut self, density: DensityType) {
        self.count += 1.0 as DensityType;
        self.sum += density;
    }

    pub fn finalize(&self) -> DensityType {
        self.sum / self.count
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
struct RegionMethAgg {
    context: HashMap<Context, MethAgg>,
    strand: HashMap<Strand, MethAgg>,
    coverage: HashMap<CountType, usize>
}

impl RegionMethAgg {
    getter_fn!(coverage, HashMap<CountType, usize>);
    getter_fn!(context, HashMap<Context, MethAgg>);
    getter_fn!(strand, HashMap<Strand, MethAgg>);
    
    fn new() -> Self {
        Self::default()
    }

    fn is_empty(&self) -> bool {
        self.context.is_empty() && self.strand.is_empty()
    }

    fn add_cytosine(&mut self, sum: DensityType, count: Option<CountType>, context: Context, strand: Strand) {
        if let Some(coverage) = count {
            self.coverage.entry(coverage).or_insert(0).add_assign(1);
        }
        let density = sum / count.unwrap_or(1) as DensityType;
        self.strand.entry(strand).or_default().add_density(density);
        self.context.entry(context).or_default().add_density(density);
    }

    fn finalize(self) -> MethylationStats {
        MethylationStats::new()
    }
}
