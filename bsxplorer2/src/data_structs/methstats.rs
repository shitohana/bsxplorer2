#![allow(unused)]
use std::ops::{
    Add,
    AddAssign,
};

use hashbrown::HashMap;
use itertools::{
    izip,
    Itertools,
};
use serde::{
    Deserialize,
    Serialize,
};

use crate::data_structs::typedef::{
    CountType,
    DensityType,
};
use crate::getter_fn;
use crate::prelude::{
    Context,
    Strand,
};

#[derive(Clone, Debug, Copy, Serialize, Deserialize)]
pub struct MethAgg {
    sum:   DensityType,
    count: DensityType,
}

impl Default for MethAgg {
    fn default() -> Self {
        Self::new()
    }
}

impl AddAssign for MethAgg {
    fn add_assign(
        &mut self,
        rhs: Self,
    ) {
        self.sum += rhs.sum;
        self.count += rhs.count;
    }
}

impl Add for MethAgg {
    type Output = MethAgg;

    fn add(
        mut self,
        rhs: Self,
    ) -> Self::Output {
        self.sum += rhs.sum;
        self.count += rhs.count;
        self
    }
}

impl From<(DensityType, DensityType)> for MethAgg {
    fn from(value: (DensityType, DensityType)) -> Self {
        Self {
            sum:   value.0,
            count: value.1,
        }
    }
}

impl MethAgg {
    getter_fn!(*sum, DensityType);

    getter_fn!(*count, DensityType);

    pub fn new() -> Self {
        Self {
            sum:   0.0,
            count: 0.0,
        }
    }

    pub fn add_density(
        &mut self,
        density: DensityType,
    ) {
        self.count += 1.0 as DensityType;
        self.sum += density;
    }

    pub fn finalize(&self) -> DensityType {
        self.sum / self.count
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct RegionMethAgg {
    context:  HashMap<Context, MethAgg>,
    strand:   HashMap<Strand, MethAgg>,
    coverage: HashMap<CountType, usize>,
}

macro_rules! union_hashmap {
    ($left:expr, $right:expr) => {
        $left.keys().chain($right.keys()).unique().map(|&k| {
            (
                k,
                $left.get(&k).cloned().unwrap_or_default()
                    + $right.get(&k).cloned().unwrap_or_default(),
            )
        })
    };
}

impl AddAssign for RegionMethAgg {
    fn add_assign(
        &mut self,
        rhs: Self,
    ) {
        self.context = union_hashmap!(self.context, rhs.context).collect();
        self.strand = union_hashmap!(self.strand, rhs.strand).collect();
        self.coverage = union_hashmap!(self.coverage, rhs.coverage).collect();
    }
}

impl Add for RegionMethAgg {
    type Output = RegionMethAgg;

    fn add(
        mut self,
        rhs: Self,
    ) -> Self::Output {
        self += rhs;
        self
    }
}

impl RegionMethAgg {
    getter_fn!(coverage, HashMap<CountType, usize>);

    getter_fn!(context, HashMap<Context, MethAgg>);

    getter_fn!(strand, HashMap<Strand, MethAgg>);

    pub fn new() -> Self {
        Self::default()
    }

    pub fn full(
        coverage: HashMap<CountType, usize>,
        context: HashMap<Context, MethAgg>,
        strand: HashMap<Strand, MethAgg>,
    ) -> RegionMethAgg {
        Self {
            coverage,
            context,
            strand,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.context.is_empty() && self.strand.is_empty()
    }

    pub fn finalize_context(&self) -> HashMap<Context, DensityType> {
        self.context
            .iter()
            .map(|(&k, v)| (k, v.finalize()))
            .collect()
    }

    pub fn finalize_strand(&self) -> HashMap<Strand, DensityType> {
        self.strand
            .iter()
            .map(|(&k, v)| (k, v.finalize()))
            .collect()
    }

    pub fn add_cytosine(
        &mut self,
        sum: DensityType,
        count: Option<CountType>,
        context: Context,
        strand: Strand,
    ) {
        if let Some(coverage) = count {
            self.coverage.entry(coverage).or_insert(0).add_assign(1);
        }
        let density = sum / count.unwrap_or(1) as DensityType;
        self.strand.entry(strand).or_default().add_density(density);
        self.context
            .entry(context)
            .or_default()
            .add_density(density);
    }
}
