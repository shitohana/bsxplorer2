use itertools::Itertools;
use polars::prelude::*;
use std::fmt::Display;
use std::ops::{Mul, Sub};
use crate::utils::types::BSXResult;

#[derive(Copy, Clone, Debug, PartialEq, Hash)]
pub enum FillNullStrategy {
    /// mean value in array
    Mean,
    /// minimal value in array
    Min,
    /// maximum value in array
    Max,
    /// replace with the value zero
    Zero,
    /// replace with the value one
    One,
    /// replace with the maximum value of that data type
    MaxBound,
    /// replace with the minimal value of that data type
    MinBound,
}

#[derive(Clone)]
pub struct NullHandleStrategy {
    skip: bool,
    fill: FillNullStrategy,
}

impl Default for NullHandleStrategy {
    fn default() -> Self {
        Self {
            skip: false,
            fill: FillNullStrategy::Zero,
        }
    }
}

impl NullHandleStrategy {
    pub fn new(skip: bool, fill: FillNullStrategy) -> Self {
        Self { skip, fill }
    }
    pub fn with_fill(mut self, fill: FillNullStrategy) -> Self {
        self.fill = fill;
        self
    }
    pub fn with_skip(mut self, skip: bool) -> Self {
        self.skip = skip;
        self
    }
    fn to_polars(&self) -> Option<polars::chunked_array::ops::FillNullStrategy> {
        if !self.skip {
            None
        } else {
            let strategy = match self.fill {
                FillNullStrategy::Zero => polars::chunked_array::ops::FillNullStrategy::Zero,
                FillNullStrategy::One => polars::chunked_array::ops::FillNullStrategy::One,
                FillNullStrategy::Max => polars::chunked_array::ops::FillNullStrategy::Max,
                FillNullStrategy::Min => polars::chunked_array::ops::FillNullStrategy::Min,
                FillNullStrategy::Mean => polars::chunked_array::ops::FillNullStrategy::Mean,
                FillNullStrategy::MaxBound => {
                    polars::chunked_array::ops::FillNullStrategy::MaxBound
                }
                FillNullStrategy::MinBound => {
                    polars::chunked_array::ops::FillNullStrategy::MinBound
                }
            };
            Some(strategy)
        }
    }
}

pub struct RegionData {
    pub data: DataFrame,
    pub coordinates: RegionCoordinates,
}

impl RegionData {
    pub fn new(data: DataFrame, coordinates: RegionCoordinates) -> Self {
        Self { data, coordinates }
    }

    pub fn drop_null(&mut self, null_strategy: NullHandleStrategy) -> Option<()> {
        for col in ["density", "count_m", "count_total"] {
            if self.data.column(col).unwrap().null_count() > 0 && null_strategy.skip {
                return None;
            } else if !null_strategy.skip {
                self.data = self.data.fill_null(null_strategy.to_polars()?).unwrap();
            }
        }
        Some(())
    }

    pub fn get_coordinates(&self) -> &RegionCoordinates {
        &self.coordinates
    }

    pub(crate) fn from_parts(
        parts: Vec<DataFrame>,
        coordinates: RegionCoordinates,
    ) -> Option<Self> {
        let mut sorted_iter = parts
            .iter()
            .filter(|df| df.height() > 0)
            .sorted_by_key(|df| {
                df.column("position")
                    .unwrap()
                    .u32()
                    .unwrap()
                    .first()
                    .unwrap()
            });
        let mut res = match sorted_iter.next() {
            Some(df) => df.clone(),
            None => return None,
        };
        for df in sorted_iter {
            res.extend(df).unwrap();
        }
        res.align_chunks_par();
        Some(Self {
            data: res,
            coordinates,
        })
    }

    pub fn discretize(&self, resolution: usize) -> DataFrame {
        self.data
            .clone()
            .lazy()
            .with_column(
                col("position")
                    .sub(lit(self.coordinates.start()))
                    .mul(lit(resolution as f64))
                    .floor_div(lit(self.coordinates.length() as f64 + 0.5))
                    .cast(DataType::UInt32)
                    .alias("fragment"),
            )
            .collect()
            .unwrap()
    }
}

impl Display for RegionData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Region {}:{}-{} data with {} rows.",
            self.coordinates.chr(),
            self.coordinates.start(),
            self.coordinates.length(),
            self.data.height()
        )
    }
}

#[derive(Clone, Hash, PartialEq, Eq, Debug)]
pub struct RegionCoordinates {
    pub(crate) chr: String,
    pub(crate) start: u32,
    pub(crate) end: u32,
}

impl RegionCoordinates {
    pub fn new(chr: String, start: u32, end: u32) -> Self {
        Self { chr, start, end }
    }
    pub fn chr(&self) -> &str {
        self.chr.as_str()
    }
    pub fn start(&self) -> u32 {
        self.start
    }
    pub fn end(&self) -> u32 {
        self.end
    }
    pub fn length(&self) -> u32 {
        self.end - self.start
    }
    
    // TODO: Add methods to modify start and end
    
    /// Expands limits of the regions by specified value. If length of 
    /// the resulting region is more than max_length, [Err] is returned.
    /// Also, if start is less than value [Err] is returned.
    /// Otherwise bounds of the region will be modified themselves
    pub fn expand(mut self, value: u32, max_length: u32) -> Self {
        if self.end + value > max_length {
            self.end = max_length;
        } else {
            self.end += value
        }
        if self.start > value {
            self.start -= value;
        }
        self
    }
}

impl Display for RegionCoordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chr, self.start, self.end)
    }
}