use log::warn;
use num::Unsigned;
use polars::export::num::PrimInt;
use serde::Serialize;
use std::cmp::Ordering;
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::ops::{Add, Shr, Sub};

#[derive(Debug, Clone, Serialize)]
pub struct GenomicPosition<N: PosNum> {
    chr: String,
    position: N,
}

impl<N: PosNum> Display for GenomicPosition<N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}", self.chr, self.position)
    }
}

impl<N: PosNum> GenomicPosition<N> {
    pub fn new(chr: String, position: N) -> GenomicPosition<N> {
        GenomicPosition { chr, position }
    }

    pub fn chr(&self) -> &str {
        &self.chr
    }

    pub fn position(&self) -> N {
        self.position
    }
}

impl<N: PosNum> Add<N> for GenomicPosition<N> {
    type Output = GenomicPosition<N>;
    fn add(self, rhs: N) -> Self::Output {
        GenomicPosition::new(self.chr, self.position + rhs)
    }
}

impl<N: PosNum> Sub<N> for GenomicPosition<N> {
    type Output = GenomicPosition<N>;
    fn sub(self, rhs: N) -> Self::Output {
        GenomicPosition::new(self.chr, self.position - rhs)
    }
}

impl<N: PosNum> PartialEq for GenomicPosition<N> {
    fn eq(&self, other: &Self) -> bool {
        self.chr == other.chr && self.position == other.position
    }
}

impl<N: PosNum> PartialOrd for GenomicPosition<N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.chr == other.chr {
            Some(self.position.cmp(&other.position))
        } else {
            None
        }
    }
}

impl<N: PosNum> Shr for GenomicPosition<N> {
    type Output = Option<RegionCoordinates<N>>;

    fn shr(self, rhs: Self) -> Self::Output {
        let chr_cmp = self.partial_cmp(&rhs);
        if chr_cmp.is_some() && chr_cmp.unwrap() == Ordering::Greater {
            warn!(
                "Bad range operation for {} and {} regions. Lhs must be less than rhs.",
                self, rhs
            );
            Some(RegionCoordinates::new(
                self.chr,
                rhs.position,
                self.position,
                Strand::None,
            ))
        } else if chr_cmp.is_some() && chr_cmp.unwrap() == Ordering::Less {
            Some(RegionCoordinates::new(
                self.chr,
                self.position,
                rhs.position,
                Strand::None,
            ))
        } else {
            None
        }
    }
}

#[derive(Clone, Hash, PartialEq, Eq, Debug, Serialize)]
pub struct RegionCoordinates<N: PosNum> {
    pub chr: String,
    pub start: N,
    pub end: N,
    pub strand: Strand,
}

impl<N: PosNum> RegionCoordinates<N> {
    pub fn new(chr: String, start: N, end: N, strand: Strand) -> Self {
        assert!(
            start <= end,
            "End position can not be less than start! {}:{}-{}",
            chr,
            start,
            end
        );
        Self {
            chr,
            start,
            end,
            strand,
        }
    }
    pub fn try_from_position(start: GenomicPosition<N>, end: GenomicPosition<N>) -> Option<Self> {
        start >> end
    }

    pub fn into_positions(self) -> (GenomicPosition<N>, GenomicPosition<N>) {
        (
            GenomicPosition::new(self.chr.clone(), self.start),
            GenomicPosition::new(self.chr, self.end),
        )
    }
    pub fn chr(&self) -> &str {
        self.chr.as_str()
    }
    pub fn start(&self) -> N {
        self.start
    }
    pub fn strand(&self) -> Strand {
        self.strand
    }
    pub fn end(&self) -> N {
        self.end
    }
    pub fn start_gpos(&self) -> GenomicPosition<N> {
        GenomicPosition::new(self.chr.clone(), self.start)
    }
    pub fn end_gpos(&self) -> GenomicPosition<N> {
        GenomicPosition::new(self.chr.clone(), self.end)
    }
    pub fn length(&self) -> N {
        self.end - self.start
    }

    pub fn set_start(&mut self, start: N) -> Result<(), Box<dyn Error>> {
        if start < self.end {
            self.start = start;
            Ok(())
        } else {
            Err(Box::from(format!(
                "Could not modify start value. Provided {} start value is greater than end {}.",
                start, self.end
            )))
        }
    }
    pub fn set_end(&mut self, end: N) -> Result<(), Box<dyn Error>> {
        if end < self.start {
            self.end = end;
            Ok(())
        } else {
            Err(Box::from(format!(
                "Could not modify end value. Provided {} end value is less than start {}.",
                end, self.start
            )))
        }
    }
    pub fn set_chromosome(&mut self, chr: String) {
        self.chr = chr;
    }

    /// Expands limits of the regions by specified value. If length of
    /// the resulting region is more than max_length, [Err] is returned.
    /// Also, if start is less than value [Err] is returned.
    /// Otherwise bounds of the region will be modified themselves
    pub fn expand(mut self, value: N, max_length: Option<N>) -> Self {
        let max_length = max_length.unwrap_or(N::max_value());

        if self.end + value > max_length {
            self.end = max_length;
        } else {
            self.end = self.end + value
        }
        if self.start > value {
            self.start = self.start - value;
        }
        self
    }
}

impl<N: PosNum> Display for RegionCoordinates<N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chr, self.start, self.end)
    }
}

#[cfg(feature = "python")]
mod python {
    use crate::data_structs::region::{GenomicPosition, RegionCoordinates};
    use pyo3::prelude::*;

    #[pyclass(name = "GenomicPosition")]
    #[derive(Clone)]
    pub struct PyGenomicPosition {
        pub chr: String,
        pub position: u64,
    }

    impl From<GenomicPosition<u64>> for PyGenomicPosition {
        fn from(pos: GenomicPosition<u64>) -> Self {
            PyGenomicPosition {
                chr: pos.chr,
                position: pos.position,
            }
        }
    }

    impl From<PyGenomicPosition> for GenomicPosition<u64> {
        fn from(pos: PyGenomicPosition) -> Self {
            GenomicPosition {
                chr: pos.chr,
                position: pos.position,
            }
        }
    }

    #[pymethods]
    impl PyGenomicPosition {
        #[new]
        fn new(chr: String, position: u64) -> PyGenomicPosition {
            PyGenomicPosition { chr, position }
        }

        fn __add__(&self, other: &PyGenomicPosition) -> PyGenomicPosition {
            PyGenomicPosition {
                chr: self.chr.clone(),
                position: self.position + other.position,
            }
        }

        fn __sub__(&self, other: &PyGenomicPosition) -> PyGenomicPosition {
            PyGenomicPosition {
                chr: self.chr.clone(),
                position: self.position - other.position,
            }
        }

        fn __rshift__(&self, other: &PyGenomicPosition) -> Option<PyRegionCoordinates> {
            let rhs = GenomicPosition::from(self.clone());
            let lhs = GenomicPosition::from(other.clone());
            (rhs >> lhs).map(PyRegionCoordinates::from)
        }

        #[getter]
        fn get_chr(&self) -> &str {
            self.chr.as_str()
        }
        #[getter]
        fn get_position(&self) -> u64 {
            self.position
        }
    }

    #[pymethods]
    impl PyRegionCoordinates {
        #[new]
        fn new(chr: &str, start: u64, end: u64) -> Self {
            PyRegionCoordinates {
                chr: chr.to_string(),
                start,
                end,
            }
        }

        #[getter]
        fn get_chr(&self) -> &str {
            self.chr.as_str()
        }
        #[getter]
        fn get_start(&self) -> u64 {
            self.start
        }
        #[getter]
        fn get_end(&self) -> u64 {
            self.end
        }
        #[getter]
        fn start_gpos(&self) -> PyGenomicPosition {
            PyGenomicPosition {
                chr: self.chr.clone(),
                position: self.start,
            }
        }
        #[getter]
        fn end_gpos(&self) -> PyGenomicPosition {
            PyGenomicPosition {
                chr: self.chr.clone(),
                position: self.end,
            }
        }
        fn __len__(&self) -> usize {
            self.end as usize - self.start as usize
        }
    }

    #[pyclass(name = "RegionCoordinates")]
    #[derive(Clone)]
    pub struct PyRegionCoordinates {
        pub chr: String,
        pub start: u64,
        pub end: u64,
    }

    impl From<RegionCoordinates<u64>> for PyRegionCoordinates {
        fn from(region: RegionCoordinates<u64>) -> Self {
            PyRegionCoordinates {
                chr: region.chr,
                start: region.start,
                end: region.end,
            }
        }
    }

    impl From<PyRegionCoordinates> for RegionCoordinates<u64> {
        fn from(region: PyRegionCoordinates) -> Self {
            RegionCoordinates {
                chr: region.chr,
                start: region.start,
                end: region.end,
            }
        }
    }
}

use crate::utils::types::PosNum;
use crate::utils::types::Strand;
#[cfg(feature = "python")]
pub(crate) use python::{PyGenomicPosition, PyRegionCoordinates};

#[cfg(test)]
mod tests {
    // TODO
}
