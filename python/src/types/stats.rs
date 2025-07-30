use std::collections::HashMap;
use std::str::FromStr;

use bsxplorer2::data_structs::typedef::DensityType;
use bsxplorer2::data_structs::{
    Context,
    MethylationStats,
    Strand,
};
use itertools::Itertools;
use pyo3::prelude::*;
use pyo3::types::PyList;

#[pyclass(name = "MethylationStats")]
#[derive(Debug, Clone)]
pub struct PyMethylationStats {
    inner: MethylationStats,
}

impl From<MethylationStats> for PyMethylationStats {
    fn from(value: MethylationStats) -> Self {
        Self { inner: value }
    }
}

impl From<PyMethylationStats> for MethylationStats {
    fn from(value: PyMethylationStats) -> Self {
        value.inner
    }
}

#[pymethods]
impl PyMethylationStats {
    #[new]
    pub fn new() -> Self {
        PyMethylationStats {
            inner: MethylationStats::new(),
        }
    }

    #[staticmethod]
    pub fn from_data(
        mean_methylation: DensityType,
        variance_methylation: DensityType,
        coverage_distribution: HashMap<u16, u32>,
        context_methylation: HashMap<String, (DensityType, u32)>,
        strand_methylation: HashMap<String, (DensityType, u32)>,
    ) -> Self {
        use hashbrown::HashMap as HashbrownMap;
        let coverage_distribution: HashbrownMap<u16, u32> =
            HashbrownMap::from_iter(coverage_distribution.into_iter());
        let context_distribution: HashbrownMap<Context, (DensityType, u32)> =
            HashbrownMap::from_iter(
                context_methylation
                    .into_iter()
                    .map(|(key, value)| (Context::from_str(&key).unwrap(), value)),
            );
        let strand_distribution: HashbrownMap<Strand, (DensityType, u32)> =
            HashbrownMap::from_iter(
                strand_methylation
                    .into_iter()
                    .map(|(key, value)| (Strand::from_str(&key).unwrap(), value)),
            );
        PyMethylationStats {
            inner: MethylationStats::from_data(
                mean_methylation,
                coverage_distribution,
                context_distribution,
                strand_distribution,
            ),
        }
    }

    pub fn finalize_methylation(&mut self) {
        self.inner.finalize_methylation();
    }

    pub fn total_coverage(&self) -> u32 {
        self.inner.total_coverage()
    }

    pub fn mean_methylation(&self) -> DensityType {
        self.inner.mean_methylation()
    }

    pub fn coverage_distribution(&self) -> HashMap<u16, u32> {
        self.inner
            .coverage_distribution()
            .into_iter()
            .map(|(key, value)| (*key, *value))
            .collect()
    }

    pub fn context_methylation(&self) -> HashMap<String, (DensityType, u32)> {
        self.inner
            .context_methylation()
            .into_iter()
            .map(|(key, value)| (key.to_string(), *value))
            .collect()
    }

    pub fn strand_methylation(&self) -> HashMap<String, (DensityType, u32)> {
        self.inner
            .strand_methylation()
            .into_iter()
            .map(|(key, value)| (key.to_string(), *value))
            .collect()
    }
}
