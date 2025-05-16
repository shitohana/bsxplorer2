use std::collections::HashMap;
use std::str::FromStr;

use bsxplorer2::data_structs::enums::{Context, Strand};
use bsxplorer2::data_structs::methstats::MethylationStats;
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
        mean_methylation: f64,
        variance_methylation: f64,
        coverage_distribution: HashMap<u16, u32>,
        context_methylation: HashMap<String, (f64, u32)>,
        strand_methylation: HashMap<String, (f64, u32)>,
    ) -> Self {
        use hashbrown::HashMap as HashbrownMap;
        let coverage_distribution: HashbrownMap<u16, u32> =
            HashbrownMap::from_iter(coverage_distribution.into_iter());
        let context_distribution: HashbrownMap<Context, (f64, u32)> =
            HashbrownMap::from_iter(
                context_methylation
                    .into_iter()
                    .map(|(key, value)| (Context::from_str(&key).unwrap(), value)),
            );
        let strand_distribution: HashbrownMap<Strand, (f64, u32)> =
            HashbrownMap::from_iter(
                strand_methylation
                    .into_iter()
                    .map(|(key, value)| (Strand::from_str(&key).unwrap(), value)),
            );
        PyMethylationStats {
            inner: MethylationStats::from_data(
                mean_methylation,
                variance_methylation,
                coverage_distribution,
                context_distribution,
                strand_distribution,
            ),
        }
    }

    pub fn merge(
        &mut self,
        other: &PyMethylationStats,
    ) {
        self.inner.merge(&other.inner);
    }

    pub fn finalize_methylation(&mut self) {
        self.inner.finalize_methylation();
    }

    pub fn total_coverage(&self) -> u32 {
        self.inner.total_coverage()
    }

    pub fn mean_methylation(&self) -> f64 {
        self.inner.mean_methylation()
    }

    #[staticmethod]
    pub fn merge_multiple(
        py: Python,
        stats_list: &Bound<PyAny>,
    ) -> PyResult<Self> {
        let list: &Bound<PyList> = stats_list.downcast::<PyList>()?;
        let stats: Vec<PyMethylationStats> = list.extract()?;
        let inner_stats = stats.into_iter().map(|s| s.inner).collect_vec();
        Ok(PyMethylationStats {
            inner: MethylationStats::merge_multiple(&inner_stats),
        })
    }

    pub fn coverage_distribution(&self) -> HashMap<u16, u32> {
        self.inner
            .coverage_distribution()
            .into_iter()
            .map(|(key, value)| (*key, *value))
            .collect()
    }

    pub fn methylation_var(&self) -> f64 {
        self.inner.methylation_var()
    }

    pub fn context_methylation(&self) -> HashMap<String, (f64, u32)> {
        self.inner
            .context_methylation()
            .into_iter()
            .map(|(key, value)| (key.to_string(), *value))
            .collect()
    }

    pub fn strand_methylation(&self) -> HashMap<String, (f64, u32)> {
        self.inner
            .strand_methylation()
            .into_iter()
            .map(|(key, value)| (key.to_string(), *value))
            .collect()
    }
}
