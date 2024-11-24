use crate::utils::{Context, Strand};
use bsx_rs::io::bsx::bsx_reader::BSXReader as BSXReaderRust;
use bsx_rs::io::bsx::region_reader::ReadFilters as ReadFiltersRust;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::fs::File;

#[pyclass]
pub struct BSXReader {
    inner: BSXReaderRust,
}

#[pymethods]
impl BSXReader {
    #[new]
    pub fn new(file_path: &str, projection: Option<Vec<usize>>) -> PyResult<Self> {
        let handle = File::open(file_path)?;
        let inner = BSXReaderRust::new(handle, projection);
        Ok(Self { inner })
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    pub fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyDataFrame> {
        match slf.inner.next() {
            Some(data) => Some(PyDataFrame(data.into())),
            None => None,
        }
    }
}

#[pyclass]
#[derive(Clone)]
pub struct ReadFilters {
    #[pyo3(get, set)]
    context: Option<Context>, // Wrapper enum for ContextRust
    #[pyo3(get, set)]
    strand: Option<Strand>, // Wrapper enum for StrandRust
}

#[pymethods]
impl ReadFilters {
    #[new]
    pub fn new(context: Option<Context>, strand: Option<Strand>) -> PyResult<Self> {
        Ok(Self { context, strand })
    }

    /// Constructor for no filtering
    #[staticmethod]
    pub fn empty() -> PyResult<Self> {
        Self::new(None, None)
    }
}

impl ReadFilters {
    /// Convert Python ReadFilters to native Rust ReadFilters
    fn into_rust(self) -> ReadFiltersRust {
        ReadFiltersRust::new(
            self.context.map(|c| c.into()), // Convert Option<Context> to Option<ContextRust>
            self.strand.map(|s| s.into()),  // Convert Option<Strand> to Option<StrandRust>
        )
    }
}
