use std::fmt::Debug;

use bsxplorer2::data_structs::batch::{BsxBatch, LazyBsxBatch};
use pyo3::prelude::*;
use pyo3_polars::error::PyPolarsErr;

use super::batch::PyBsxBatch;
use super::utils::{PyContext, PyStrand};

#[pyclass(name = "LazyBsxBatch")]
#[derive(Clone)]
pub struct PyLazyBsxBatch {
    inner: LazyBsxBatch,
}

impl Debug for PyLazyBsxBatch {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        write!(f, "PyLazyBsxBatch")
    }
}

impl From<PyLazyBsxBatch> for LazyBsxBatch {
    fn from(value: PyLazyBsxBatch) -> Self {
        value.inner
    }
}

impl From<LazyBsxBatch> for PyLazyBsxBatch {
    fn from(value: LazyBsxBatch) -> Self {
        Self { inner: value }
    }
}

#[pymethods]
impl PyLazyBsxBatch {
    #[new]
    pub fn new(batch: PyBsxBatch) -> Self {
        Self {
            inner: BsxBatch::from(batch).lazy(),
        }
    }

    // Based on LazyBsxBatch::collect -> PolarsResult<BsxBatch>
    pub fn collect(&self) -> PyResult<PyBsxBatch> {
        self.inner
            .clone()
            .collect()
            .map(|val: BsxBatch| -> PyBsxBatch { val.into() })
            .map_err(|e| PyPolarsErr::Polars(e).into())
    }

    // Based on LazyBsxBatch::filter_pos_lt
    pub fn filter_pos_lt(
        &self,
        pos: u32,
    ) -> Self {
        Self {
            inner: self.inner.clone().filter_pos_lt(pos),
        }
    }

    // Based on LazyBsxBatch::filter_pos_gt
    pub fn filter_pos_gt(
        &self,
        pos: u32,
    ) -> Self {
        Self {
            inner: self.inner.clone().filter_pos_gt(pos),
        }
    }

    // Based on LazyBsxBatch::filter_coverage_lt
    pub fn filter_coverage_lt(
        &self,
        coverage: u32,
    ) -> Self {
        Self {
            inner: self.inner.clone().filter_coverage_lt(coverage),
        }
    }

    // Based on LazyBsxBatch::filter_strand
    pub fn filter_strand(
        &self,
        strand: PyStrand,
    ) -> Self {
        Self {
            inner: self.inner.clone().filter_strand(strand.into()),
        }
    }

    // Based on LazyBsxBatch::filter_context
    pub fn filter_context(
        &self,
        context: PyContext,
    ) -> Self {
        Self {
            inner: self.inner.clone().filter_context(context.into()),
        }
    }
}
