use std::fs::File;
use polars::io::mmap::MmapBytesReader;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use bsx_rs::io::report::read::{
    ReportReader as ReportReaderRust,
};
use bsx_rs::io::report::types::{
    ReportType as ReportTypeRust,
};

#[pyclass]
#[derive(Clone)]
pub enum ReportType {
    BISMARK,
    CGMAP,
    BEDGRAPH,
    COVERAGE
}

impl ReportType {
    pub(crate) fn into_rust(self) -> ReportTypeRust {
        match self {
            Self::BEDGRAPH => ReportTypeRust::BEDGRAPH,
            Self::BISMARK => ReportTypeRust::BISMARK,
            Self::COVERAGE => ReportTypeRust::COVERAGE,
            Self::CGMAP => ReportTypeRust::CGMAP,
        }
    }

    pub(crate) fn from_rust(rust_type: &ReportTypeRust) -> Self {
        match rust_type {
            ReportTypeRust::BISMARK => Self::BISMARK,
            ReportTypeRust::BEDGRAPH => Self::BEDGRAPH,
            ReportTypeRust::COVERAGE => Self::COVERAGE,
            ReportTypeRust::CGMAP => Self::CGMAP,
        }
    }
}

#[pyclass]
pub struct ReportReader {
    inner: ReportReaderRust
}

#[pymethods]
impl ReportReader {
    #[new]
    pub fn new(
        report_type: ReportType,
        report_path: String,
        chunk_size: Option<usize>,
        low_memory: Option<bool>,
        n_threads: Option<usize>,
        batch_per_read: Option<usize>,
        batch_size: Option<usize>,
        fa_path: Option<String>,
        fai_path: Option<String>,
    ) -> PyResult<Self> {
        let report_type = report_type.into_rust();
        let handle: Box<dyn MmapBytesReader> = Box::new(File::open(&report_path)?);
        let reader = report_type.get_reader(
            handle,
            chunk_size,
            low_memory,
            n_threads,
            Some(true)
        );

        let inner = ReportReaderRust::new(
            report_type,
            reader,
            batch_per_read, batch_size,
            fa_path, fai_path
        );
        Ok(Self {inner})
    }

    #[getter]
    pub fn get_report_type(slf: PyRef<'_, Self>) -> ReportType {
        let report_type = slf.inner.get_report_type();
        ReportType::from_rust(report_type)
    }
    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    pub fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyDataFrame> {
        match slf.inner.next() {
            Some(data) => Some(PyDataFrame(data)),
            None => None
        }
    }
}