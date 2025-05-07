use std::sync::Arc;

use bsxplorer2::io::report::ReportTypeSchema as RustReportTypeSchema;
use pyo3::prelude::*;
use pyo3_polars::PySchema;

#[pyclass(name = "ReportTypeSchema", eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum PyReportTypeSchema {
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
}

impl From<RustReportTypeSchema> for PyReportTypeSchema {
    fn from(rust: RustReportTypeSchema) -> Self {
        PyReportTypeSchema::from_rust(rust)
    }
}

impl From<PyReportTypeSchema> for RustReportTypeSchema {
    fn from(py: PyReportTypeSchema) -> Self { PyReportTypeSchema::to_rust(&py) }
}

#[pymethods]
impl PyReportTypeSchema {
    pub fn col_names(&self) -> Vec<&'static str> {
        self.to_rust().col_names().to_vec()
    }

    pub fn schema(&self) -> PyResult<PySchema> {
        Ok(PySchema(Arc::new(self.to_rust().schema())))
    }

    pub fn chr_col(&self) -> &'static str { self.to_rust().chr_col() }

    pub fn position_col(&self) -> &'static str { self.to_rust().position_col() }

    pub fn context_col(&self) -> Option<&'static str> {
        self.to_rust().context_col()
    }

    pub fn strand_col(&self) -> Option<&'static str> {
        self.to_rust().strand_col()
    }

    pub fn need_align(&self) -> bool { self.to_rust().need_align() }
}

impl PyReportTypeSchema {
    pub fn from_rust(rust: RustReportTypeSchema) -> Self {
        match rust {
            RustReportTypeSchema::Bismark => PyReportTypeSchema::Bismark,
            RustReportTypeSchema::CgMap => PyReportTypeSchema::CgMap,
            RustReportTypeSchema::BedGraph => PyReportTypeSchema::BedGraph,
            RustReportTypeSchema::Coverage => PyReportTypeSchema::Coverage,
        }
    }

    pub fn to_rust(&self) -> RustReportTypeSchema {
        match self {
            PyReportTypeSchema::Bismark => RustReportTypeSchema::Bismark,
            PyReportTypeSchema::CgMap => RustReportTypeSchema::CgMap,
            PyReportTypeSchema::BedGraph => RustReportTypeSchema::BedGraph,
            PyReportTypeSchema::Coverage => RustReportTypeSchema::Coverage,
        }
    }
}
