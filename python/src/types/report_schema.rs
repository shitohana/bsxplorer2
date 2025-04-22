use std::sync::Arc;

use pyo3::prelude::*;
use pyo3_polars::{PySchema};
use bsxplorer2::io::report::ReportTypeSchema as RustReportTypeSchema;

/// Represents different input/output file formats for methylation data.
///
/// Attributes
/// ----------
/// Bismark
///     Represents the Bismark format (cytosine report).
/// CgMap
///     Represents the CGmap format.
/// BedGraph
///     Represents the BedGraph format.
/// Coverage
///     Represents the Bismark coverage format.
#[pyclass(name="ReportTypeSchema", eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum PyReportTypeSchema {
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
}

#[pymethods]
impl PyReportTypeSchema {
    /// Get the list of column names for this report format.
    ///
    /// Returns
    /// -------
    /// list[str]
    ///     A list of column names.
    pub fn col_names(&self) -> Vec<&'static str> {
        self.to_rust().col_names().to_vec()
    }

    /// Get the Polars schema for this report format.
    ///
    /// Returns
    /// -------
    /// Schema
    ///     The Polars schema definition.
    pub fn schema(&self) -> PyResult<PySchema> {
        Ok(PySchema(Arc::new(self.to_rust().schema())))
    }

    /// Get the name of the chromosome column for this format.
    ///
    /// Returns
    /// -------
    /// str
    ///     The chromosome column name.
    pub fn chr_col(&self) -> &'static str {
        self.to_rust().chr_col()
    }

    /// Get the name of the position column for this format.
    ///
    /// Notes
    /// -----
    /// For some formats like BedGraph, this might be 'start'.
    ///
    /// Returns
    /// -------
    /// str
    ///     The position column name.
    pub fn position_col(&self) -> &'static str {
        self.to_rust().position_col()
    }

    /// Get the name of the context column, if available for this format.
    ///
    /// Returns
    /// -------
    /// Optional[str]
    ///     The context column name or None if not applicable.
    pub fn context_col(&self) -> Option<&'static str> {
        self.to_rust().context_col()
    }

    /// Get the name of the strand column, if available for this format.
    ///
    /// Returns
    /// -------
    /// Optional[str]
    ///     The strand column name or None if not applicable.
    pub fn strand_col(&self) -> Option<&'static str> {
        self.to_rust().strand_col()
    }

    /// Check if this report format typically requires alignment with context data.
    ///
    /// Notes
    /// -----
    /// Formats like BedGraph or Coverage often lack explicit strand/context
    /// information and benefit from alignment with genomic context.
    ///
    /// Returns
    /// -------
    /// bool
    ///     True if alignment is generally needed, False otherwise.
    pub fn need_align(&self) -> bool {
        self.to_rust().need_align()
    }
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