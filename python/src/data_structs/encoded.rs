use pyo3::prelude::*;
use pyo3_polars::{PyDataFrame, PyDataType};
use bsxplorer2::data_structs::batch::BsxBatch as RsBsxBatch;
use bsxplorer2::data_structs::batch::EncodedBsxBatch as RsEncodedBsxBatch;
use bsxplorer2::io::report::ReportTypeSchema as RustReportTypeSchema;
use bsxplorer2::data_structs::context_data::ContextData as RustContextData;
use bsxplorer2::utils::types::{Context as RsContext, Strand as RsStrand};

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
enum ReportTypeSchema {
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
}

impl ReportTypeSchema {
    fn from_rust(rust: RustReportTypeSchema) -> Self {
        match rust {
            RustReportTypeSchema::Bismark => ReportTypeSchema::Bismark,
            RustReportTypeSchema::CgMap => ReportTypeSchema::CgMap,
            RustReportTypeSchema::BedGraph => ReportTypeSchema::BedGraph,
            RustReportTypeSchema::Coverage => ReportTypeSchema::Coverage,
        }
    }

    fn to_rust(&self) -> RustReportTypeSchema {
        match self {
            ReportTypeSchema::Bismark => RustReportTypeSchema::Bismark,
            ReportTypeSchema::CgMap => RustReportTypeSchema::CgMap,
            ReportTypeSchema::BedGraph => RustReportTypeSchema::BedGraph,
            ReportTypeSchema::Coverage => RustReportTypeSchema::Coverage,
        }
    }
}

#[pyclass(eq, eq_int)]
#[derive(PartialEq)]
enum Strand {
    Forward,
    Reverse,
    None,
}

impl Strand {
    fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "forward" => Some(Strand::Forward),
            "reverse" => Some(Strand::Reverse),
            "none" => Some(Strand::None),
            _ => None,
        }
    }
}

impl From<RsStrand> for Strand {
    fn from(s: RsStrand) -> Self {
        match s {
            RsStrand::Forward => Strand::Forward,
            RsStrand::Reverse => Strand::Reverse,
            RsStrand::None => Strand::None,
        }
    }
}

impl From<Strand> for RsStrand {
    fn from(s: Strand) -> Self {
        match s {
            Strand::Forward => RsStrand::Forward,
            Strand::Reverse => RsStrand::Reverse,
            Strand::None => RsStrand::None,
        }
    }
}

#[pyclass(eq, eq_int)]
#[derive(PartialEq)]
enum Context {
    CG,
    CHG,
    CHH,
}

impl Context {
    fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "cg" => Some(Context::CG),
            "chg" => Some(Context::CHG),
            "chh" => Some(Context::CHH),
            _ => None,
        }
    }
}

impl From<RsContext> for Context {
    fn from(s: RsContext) -> Self {
        match s {
            RsContext::CG => Context::CG,
            RsContext::CHG => Context::CHG,
            RsContext::CHH => Context::CHH,
        }
    }
}

impl From<Context> for RsContext {
    fn from(s: Context) -> Self {
        match s {
            Context::CG => RsContext::CG,
            Context::CHG => RsContext::CHG,
            Context::CHH => RsContext::CHH,
        }
    }
}

#[pyclass]
#[derive(Clone)]
struct ContextData {
    data: RustContextData,
}   

#[pymethods]
impl ContextData {
    #[new]
    fn new(sequence: Vec<u8>) -> Self {
        Self {
            data: RustContextData::from_sequence(&sequence),
        }
    }

    fn to_decoded_df(&self) -> PyResult<PyDataFrame> {
        let df = self.data.clone().to_df::<RsBsxBatch>();
        Ok(PyDataFrame(df))
    }

    pub fn to_encoded_df(&self) -> PyResult<PyDataFrame> {
        let df = self.data.clone().to_df::<RsEncodedBsxBatch>();
        Ok(PyDataFrame(df))
    }
}

#[pyclass]
#[derive(Clone)]
struct BsxBatch {
    inner: RsBsxBatch
}


