use bsxplorer2::data_structs::enums::{Context as RsContext, Strand as RsStrand};
use pyo3::prelude::*;

#[pyclass(name = "Strand", eq, eq_int, hash, frozen)]
#[derive(PartialEq, Debug, Clone, Copy, Hash, Eq)]
pub enum PyStrand {
    Forward,
    Reverse,
    Null,
}

impl PyStrand {
    fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "forward" => Some(PyStrand::Forward),
            "reverse" => Some(PyStrand::Reverse),
            "null" => Some(PyStrand::Null),
            _ => None,
        }
    }
}

#[pymethods]
impl PyStrand {
    #[getter]
    fn name(&self) -> &str {
        match self {
            PyStrand::Forward => "Forward",
            PyStrand::Reverse => "Reverse",
            PyStrand::Null => "Null",
        }
    }
}

impl From<RsStrand> for PyStrand {
    fn from(s: RsStrand) -> Self {
        match s {
            RsStrand::Forward => PyStrand::Forward,
            RsStrand::Reverse => PyStrand::Reverse,
            RsStrand::None => PyStrand::Null,
        }
    }
}

impl From<PyStrand> for RsStrand {
    fn from(s: PyStrand) -> Self {
        match s {
            PyStrand::Forward => RsStrand::Forward,
            PyStrand::Reverse => RsStrand::Reverse,
            PyStrand::Null => RsStrand::None,
        }
    }
}

#[pyclass(name = "Context", eq, eq_int, hash, frozen)]
#[derive(PartialEq, Debug, Clone, Copy, Hash, Eq)]
pub enum PyContext {
    CG,
    CHG,
    CHH,
}

impl PyContext {
    fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "cg" => Some(PyContext::CG),
            "chg" => Some(PyContext::CHG),
            "chh" => Some(PyContext::CHH),
            _ => None,
        }
    }
}

#[pymethods]
impl PyContext {
    #[getter]
    fn name(&self) -> &str {
        match self {
            PyContext::CG => "CG",
            PyContext::CHG => "CHG",
            PyContext::CHH => "CHH",
        }
    }
}

impl From<RsContext> for PyContext {
    fn from(s: RsContext) -> Self {
        match s {
            RsContext::CG => PyContext::CG,
            RsContext::CHG => PyContext::CHG,
            RsContext::CHH => PyContext::CHH,
        }
    }
}

impl From<PyContext> for RsContext {
    fn from(s: PyContext) -> Self {
        match s {
            PyContext::CG => RsContext::CG,
            PyContext::CHG => RsContext::CHG,
            PyContext::CHH => RsContext::CHH,
        }
    }
}
