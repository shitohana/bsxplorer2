#[cfg(feature = "python")]
use pyo3::pyclass;
use serde::{Deserialize, Serialize};

pub trait IPCEncodedEnum {
    fn from_bool(value: Option<bool>) -> Self;
    fn to_bool(&self) -> Option<bool>;
    fn from_str(value: &str) -> Self;
    fn to_string(&self) -> String;
}

pub type BSXResult<T> = Result<T, Box<dyn std::error::Error>>;
#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "python", pyclass)]
pub enum Context {
    CG,
    CHG,
    CHH,
}

impl IPCEncodedEnum for Context {
    fn from_bool(value: Option<bool>) -> Context {
        match value {
            Some(true) => Context::CG,
            Some(false) => Context::CHG,
            None => Context::CHH,
        }
    }
    fn to_bool(&self) -> Option<bool> {
        match self {
            Context::CG => Some(true),
            Context::CHG => Some(false),
            Context::CHH => None,
        }
    }

    fn from_str(value: &str) -> Self {
        match value.to_uppercase().as_str() {
            "CG" => Context::CG,
            "CHG" => Context::CHG,
            "CHH" => Context::CHH,
            _ => panic!(),
        }
    }

    fn to_string(&self) -> String {
        match self {
            Context::CG => String::from("CG"),
            Context::CHG => String::from("CHG"),
            Context::CHH => String::from("CHH"),
        }
    }
}

#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "python", pyclass)]
pub enum Strand {
    Forward,
    Reverse,
    None,
}

impl IPCEncodedEnum for Strand {
    fn from_bool(value: Option<bool>) -> Strand {
        match value {
            Some(true) => Strand::Forward,
            Some(false) => Strand::Reverse,
            None => Strand::None,
        }
    }
    fn to_bool(&self) -> Option<bool> {
        match self {
            Strand::Forward => Some(true),
            Strand::Reverse => Some(false),
            Strand::None => None,
        }
    }
    fn from_str(value: &str) -> Self {
        match value.to_lowercase().as_str() {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::None,
        }
    }
    fn to_string(&self) -> String {
        match self {
            Strand::Forward => String::from("+"),
            Strand::Reverse => String::from("-"),
            Strand::None => String::from("."),
        }
    }
}
