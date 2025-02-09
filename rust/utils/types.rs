use num::{PrimInt, Unsigned};
#[cfg(feature = "python")]
use pyo3::pyclass;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::fmt::Display;
use std::hash::Hash;

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

impl PartialOrd<Self> for Context {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Context {
    fn cmp(&self, other: &Self) -> Ordering {
        if self == other {
            Ordering::Equal
        } else {
            match (self, other) {
                (Context::CG, _) => Ordering::Greater,
                (_, Context::CG) => Ordering::Less,
                (Context::CHG, _) => Ordering::Greater,
                (_, Context::CHG) => Ordering::Less,
                (Context::CHH, _) => Ordering::Greater,
                (_, Context::CHH) => Ordering::Less,
            }
        }
    }
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
            other => unimplemented!("Context {} is not yet supported", other),
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

impl PartialOrd<Self> for Strand {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Strand {
    fn cmp(&self, other: &Self) -> Ordering {
        if self == other {
            Ordering::Equal
        } else {
            match (self, other) {
                (Strand::Forward, _) => Ordering::Greater,
                (_, Strand::Forward) => Ordering::Less,
                (Strand::Reverse, _) => Ordering::Greater,
                (_, Strand::Reverse) => Ordering::Less,
                (Strand::None, _) => Ordering::Greater,
                (_, Strand::None) => Ordering::Less,
            }
        }
    }
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

pub trait PosNum: PrimInt + Unsigned + Clone + Serialize + Display + Sync {}

impl<T> PosNum for T where T: PrimInt + Unsigned + Clone + Serialize + Display + Sync {}

pub trait RefId: Eq + Hash + From<String> + Clone + Default {}

impl<T> RefId for T where T: Eq + Hash + From<String> + Clone + Default {}

pub trait Data: Sized + Clone {}

impl<T> Data for T where T: Sized + Clone {}
