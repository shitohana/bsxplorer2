use std::fmt::Display;
use std::hash::Hash;

use serde::{Deserialize, Serialize};

pub trait IPCEncodedEnum {
    /// Converts an optional boolean value to the enum.
    fn from_bool(value: Option<bool>) -> Self;
    /// Converts the enum to an optional boolean value.
    fn to_bool(&self) -> Option<bool>;
    /// Converts a string slice to the enum.
    fn from_str(value: &str) -> Self;
}

pub type BSXResult<T> = Result<T, Box<dyn std::error::Error>>;

#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug, PartialOrd, Ord)]
#[cfg_attr(feature = "console", derive(clap::ValueEnum))]
pub enum Context {
    /// CG context.
    CG,
    /// CHG context.
    CHG,
    /// CHH context.
    CHH,
}

impl Display for Context {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            Context::CG => write!(f, "CG"),
            Context::CHG => write!(f, "CHG"),
            Context::CHH => write!(f, "CHH"),
        }
    }
}

impl std::str::FromStr for Context {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(<Self as IPCEncodedEnum>::from_str(s))
    }
}

impl Serialize for Context {
    fn serialize<S>(
        &self,
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer, {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for Context {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>, {
        let s = String::deserialize(deserializer)?;
        std::str::FromStr::from_str(&s).map_err(serde::de::Error::custom)
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
}

#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug, PartialOrd, Ord)]
#[cfg_attr(feature = "console", derive(clap::ValueEnum))]
pub enum Strand {
    /// Forward strand.
    Forward,
    /// Reverse strand.
    Reverse,
    /// No strand.
    None,
}

impl From<Strand> for bool {
    fn from(value: Strand) -> bool {
        match value {
            Strand::Forward => true,
            Strand::Reverse => false,
            Strand::None => unimplemented!(),
        }
    }
}

impl From<bool> for Strand {
    fn from(value: bool) -> Self {
        match value {
            true => Strand::Forward,
            false => Strand::Reverse,
        }
    }
}

impl From<Strand> for char {
    fn from(value: Strand) -> Self {
        match value {
            Strand::Forward => '+',
            Strand::Reverse => '-',
            Strand::None => '.',
        }
    }
}


impl Display for Strand {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        write!(f, "{}", match self {
            Strand::Forward => String::from("+"),
            Strand::Reverse => String::from("-"),
            Strand::None => String::from("."),
        })
    }
}

impl std::str::FromStr for Strand {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(<Self as IPCEncodedEnum>::from_str(s))
    }
}

impl Serialize for Strand {
    fn serialize<S>(
        &self,
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer, {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for Strand {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>, {
        let s = String::deserialize(deserializer)?;
        std::str::FromStr::from_str(&s).map_err(serde::de::Error::custom)
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
            _ => Strand::None, /* Default to None for any other string
                                * including "." */
        }
    }
}
