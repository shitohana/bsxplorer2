use std::convert::Infallible;
use std::fmt::Display;
use std::hash::Hash;
use std::str::FromStr;

use serde::{
    Deserialize,
    Serialize,
};

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

impl From<Option<bool>> for Context {
    fn from(value: Option<bool>) -> Self {
        match value {
            Some(true) => Context::CG,
            Some(false) => Context::CHG,
            None => Context::CHH,
        }
    }
}

impl From<Context> for Option<bool> {
    fn from(value: Context) -> Self {
        match value {
            Context::CG => Some(true),
            Context::CHG => Some(false),
            Context::CHH => None,
        }
    }
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
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "CG" => Ok(Context::CG),
            "CHG" => Ok(Context::CHG),
            _ => Ok(Context::CHH),
        }
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

impl FromStr for Strand {
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            _ => Ok(Strand::None),
        }
    }
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

impl From<Strand> for Option<bool> {
    fn from(value: Strand) -> Option<bool> {
        match value {
            Strand::Forward => Some(true),
            Strand::Reverse => Some(false),
            Strand::None => None,
        }
    }
}

impl From<Option<bool>> for Strand {
    fn from(value: Option<bool>) -> Strand {
        match value {
            Some(true) => Strand::Forward,
            Some(false) => Strand::Reverse,
            None => Strand::None,
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
