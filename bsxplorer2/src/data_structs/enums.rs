use std::cmp::Ordering;
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

#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug)]
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

impl PartialOrd<Self> for Context {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Context {
    fn cmp(
        &self,
        other: &Self,
    ) -> Ordering {
        // Define an explicit order
        let self_val = match self {
            Context::CG => 2,
            Context::CHG => 1,
            Context::CHH => 0,
        };
        let other_val = match other {
            Context::CG => 2,
            Context::CHG => 1,
            Context::CHH => 0,
        };
        self_val.cmp(&other_val)
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

#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug)]
#[cfg_attr(feature = "console", derive(clap::ValueEnum))]
pub enum Strand {
    /// Forward strand.
    Forward,
    /// Reverse strand.
    Reverse,
    /// No strand.
    None,
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

impl<T: AsRef<str>> From<T> for Strand {
    fn from(value: T) -> Self {
        match value.as_ref() {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::None,
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

impl PartialOrd<Self> for Strand {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Strand {
    fn cmp(
        &self,
        other: &Self,
    ) -> Ordering {
        // Define an explicit order
        let self_val = match self {
            Strand::Forward => 2,
            Strand::Reverse => 1,
            Strand::None => 0,
        };
        let other_val = match other {
            Strand::Forward => 2,
            Strand::Reverse => 1,
            Strand::None => 0,
        };
        self_val.cmp(&other_val)
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

#[cfg(test)]
mod tests {
    use std::panic;

    use super::*;

    // --- Context Tests ---

    #[test]
    fn test_context_from_str() {
        assert_eq!(Context::from_str("CG"), Context::CG);
        assert_eq!(Context::from_str("cg"), Context::CG);
        assert_eq!(Context::from_str("CHG"), Context::CHG);
        assert_eq!(Context::from_str("chg"), Context::CHG);
        assert_eq!(Context::from_str("CHH"), Context::CHH);
        assert_eq!(Context::from_str("chh"), Context::CHH);
    }

    #[test]
    fn test_context_from_str_unimplemented() {
        let result = panic::catch_unwind(|| Context::from_str("XYZ"));
        assert!(result.is_err()); // Expecting a panic
    }

    #[test]
    fn test_context_ord() {
        assert_eq!(Context::CG.cmp(&Context::CG), Ordering::Equal);
        assert_eq!(Context::CHG.cmp(&Context::CHG), Ordering::Equal);
        assert_eq!(Context::CHH.cmp(&Context::CHH), Ordering::Equal);

        assert_eq!(Context::CG.cmp(&Context::CHG), Ordering::Greater);
        assert_eq!(Context::CG.cmp(&Context::CHH), Ordering::Greater);
        assert_eq!(Context::CHG.cmp(&Context::CHH), Ordering::Greater);

        assert_eq!(Context::CHG.cmp(&Context::CG), Ordering::Less);
        assert_eq!(Context::CHH.cmp(&Context::CG), Ordering::Less);
        assert_eq!(Context::CHH.cmp(&Context::CHG), Ordering::Less);
    }

    #[test]
    fn test_context_ipc_encoded() {
        assert_eq!(Context::from_bool(Some(true)), Context::CG);
        assert_eq!(Context::from_bool(Some(false)), Context::CHG);
        assert_eq!(Context::from_bool(None), Context::CHH);

        assert_eq!(Context::CG.to_bool(), Some(true));
        assert_eq!(Context::CHG.to_bool(), Some(false));
        assert_eq!(Context::CHH.to_bool(), None);
    }

    // --- Strand Tests ---

    #[test]
    fn test_strand_from_str() {
        assert_eq!(Strand::from_str("+"), Strand::Forward);
        assert_eq!(Strand::from_str("-"), Strand::Reverse);
        assert_eq!(Strand::from_str("."), Strand::None);
        assert_eq!(Strand::from_str("forward"), Strand::None); // Defaults to None
        assert_eq!(Strand::from_str(""), Strand::None);
        assert_eq!(Strand::from_str("AnythingElse"), Strand::None);
    }

    #[test]
    fn test_strand_ord() {
        assert_eq!(Strand::Forward.cmp(&Strand::Forward), Ordering::Equal);
        assert_eq!(Strand::Reverse.cmp(&Strand::Reverse), Ordering::Equal);
        assert_eq!(Strand::None.cmp(&Strand::None), Ordering::Equal);

        assert_eq!(Strand::Forward.cmp(&Strand::Reverse), Ordering::Greater);
        assert_eq!(Strand::Forward.cmp(&Strand::None), Ordering::Greater);
        assert_eq!(Strand::Reverse.cmp(&Strand::None), Ordering::Greater);

        assert_eq!(Strand::Reverse.cmp(&Strand::Forward), Ordering::Less);
        assert_eq!(Strand::None.cmp(&Strand::Forward), Ordering::Less);
        assert_eq!(Strand::None.cmp(&Strand::Reverse), Ordering::Less);
    }

    #[test]
    fn test_strand_ipc_encoded() {
        assert_eq!(Strand::from_bool(Some(true)), Strand::Forward);
        assert_eq!(Strand::from_bool(Some(false)), Strand::Reverse);
        assert_eq!(Strand::from_bool(None), Strand::None);

        assert_eq!(Strand::Forward.to_bool(), Some(true));
        assert_eq!(Strand::Reverse.to_bool(), Some(false));
        assert_eq!(Strand::None.to_bool(), None);
    }
}
