use bsxplorer2::utils::types::{Context as RsContext, Strand as RsStrand};
use pyo3::prelude::*;

/// Represents DNA strand information.
///
/// Attributes
/// ----------
/// Forward
///     Represents the forward strand (+).
/// Reverse
///     Represents the reverse strand (-).
/// None
///     Represents unknown or N/A strand (.).
#[pyclass(eq, eq_int)]
#[derive(PartialEq)]
pub enum Strand {
    /// Forward strand (+).
    Forward,
    /// Reverse strand (-).
    Reverse,
    /// Unknown or N/A strand (.).
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

/// Represents methylation context (sequence context).
///
/// Attributes
/// ----------
/// CG
///     CpG context.
/// CHG
///     CHG context (H = A, C, or T).
/// CHH
///     CHH context (H = A, C, or T).
#[pyclass(eq, eq_int)]
#[derive(PartialEq)]
pub enum Context {
    /// CpG context.
    CG,
    /// CHG context (H = A, C, or T).
    CHG,
    /// CHH context (H = A, C, or T).
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
