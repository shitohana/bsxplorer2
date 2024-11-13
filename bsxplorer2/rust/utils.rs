use pyo3::prelude::*;
use bsx_rs::utils::types::{
    Context as ContextRust,
    Strand as StrandRust
};

#[pyclass]
#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug)]
pub enum Context {
    CG,
    CHG,
    CHH,
    ALL,
}

#[pyclass]
#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug)]
pub enum Strand {
    Forward,
    Reverse,
    None,
}

// Implement conversion from the Python enum Context to the Rust enum ContextRust
impl From<Context> for ContextRust {
    fn from(value: Context) -> Self {
        match value {
            Context::CG => ContextRust::CG,
            Context::CHG => ContextRust::CHG,
            Context::CHH => ContextRust::CHH,
            Context::ALL => ContextRust::ALL,
        }
    }
}

// Implement conversion from the Python enum Strand to the Rust enum StrandRust
impl From<Strand> for StrandRust {
    fn from(value: Strand) -> Self {
        match value {
            Strand::Forward => StrandRust::Forward,
            Strand::Reverse => StrandRust::Reverse,
            Strand::None => StrandRust::None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::mem::variant_count;

    #[test]
    fn test_context_conversion() {
        // Verify the number of variants to ensure exhaustive testing
        assert_eq!(variant_count::<Context>(), 4);
        assert_eq!(variant_count::<ContextRust>(), 4);

        // Test each variant of Context converting correctly to ContextRust
        assert_eq!(ContextRust::from(Context::CG), ContextRust::CG);
        assert_eq!(ContextRust::from(Context::CHG), ContextRust::CHG);
        assert_eq!(ContextRust::from(Context::CHH), ContextRust::CHH);
        assert_eq!(ContextRust::from(Context::ALL), ContextRust::ALL);
    }

    #[test]
    fn test_strand_conversion() {
        // Verify the number of variants to ensure exhaustive testing
        assert_eq!(variant_count::<Strand>(), 3);
        assert_eq!(variant_count::<StrandRust>(), 3);

        // Test each variant of Strand converting correctly to StrandRust
        assert_eq!(StrandRust::from(Strand::Forward), StrandRust::Forward);
        assert_eq!(StrandRust::from(Strand::Reverse), StrandRust::Reverse);
        assert_eq!(StrandRust::from(Strand::None), StrandRust::None);
    }
}
