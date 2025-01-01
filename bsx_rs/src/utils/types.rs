
pub trait IPCEncodedEnum {
    fn from_bool(value: Option<bool>) -> Self;
    fn to_bool(&self) -> Option<bool>;
    fn from_str(value: &str) -> Self;
    fn to_string(&self) -> String;
}

pub type BSXResult<T> = Result<T, Box<dyn std::error::Error>>;
#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug)]
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

#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug)]
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_context_from_bool() {
        assert_eq!(Context::from_bool(Some(true)), Context::CG);
        assert_eq!(Context::from_bool(Some(false)), Context::CHG);
        assert_eq!(Context::from_bool(None), Context::CHH);
    }

    #[test]
    fn test_context_to_bool() {
        assert_eq!(Context::CG.to_bool(), Some(true));
        assert_eq!(Context::CHG.to_bool(), Some(false));
        assert_eq!(Context::CHH.to_bool(), None);
    }

    #[test]
    fn test_context_from_str() {
        assert_eq!(Context::from_str("CG"), Context::CG);
        assert_eq!(Context::from_str("CHG"), Context::CHG);
        assert_eq!(Context::from_str("CHH"), Context::CHH);
    }

    #[test]
    #[should_panic]
    fn test_context_from_str_invalid() {
        Context::from_str("invalid");
    }

    #[test]
    fn test_context_to_string() {
        assert_eq!(Context::CG.to_string(), "CG");
        assert_eq!(Context::CHG.to_string(), "CHG");
        assert_eq!(Context::CHH.to_string(), "CHH");
    }

    #[test]
    fn test_strand_from_bool() {
        assert_eq!(Strand::from_bool(Some(true)), Strand::Forward);
        assert_eq!(Strand::from_bool(Some(false)), Strand::Reverse);
        assert_eq!(Strand::from_bool(None), Strand::None);
    }

    #[test]
    fn test_strand_to_bool() {
        assert_eq!(Strand::Forward.to_bool(), Some(true));
        assert_eq!(Strand::Reverse.to_bool(), Some(false));
        assert_eq!(Strand::None.to_bool(), None);
    }

    #[test]
    fn test_strand_from_str() {
        assert_eq!(Strand::from_str("+"), Strand::Forward);
        assert_eq!(Strand::from_str("-"), Strand::Reverse);
        assert_eq!(Strand::from_str("."), Strand::None);
    }

    #[test]
    fn test_strand_from_str_case_insensitive() {
        assert_eq!(Strand::from_str("+"), Strand::Forward);
        assert_eq!(Strand::from_str("-"), Strand::Reverse);
        assert_eq!(Strand::from_str("."), Strand::None);
    }

    #[test]
    fn test_strand_to_string() {
        assert_eq!(Strand::Forward.to_string(), "+");
        assert_eq!(Strand::Reverse.to_string(), "-");
        assert_eq!(Strand::None.to_string(), ".");
    }
}
