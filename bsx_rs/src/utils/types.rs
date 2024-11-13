#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug)]
pub enum Context {
    CG,
    CHG,
    CHH,
    ALL,
}

impl Context {
    pub fn from_bool(value: Option<bool>) -> Self {
        match value {
            Some(true) => Context::CG,
            Some(false) => Context::CHG,
            None => Context::CHH,
        }
    }

    pub fn from_str(value: &str) -> Self {
        match value.to_lowercase().as_str() {
            "all" => Context::ALL,
            "cg" => Context::CG,
            "chg" => Context::CHG,
            "chh" => Context::CHH,
            _ => Context::ALL,
        }
    }
}

#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug)]
pub enum Strand {
    Forward,
    Reverse,
    None,
}

impl Strand {
    pub fn from_bool(value: Option<bool>) -> Self {
        match value {
            Some(true) => Strand::Forward,
            Some(false) => Strand::Reverse,
            None => Strand::None,
        }
    }

    pub fn from_str(value: &str) -> Self {
        match value {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::None,
        }
    }
}
