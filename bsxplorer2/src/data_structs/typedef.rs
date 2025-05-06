use smallstr::SmallString;

pub(crate) const SMALLSTR_SIZE: usize = 20;
pub(crate) type BsxSmallStr = SmallString<[u8; SMALLSTR_SIZE]>;
