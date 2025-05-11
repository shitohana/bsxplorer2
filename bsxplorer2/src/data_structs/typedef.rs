use std::hash::Hash;

use num::{PrimInt, Unsigned};
use smallstr::SmallString;

pub const SMALLSTR_SIZE: usize = 20;
pub type BsxSmallStr = SmallString<[u8; SMALLSTR_SIZE]>;

pub trait SeqNameStr:
    for<'a> From<&'a str> + AsRef<str> + Clone + Eq + PartialEq + Hash {
}

pub trait SeqPosNum: Unsigned + PrimInt {}

impl<T> SeqPosNum for T where T: Unsigned + PrimInt {}

impl<T> SeqNameStr for T where
    T: AsRef<str> + Clone + Eq + PartialEq + Hash + for<'a> From<&'a str>
{
}
