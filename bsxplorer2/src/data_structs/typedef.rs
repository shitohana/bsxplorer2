use std::hash::Hash;

use arcstr::ArcStr;
use num::{
    PrimInt,
    Unsigned,
};

// TODO:    Test, whether using ArcStr improves performance
//          if no: go back to SmallStr
pub type BsxSmallStr = ArcStr;
pub type PosType = u32;
pub type CountType = u16;
pub type DensityType = f32;

pub trait SeqNameStr:
    for<'a> From<&'a str> + AsRef<str> + Clone + Eq + PartialEq + Hash {
}

pub trait SeqPosNum: Unsigned + PrimInt {}

impl<T> SeqPosNum for T where T: Unsigned + PrimInt {}

impl<T> SeqNameStr for T where
    T: AsRef<str> + Clone + Eq + PartialEq + Hash + for<'a> From<&'a str>
{
}
