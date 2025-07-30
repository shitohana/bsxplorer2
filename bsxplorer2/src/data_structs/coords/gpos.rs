use std::fmt::Display;
use std::ops::{
    Add,
    Sub,
};

use bio::bio_types::annot::loc::Loc;
use bio::bio_types::strand::NoStrand;
use serde::{
    Deserialize,
    Serialize,
};

use crate::data_structs::typedef::{
    BsxSmallStr,
    PosType,
};

/// Represents a genomic position with a sequence name and a position.
#[derive(Debug, Clone, Serialize, Deserialize, Default, Eq, PartialEq, Hash)]
pub struct GenomicPosition {
    seqname:  BsxSmallStr,
    position: PosType,
}

#[allow(clippy::non_canonical_partial_ord_impl)]
impl PartialOrd for GenomicPosition {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<std::cmp::Ordering> {
        if self.seqname == other.seqname {
            self.position.partial_cmp(&other.position)
        }
        else {
            None
        }
    }
}

impl Ord for GenomicPosition {
    fn cmp(
        &self,
        other: &Self,
    ) -> std::cmp::Ordering {
        self.partial_cmp(other)
            .expect("GenomicPosition should have equal sequence names to be comparable")
    }
}

impl GenomicPosition {
    /// Creates a new `GenomicPosition`.
    pub fn new(
        seqname: BsxSmallStr,
        position: PosType,
    ) -> Self {
        Self { seqname, position }
    }

    /// Returns the sequence name.
    pub fn seqname(&self) -> &BsxSmallStr {
        &self.seqname
    }

    /// Returns the position.
    pub fn position(&self) -> PosType {
        self.position
    }

    pub fn is_zero(&self) -> bool {
        self.position == PosType::default() && self.seqname == BsxSmallStr::default()
    }

    pub fn shift(
        mut self,
        shift: isize,
    ) -> GenomicPosition {
        if shift > 0 {
            self.position += shift.unsigned_abs() as PosType
        }
        else {
            self.position -= shift.unsigned_abs() as PosType
        };
        self
    }
}

impl From<bio::bio_types::annot::pos::SeqPosUnstranded> for GenomicPosition {
    /// Converts from `bio_types::annot::pos::SeqPosUnstranded`.
    fn from(value: bio::bio_types::annot::pos::SeqPosUnstranded) -> Self {
        Self {
            seqname:  BsxSmallStr::from(value.refid().clone()),
            position: value.pos() as PosType, // Convert isize to PosType (u32)
        }
    }
}

impl From<GenomicPosition> for bio::bio_types::annot::pos::SeqPosUnstranded {
    /// Converts into `bio_types::annot::pos::SeqPosUnstranded`.\
    fn from(value: GenomicPosition) -> Self {
        bio::bio_types::annot::pos::SeqPosUnstranded::new(
            value.seqname.to_string(),
            value.position as isize, // Convert PosType (u32) to isize
            NoStrand::Unknown,
        )
    }
}

impl Sub for GenomicPosition {
    type Output = Option<Self>;

    /// Subtracts two `GenomicPosition`s.
    ///
    /// Returns `None` if the sequence names are different or if the right-hand
    /// side position is greater than the left-hand side position.
    fn sub(
        self,
        rhs: Self,
    ) -> Self::Output {
        if self.seqname != rhs.seqname {
            return None;
        }
        if self.position < rhs.position {
            return None;
        }
        let new_position = self.position - rhs.position;
        Some(Self::new(self.seqname, new_position))
    }
}

impl Add for GenomicPosition {
    type Output = Option<Self>;

    /// Adds two `GenomicPosition`s.
    ///
    /// Returns `None` if the sequence names are different.\
    fn add(
        self,
        rhs: Self,
    ) -> Self::Output {
        if self.seqname != rhs.seqname {
            return None;
        }
        let new_position = self.position + rhs.position;
        Some(Self::new(self.seqname, new_position))
    }
}

impl Display for GenomicPosition {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        write!(f, "{}:{}", self.seqname, self.position)
    }
}
