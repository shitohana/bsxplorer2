use std::fmt::Display;
use std::ops::{Add, Sub};
use std::str::FromStr;

use bio::bio_types::annot::loc::Loc;
use bio::bio_types::strand::NoStrand;
use serde::{Deserialize, Serialize};

use crate::data_structs::typedef::{SeqNameStr, SeqPosNum};

/// Represents a genomic position with a sequence name and a position.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomicPosition<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum, {
    seqname:  R,
    position: P,
}

impl<R, P> GenomicPosition<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    /// Creates a new `GenomicPosition`.
    pub fn new(
        seqname: R,
        position: P,
    ) -> Self {
        Self { seqname, position }
    }

    /// Returns the sequence name.
    pub fn seqname(&self) -> R {
        self.seqname.clone()
    }

    /// Returns the position.
    pub fn position(&self) -> P {
        self.position
    }

    pub fn is_zero(&self) -> bool {
        self.position == P::zero() && self.seqname.as_ref() == ""
    }
}

impl<P, S> From<bio::bio_types::annot::pos::SeqPosUnstranded>
    for GenomicPosition<S, P>
where
    S: SeqNameStr + FromStr,
    P: SeqPosNum,
    <S as FromStr>::Err: std::fmt::Debug,
{
    /// Converts from `bio_types::annot::pos::SeqPosUnstranded`.
    fn from(value: bio::bio_types::annot::pos::SeqPosUnstranded) -> Self {
        Self {
            seqname:  S::from_str(value.refid()).unwrap(),
            position: P::from(value.pos())
                .expect("Failed to convert position to P"),
        }
    }
}

impl<R, P> From<GenomicPosition<R, P>>
    for bio::bio_types::annot::pos::SeqPosUnstranded
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    /// Converts into `bio_types::annot::pos::SeqPosUnstranded`.
    fn from(value: GenomicPosition<R, P>) -> Self {
        bio::bio_types::annot::pos::SeqPosUnstranded::new(
            value.seqname.as_ref().to_string(),
            value
                .position
                .to_isize()
                .expect("Failed to convert position to isize"),
            NoStrand::Unknown,
        )
    }
}

impl<R, P> Sub for GenomicPosition<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    type Output = Option<Self>;

    /// Subtracts two `GenomicPosition`s.
    ///
    /// Returns `None` if the sequence names are different or if the right-hand
    /// side position is greater than the left-hand side position.
    fn sub(
        self,
        rhs: Self,
    ) -> Self::Output {
        if self.seqname.as_ref() != rhs.seqname.as_ref() {
            return None;
        }
        if self.position < rhs.position {
            return None;
        }
        let new_position = self.position - rhs.position;
        Some(Self::new(self.seqname, new_position))
    }
}

impl<R, P> Add for GenomicPosition<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    type Output = Option<Self>;

    /// Adds two `GenomicPosition`s.
    ///
    /// Returns `None` if the sequence names are different.
    fn add(
        self,
        rhs: Self,
    ) -> Self::Output {
        if self.seqname.as_ref() != rhs.seqname.as_ref() {
            return None;
        }
        let new_position = self.position + rhs.position;
        Some(Self::new(self.seqname, new_position))
    }
}

impl<R, P> Ord for GenomicPosition<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    fn cmp(
        &self,
        other: &Self,
    ) -> std::cmp::Ordering {
        self.partial_cmp(other)
            .expect("Cannot compare genomic positions with different seqnames")
    }
}

impl<R, P> Eq for GenomicPosition<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
}

impl<R, P> PartialEq for GenomicPosition<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.seqname.as_ref() == other.seqname.as_ref()
            && self.position == other.position
    }
}

#[allow(clippy::non_canonical_partial_ord_impl)]
impl<R, P> PartialOrd for GenomicPosition<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<std::cmp::Ordering> {
        if self.seqname.as_ref() == other.seqname.as_ref() {
            self.position
                .to_usize()
                .expect("Failed to convert position to usize")
                .partial_cmp(
                    &other
                        .position
                        .to_usize()
                        .expect("Failed to convert position to usize"),
                )
        }
        else {
            None
        }
    }
}

impl<R, P> Display for GenomicPosition<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        write!(
            f,
            "{}:{}",
            self.seqname.as_ref(),
            self.position
                .to_usize()
                .expect("Failed to convert position to usize")
        )
    }
}
