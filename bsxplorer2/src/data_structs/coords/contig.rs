use std::fmt::Display;
use std::ops::Range;
use std::str::FromStr;

use bio::bio_types::annot::loc::Loc;
use bio::bio_types::strand::ReqStrand;
use serde::{Deserialize, Serialize};

use super::GenomicPosition;
use crate::data_structs::enums::Strand;
use crate::data_structs::typedef::{SeqNameStr, SeqPosNum};

/// Represents a contig with a sequence name, start position, end position, and
/// strand.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Contig<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum, {
    seqname: R,
    start:   P,
    end:     P,
    strand:  Strand,
}

impl<R, P> Contig<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    /// Creates a new `Contig`.
    pub fn new(
        seqname: R,
        start: P,
        end: P,
        strand: Strand,
    ) -> Self {
        assert!(
            start <= end,
            "Start position must be less than or equal to end position"
        );
        Self {
            seqname,
            start,
            end,
            strand,
        }
    }

    /// Returns the start position.
    pub fn start(&self) -> P {
        self.start
    }

    /// Returns the end position.
    pub fn end(&self) -> P {
        self.end
    }

    /// Returns the start position of the contig.
    pub fn start_gpos(&self) -> GenomicPosition<R, P> {
        GenomicPosition::new(self.seqname.clone(), self.start)
    }

    /// Returns the end position of the contig.
    pub fn end_gpos(&self) -> GenomicPosition<R, P> {
        GenomicPosition::new(self.seqname.clone(), self.end)
    }

    /// Returns the strand of the contig.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns the sequence name of the contig.
    pub fn seqname(&self) -> &R {
        &self.seqname
    }

    /// Returns the length of the contig.
    pub fn length(&self) -> P {
        self.end - self.start
    }

    /// Extends the contig upstream by a given length.
    pub fn extend_upstream(
        &mut self,
        length: P,
    ) {
        self.start = self.start.saturating_sub(length);
    }

    /// Extends the contig downstream by a given length.
    pub fn extend_downstream(
        &mut self,
        length: P,
    ) {
        self.end = self.end.saturating_add(length);
    }

    /// Sets the start position of the contig.
    pub fn set_start(
        &mut self,
        start: P,
    ) {
        self.start = start;
    }

    /// Sets the end position of the contig.
    pub fn set_end(
        &mut self,
        end: P,
    ) {
        self.end = end;
    }

    /// Checks if this contig is fully contained within another contig.
    pub fn is_in(
        &self,
        other: &Self,
    ) -> bool {
        self.seqname.as_ref() == other.seqname.as_ref()
            && self.start >= other.start
            && self.end <= other.end
    }

    /// Casts the contig to a new type.
    pub fn cast<R2: SeqNameStr, P2: SeqPosNum>(
        self,
        seqname_fn: fn(R) -> R2,
        pos_fn: fn(P) -> P2,
    ) -> Contig<R2, P2> {
        Contig {
            seqname: seqname_fn(self.seqname),
            start:   pos_fn(self.start),
            end:     pos_fn(self.end),
            strand:  self.strand,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.start == self.end
            && self.start == P::zero()
            && self.seqname.as_ref() == ""
    }
}

impl<R, P> From<Range<GenomicPosition<R, P>>> for Contig<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    /// Converts from a range of `GenomicPosition`s.
    fn from(value: Range<GenomicPosition<R, P>>) -> Self {
        assert_eq!(
            value.start.seqname().as_ref(),
            value.end.seqname().as_ref(),
            "Start and end positions must have the same sequence name"
        );
        Self {
            seqname: value.start.seqname(),
            start:   value.start.position(),
            end:     value.end.position(),
            strand:  Strand::None,
        }
    }
}

impl<R, P> From<Contig<R, P>> for Range<GenomicPosition<R, P>>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    /// Converts into a range of `GenomicPosition`s.
    fn from(value: Contig<R, P>) -> Self {
        value.start_gpos()..value.end_gpos()
    }
}

impl<S, P> From<bio::io::bed::Record> for Contig<S, P>
where
    S: SeqNameStr + FromStr,
    P: SeqPosNum,
    <S as FromStr>::Err: std::fmt::Debug,
{
    /// Converts from a `bio::io::bed::Record`.
    fn from(value: bio::io::bed::Record) -> Self {
        Self {
            seqname: S::from_str(value.chrom()).unwrap(),
            start:   P::from(value.start()).expect("Failed to convert start to P"),
            end:     P::from(value.end()).expect("Failed to convert end to P"),
            strand:  match value.strand() {
                Some(bio::bio_types::strand::Strand::Forward) => Strand::Forward,
                Some(bio::bio_types::strand::Strand::Reverse) => Strand::Reverse,
                Some(bio::bio_types::strand::Strand::Unknown) => Strand::None,
                None => Strand::None,
            },
        }
    }
}

impl<R, P> From<Contig<R, P>> for bio::io::bed::Record
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    /// Converts into a `bio::io::bed::Record`.
    fn from(value: Contig<R, P>) -> Self {
        let mut record = bio::io::bed::Record::new();
        record.set_chrom(value.seqname.as_ref());
        record.set_start(
            value
                .start
                .to_u64()
                .expect("Failed to convert start to u64"),
        );
        record.set_end(value.end.to_u64().expect("Failed to convert end to u64"));
        record
    }
}

impl<S, P> From<bio::io::gff::Record> for Contig<S, P>
where
    S: SeqNameStr + FromStr,
    P: SeqPosNum,
    <S as FromStr>::Err: std::fmt::Debug,
{
    /// Converts from a `bio::io::gff::Record`.
    fn from(value: bio::io::gff::Record) -> Self {
        Self {
            seqname: S::from_str(value.seqname()).unwrap(),
            start:   P::from(value.start().to_owned())
                .expect("Failed to convert start to P"),
            end:     P::from(value.end().to_owned())
                .expect("Failed to convert end to P"),
            strand:  match value.strand() {
                Some(bio::bio_types::strand::Strand::Forward) => Strand::Forward,
                Some(bio::bio_types::strand::Strand::Reverse) => Strand::Reverse,
                Some(bio::bio_types::strand::Strand::Unknown) => Strand::None,
                None => Strand::None,
            },
        }
    }
}

impl<R, P, S> From<bio::bio_types::annot::contig::Contig<R, S>> for Contig<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
    S: Into<Option<ReqStrand>> + Copy,
{
    /// Converts from a `bio_types::annot::contig::Contig`.
    fn from(value: bio::bio_types::annot::contig::Contig<R, S>) -> Self {
        let s: Option<ReqStrand> = value.strand().into();
        let strand = match s {
            Some(ReqStrand::Forward) => Strand::Forward,
            Some(ReqStrand::Reverse) => Strand::Reverse,
            None => Strand::None,
        };
        Self::new(
            value.refid().to_owned(),
            P::from(value.start()).expect("Failed to convert start to P"),
            P::from(value.start() + value.length() as isize)
                .expect("Failed to convert position to P"),
            strand,
        )
    }
}

impl<R, P> From<Contig<R, P>>
    for bio::bio_types::annot::contig::Contig<R, Option<ReqStrand>>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    /// Converts into a `bio_types::annot::contig::Contig`.
    fn from(value: Contig<R, P>) -> Self {
        let strand = match value.strand {
            Strand::Forward => Some(ReqStrand::Forward),
            Strand::Reverse => Some(ReqStrand::Reverse),
            Strand::None => None,
        };
        bio::bio_types::annot::contig::Contig::new(
            value.seqname.to_owned(),
            value
                .start
                .to_isize()
                .expect("Failed to convert start to isize"),
            value
                .length()
                .to_usize()
                .expect("Failed to convert length to usize"),
            strand,
        )
    }
}

impl<R, P> PartialOrd for Contig<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    /// Compares two `Contig`s.
    ///
    /// Returns `None` if the sequence names are different or regions intersect.
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<std::cmp::Ordering> {
        if self.seqname.as_ref() != other.seqname.as_ref() {
            return None;
        }
        if self.start >= other.end {
            return Some(std::cmp::Ordering::Greater);
        }
        if self.end <= other.start {
            return Some(std::cmp::Ordering::Less);
        }
        None
    }
}

impl<R, P> Eq for Contig<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
}

impl<R, P> PartialEq for Contig<R, P>
where
    R: SeqNameStr,
    P: SeqPosNum,
{
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.seqname.as_ref() == other.seqname.as_ref()
            && self.start == other.start
            && self.end == other.end
            && self.strand == other.strand
    }
}

impl<R, P> Display for Contig<R, P>
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
            "{}:{}-{} ({})",
            self.seqname.as_ref(),
            self.start
                .to_usize()
                .expect("Failed to convert start to usize"),
            self.end.to_usize().expect("Failed to convert end to usize"),
            self.strand
        )
    }
}
