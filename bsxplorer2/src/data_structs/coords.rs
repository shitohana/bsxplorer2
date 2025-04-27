use std::{fmt::Display, ops::{Add, Range, Sub}};

use bio_types::{annot::loc::Loc, strand::{NoStrand, ReqStrand}};
use num::{PrimInt, Unsigned};

use crate::data_structs::enums::Strand;



/// Represents a genomic position with a sequence name and a position.
#[derive(Debug, Clone)]
pub struct GenomicPosition<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    seqname: R,
    position: P
}

impl<P> From<bio_types::annot::pos::SeqPosUnstranded> for GenomicPosition<String, P>
where 
    P: Unsigned + PrimInt
{
    fn from(value: bio_types::annot::pos::SeqPosUnstranded) -> Self {
        Self {
            seqname: value.refid().to_owned(),
            position: P::from(value.pos()).expect("Failed to convert position to P")
        }
    }
}

impl<R, P> From<GenomicPosition<R, P>> for bio_types::annot::pos::SeqPosUnstranded
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    fn from(value: GenomicPosition<R, P>) -> Self {
        bio_types::annot::pos::SeqPosUnstranded::new(
            value.seqname.as_ref().to_string(), 
            value.position.to_isize().expect("Failed to convert position to isize"), 
            NoStrand::Unknown
        )
    }
}


impl<R, P> Sub for GenomicPosition<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    type Output = Option<Self>;

    /// Subtracts two `GenomicPosition`s.
    ///
    /// Returns `None` if the sequence names are different or if the right-hand side position is greater than the left-hand side position.
    fn sub(self, rhs: Self) -> Self::Output {
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
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    type Output = Option<Self>;

    /// Adds two `GenomicPosition`s.
    ///
    /// Returns `None` if the sequence names are different.
    fn add(self, rhs: Self) -> Self::Output {
        if self.seqname.as_ref() != rhs.seqname.as_ref() {
            return None;
        }
        let new_position = self.position + rhs.position;
        Some(Self::new(self.seqname, new_position))
    }
}

impl<R, P> GenomicPosition<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    /// Creates a new `GenomicPosition`.
    pub fn new(seqname: R, position: P) -> Self {
        Self { seqname, position }
    }
    
    pub fn seqname(&self) -> R {
        self.seqname.clone()
    }
    
    pub fn position(&self) -> P {
        self.position
    }
}

impl<R, P> Ord for GenomicPosition<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).expect("Cannot compare genomic positions with different seqnames")
    }
}

impl<R, P> Eq for GenomicPosition<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
}

impl<R, P> PartialEq for GenomicPosition<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    fn eq(&self, other: &Self) -> bool {
        self.seqname.as_ref() == other.seqname.as_ref() && self.position == other.position
    }
}


impl<R, P> PartialOrd for GenomicPosition<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.seqname.as_ref() == other.seqname.as_ref() {
            self.position.to_usize().expect("Failed to convert position to usize").partial_cmp(&other.position.to_usize().expect("Failed to convert position to usize"))
        } else {
            None
        }
    }
}

impl<R, P> Display for GenomicPosition<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}",
            self.seqname.as_ref(),
            self.position.to_usize().expect("Failed to convert position to usize")
        )
    }
}

/// Represents a contig with a sequence name, start position, end position, and strand.
#[derive(Debug, Clone)]
pub struct Contig<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    seqname: R,
    start: P,
    end: P,
    strand: Strand
}

impl<R, P> From<Range<GenomicPosition<R, P>>> for Contig<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    fn from(value: Range<GenomicPosition<R, P>>) -> Self {
        assert!(value.start.seqname.as_ref() == value.end.seqname.as_ref(), "Start and end positions must have the same sequence name");
        Self {
            seqname: value.start.seqname,
            start: value.start.position,
            end: value.end.position,
            strand: Strand::None
        }
    }
}

impl<R, P> From<Contig<R, P>> for Range<GenomicPosition<R, P>>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    fn from(value: Contig<R, P>) -> Self {
        value.start()..value.end()
    }
}

impl<P> From<bio::io::bed::Record> for Contig<String, P>
where 
    P: Unsigned + PrimInt
{
    fn from(value: bio::io::bed::Record) -> Self {
        Self {
            seqname: value.chrom().to_owned(),
            start: P::from(value.start()).expect("Failed to convert start to P"),
            end: P::from(value.end()).expect("Failed to convert end to P"),
            strand: match value.strand() {
                Some(bio_types::strand::Strand::Forward) => Strand::Forward,
                Some(bio_types::strand::Strand::Reverse) => Strand::Reverse,
                Some(bio_types::strand::Strand::Unknown) => Strand::None,
                None => Strand::None,
            }
        }
    }
}

impl<R, P> From<Contig<R, P>> for bio::io::bed::Record
where
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt,
{
    fn from(value: Contig<R, P>) -> Self {
        let strand = match value.strand {
            Strand::Forward => Some(bio_types::strand::Strand::Forward),
            Strand::Reverse => Some(bio_types::strand::Strand::Reverse),
            Strand::None => None,
        };
        let mut record = bio::io::bed::Record::new();
        record.set_chrom(value.seqname.as_ref());
        record.set_start(value.start.to_u64().expect("Failed to convert start to u64"));
        record.set_end(value.end.to_u64().expect("Failed to convert end to u64"));
        record
    }
}

impl<P> From<bio::io::gff::Record> for Contig<String, P>
where 
    P: Unsigned + PrimInt
{
    fn from(value: bio::io::gff::Record) -> Self {
        Self {
            seqname: value.seqname().to_owned(),
            start: P::from(value.start().to_owned()).expect("Failed to convert start to P"),
            end: P::from(value.end().to_owned()).expect("Failed to convert end to P"),
            strand: match value.strand() {
                Some(bio_types::strand::Strand::Forward) => Strand::Forward,
                Some(bio_types::strand::Strand::Reverse) => Strand::Reverse,
                Some(bio_types::strand::Strand::Unknown) => Strand::None,
                None => Strand::None,
            }
        }
    }
}

impl<R, P, S> From<bio_types::annot::contig::Contig<R, S>> for Contig<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt,
    S: Into<Option<ReqStrand>> + Copy
{
    fn from(value: bio_types::annot::contig::Contig<R, S>) -> Self {
        let s: Option<ReqStrand> = value.strand().into();
        let strand = match s {
            Some(bio_types::strand::ReqStrand::Forward) => Strand::Forward,
            Some(bio_types::strand::ReqStrand::Reverse) => Strand::Reverse,
            None => Strand::None,
        };
        Self::new(
            value.refid().to_owned(), 
            P::from(value.start()).expect("Failed to convert start to P"), 
            P::from(value.start() + value.length() as isize).expect("Failed to convert position to P"), 
            strand
        )
    }
}

impl<R, P> From<Contig<R, P>> for bio_types::annot::contig::Contig<R, Option<ReqStrand>>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt,
{
    fn from(value: Contig<R, P>) -> Self {
        let strand  = match value.strand {
            Strand::Forward => Some(bio_types::strand::ReqStrand::Forward),
            Strand::Reverse => Some(bio_types::strand::ReqStrand::Reverse),
            Strand::None => None,
        };
        bio_types::annot::contig::Contig::new(
            value.seqname.to_owned(), 
            value.start.to_isize().expect("Failed to convert start to isize"), 
            value.length().to_usize().expect("Failed to convert length to usize"), 
            strand
        )
    }
}


impl<R, P> Contig<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    /// Creates a new `Contig`.
    pub fn new(seqname: R, start: P, end: P, strand: Strand) -> Self {
        assert!(start <= end, "Start position must be less than or equal to end position");
        Self { seqname, start, end, strand }
    }

    /// Returns the start position of the contig.
    pub fn start(&self) -> GenomicPosition<R, P> {
        GenomicPosition::new(self.seqname.clone(), self.start)
    }

    /// Returns the end position of the contig.
    pub fn end(&self) -> GenomicPosition<R, P> {
        GenomicPosition::new(self.seqname.clone(), self.end)
    }

    /// Returns the strand of the contig.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns the sequence name of the contig.
    pub fn seqname(&self) -> R {
        self.seqname.clone()
    }

    /// Returns the length of the contig.
    pub fn length(&self) -> P {
        self.end - self.start
    }

    pub fn extend_upstream(&mut self, length: P) {
        self.start = self.start.saturating_sub(length);
    }
    pub fn extend_downstream(&mut self, length: P) {
        self.end = self.end.saturating_add(length);
    }
}


impl<R, P> PartialOrd for Contig<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    /// Compares two `Contig`s.
    ///
    /// Returns `None` if the sequence names are different.
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
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
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
}

impl<R, P> PartialEq for Contig<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    fn eq(&self, other: &Self) -> bool {
        self.seqname.as_ref() == other.seqname.as_ref() && self.start == other.start && self.end == other.end && self.strand == other.strand
    }
}

impl<R, P> Display for Contig<R, P>
where 
    R: AsRef<str> + Clone,
    P: Unsigned + PrimInt
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}-{} ({})",
            self.seqname.as_ref(),
            self.start.to_usize().expect("Failed to convert start to usize"),
            self.end.to_usize().expect("Failed to convert end to usize"),
            self.strand
        )
    }
}