use std::fmt::Display;
use std::ops::Range;

use bio::bio_types::annot::loc::Loc;
use bio::bio_types::strand::ReqStrand;
use serde::{
    Deserialize,
    Serialize,
};

use super::GenomicPosition;
use crate::data_structs::enums::Strand;
use crate::data_structs::typedef::{
    BsxSmallStr,
    PosType,
    SeqNameStr,
};

/// Represents a contig with a sequence name, start position, end position, and
/// strand.
#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct Contig {
    seqname: BsxSmallStr,
    start:   PosType,
    end:     PosType,
    strand:  Strand,
}

impl Contig {
    /// Creates a new `Contig`.
    pub fn new(
        seqname: BsxSmallStr,
        start: PosType,
        end: PosType,
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
    pub fn start(&self) -> PosType {
        self.start
    }

    /// Returns the end position.
    pub fn end(&self) -> PosType {
        self.end
    }

    /// Returns the start position of the contig.
    pub fn start_gpos(&self) -> GenomicPosition {
        GenomicPosition::new(self.seqname.clone(), self.start)
    }

    /// Returns the end position of the contig.
    pub fn end_gpos(&self) -> GenomicPosition {
        GenomicPosition::new(self.seqname.clone(), self.end)
    }

    /// Returns the strand of the contig.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns the sequence name of the contig.
    pub fn seqname(&self) -> &BsxSmallStr {
        &self.seqname
    }

    /// Returns the length of the contig.
    pub fn length(&self) -> PosType {
        self.end - self.start
    }

    /// Extends the contig upstream by a given length.
    pub fn extend_upstream(
        &mut self,
        length: PosType,
    ) {
        self.start = self.start.saturating_sub(length);
    }

    /// Extends the contig downstream by a given length.
    pub fn extend_downstream(
        &mut self,
        length: PosType,
    ) {
        self.end = self.end.saturating_add(length);
    }

    /// Sets the start position of the contig.
    pub fn set_start(
        &mut self,
        start: PosType,
    ) {
        self.start = start;
    }

    /// Sets the end position of the contig.
    pub fn set_end(
        &mut self,
        end: PosType,
    ) {
        self.end = end;
    }

    /// Checks if this contig is fully contained within another contig.
    pub fn is_in(
        &self,
        other: &Self,
    ) -> bool {
        self.seqname == other.seqname
            && self.start >= other.start
            && self.end <= other.end
    }

    pub fn is_empty(&self) -> bool {
        self.start == self.end && self.start == 0 && self.seqname.as_str() == ""
    }
}

impl From<Range<GenomicPosition>> for Contig {
    /// Converts from a range of `GenomicPosition`s.
    fn from(value: Range<GenomicPosition>) -> Self {
        if value.start.seqname() != value.end.seqname() {
            panic!("Start and end positions must have the same sequence name")
        }
        Self {
            seqname: value.start.seqname().clone(),
            start:   value.start.position(),
            end:     value.end.position(),
            strand:  Strand::None,
        }
    }
}

impl From<Contig> for Range<GenomicPosition> {
    /// Converts into a range of `GenomicPosition`s.
    fn from(value: Contig) -> Self {
        value.start_gpos()..value.end_gpos()
    }
}

impl From<&Contig> for Range<GenomicPosition> {
    /// Converts into a range of `GenomicPosition`s.
    fn from(value: &Contig) -> Self {
        value.start_gpos()..value.end_gpos()
    }
}

impl From<bio::io::bed::Record> for Contig {
    /// Converts from a `bio::io::bed::Record`.
    fn from(value: bio::io::bed::Record) -> Self {
        Self {
            seqname: BsxSmallStr::from(value.chrom()),
            start:   value.start() as PosType,
            end:     value.end() as PosType,
            strand:  match value.strand() {
                Some(bio::bio_types::strand::Strand::Forward) => Strand::Forward,
                Some(bio::bio_types::strand::Strand::Reverse) => Strand::Reverse,
                Some(bio::bio_types::strand::Strand::Unknown) => Strand::None,
                None => Strand::None,
            },
        }
    }
}

impl From<Contig> for bio::io::bed::Record {
    /// Converts into a `bio::io::bed::Record`.
    fn from(value: Contig) -> Self {
        let mut record = bio::io::bed::Record::new();
        record.set_chrom(value.seqname.as_ref());
        record.set_start(value.start as u64);
        record.set_end(value.end as u64);
        record
    }
}

impl From<bio::io::gff::Record> for Contig {
    /// Converts from a `bio::io::gff::Record`.
    fn from(value: bio::io::gff::Record) -> Self {
        Self {
            seqname: BsxSmallStr::from(value.seqname()),
            start:   *value.start() as PosType,
            end:     *value.end() as PosType,
            strand:  match value.strand() {
                Some(bio::bio_types::strand::Strand::Forward) => Strand::Forward,
                Some(bio::bio_types::strand::Strand::Reverse) => Strand::Reverse,
                Some(bio::bio_types::strand::Strand::Unknown) => Strand::None,
                None => Strand::None,
            },
        }
    }
}

impl<R, S> From<bio::bio_types::annot::contig::Contig<R, S>> for Contig
where
    R: SeqNameStr,
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
            BsxSmallStr::from(value.refid().as_ref()),
            value.start() as PosType,
            (value.start() + value.length() as isize) as PosType,
            strand,
        )
    }
}

impl<R> From<Contig> for bio::bio_types::annot::contig::Contig<R, Option<ReqStrand>>
where
    R: SeqNameStr + From<String>,
{
    /// Converts into a `bio_types::annot::contig::Contig`.
    fn from(value: Contig) -> Self {
        let strand = match value.strand {
            Strand::Forward => Some(ReqStrand::Forward),
            Strand::Reverse => Some(ReqStrand::Reverse),
            Strand::None => None,
        };
        bio::bio_types::annot::contig::Contig::new(
            R::from(value.seqname.to_string()),
            value.start as isize,
            (value.length()) as usize,
            strand,
        )
    }
}

impl PartialOrd for Contig {
    /// Compares two `Contig`s.
    ///
    /// Returns `None` if the sequence names are different or regions intersect.
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<std::cmp::Ordering> {
        if self.seqname != other.seqname {
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

impl Display for Contig {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self.strand {
            Strand::None => write!(f, "{}:{}-{}", self.seqname, self.start, self.end),
            Strand::Forward => {
                write!(f, "{}:{}-{} (+)", self.seqname, self.start, self.end)
            },
            Strand::Reverse => {
                write!(f, "{}:{}-{} (-)", self.seqname, self.start, self.end)
            },
        }
    }
}
