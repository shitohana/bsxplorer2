use log::warn;
use std::cmp::Ordering;
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::ops::{Add, Range, Shr, Sub};
use num::Unsigned;
use polars::export::num::PrimInt;

#[derive(Clone, Hash, PartialEq, Eq, Debug)]
pub struct RegionCoordinates<N>
where 
    N: PrimInt + Unsigned + Clone 
{
    pub(crate) chr: String,
    pub(crate) start: N,
    pub(crate) end: N,
}

#[derive(Debug, Clone)]
pub struct GenomicPosition<N>
where
    N: PrimInt + Unsigned + Clone
{
    chr: String,
    position: N,
}

impl<N> Display for GenomicPosition<N> where
    N: PrimInt + Unsigned + Clone + Display
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}", self.chr, self.position)
    }
}

impl<N> GenomicPosition<N>
where
    N: PrimInt + Unsigned + Clone
{
    pub fn new(chr: String, position: N) -> GenomicPosition<N> {
        GenomicPosition { chr, position }
    }

    pub fn chr(&self) -> &str {
        &self.chr
    }

    pub fn position(&self) -> N {
        self.position
    }
}

impl<N> Add<N> for GenomicPosition<N> where
    N: PrimInt + Unsigned + Clone
{
    type Output = GenomicPosition<N>;
    fn add(self, rhs: N) -> Self::Output {
        GenomicPosition::new(self.chr, self.position + rhs)
    }
}

impl<N> Sub<N> for GenomicPosition<N> 
where
    N: PrimInt + Unsigned + Clone
{
    type Output = GenomicPosition<N>;
    fn sub(self, rhs: N) -> Self::Output {
        GenomicPosition::new(self.chr, self.position - rhs)
    }
}

impl<N> PartialEq for GenomicPosition<N> where
    N: PrimInt + Unsigned + Clone
{
    fn eq(&self, other: &Self) -> bool {
        self.chr == other.chr && self.position == other.position
    }
}

impl<N> PartialOrd for GenomicPosition<N>
where
    N: PrimInt + Unsigned + Clone
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.chr == other.chr {
            Some(self.position.cmp(&other.position))
        } else {
            None
        }
    }
}

impl<N> Shr for GenomicPosition<N>
where
    N: PrimInt + Unsigned + Clone + Display
{
    type Output = Option<RegionCoordinates<N>>;

    fn shr(self, rhs: Self) -> Self::Output {
        let chr_cmp = self.partial_cmp(&rhs);
        if chr_cmp.is_some() && chr_cmp.unwrap() == Ordering::Greater {
            warn!(
                "Bad range operation for {} and {} regions. Lhs must be less than rhs.",
                self, rhs
            );
            Some(RegionCoordinates::new(
                self.chr,
                rhs.position,
                self.position,
            ))
        } else if chr_cmp.is_some() && chr_cmp.unwrap() == Ordering::Less {
            Some(RegionCoordinates::new(
                self.chr,
                self.position,
                rhs.position,
            ))
        } else {
            None
        }
    }
}

impl<N> RegionCoordinates<N> 
where
    N: PrimInt + Unsigned + Clone + Display
{
    pub fn new(chr: String, start: N, end: N) -> Self {
        assert!(
            start <= end,
            "End position can not be less than start! {}:{}-{}",
            chr,
            start,
            end
        );
        Self { chr, start, end }
    }
    pub fn try_from_position(start: GenomicPosition<N>, end: GenomicPosition<N>) -> Option<Self> {
        start >> end
    }

    pub fn into_positions(self) -> (GenomicPosition<N>, GenomicPosition<N>) {
        (
            GenomicPosition::new(self.chr.clone(), self.start),
            GenomicPosition::new(self.chr, self.end),
        )
    }
    pub fn chr(&self) -> &str {
        self.chr.as_str()
    }
    pub fn start(&self) -> N {
        self.start
    }
    pub fn end(&self) -> N {
        self.end
    }
    pub fn start_gpos(&self) -> GenomicPosition<N> {
        GenomicPosition::new(self.chr.clone(), self.start)
    }
    pub fn end_gpos(&self) -> GenomicPosition<N> {
        GenomicPosition::new(self.chr.clone(), self.end)
    }
    pub fn length(&self) -> N {
        self.end - self.start
    }

    pub fn intersect(&self, other: &Self) -> Range<N> {
        // Determine longer
        let (lhs, rhs) = if self.length() < other.length() {
            (self, other)
        } else {
            (other, self)
        };

        // Fully covers
        // rhs | |#####|------->
        // lhs | |-----|
        // OR
        // rhs | <-|#####|----->
        // lhs |   |-----|
        // OR
        // rhs | <-------|#####|
        // lhs |         |-----|
        if rhs.start <= lhs.start && rhs.end >= lhs.end {
            rhs.start..rhs.end + N::from(1).unwrap()
        }
        // rhs |      |###|------>
        // lhs | <----|---|
        // OR
        // rhs |      #---------->
        // lhs | <----|
        else if rhs.start > lhs.start && rhs.end >= lhs.start {
            rhs.start..lhs.end + N::from(1).unwrap()
        }
        // rhs | <-------|##|
        // lhs |         |------->
        // OR
        // rhs | <----------#
        // lhs |            |---->
        else if rhs.end >= lhs.start && rhs.end > lhs.end {
            lhs.start..rhs.end + N::from(1).unwrap()
        } else {
            N::from(0).unwrap()..N::from(0).unwrap()
        }
    }

    pub fn set_start(&mut self, start: N) -> Result<(), Box<dyn Error>> {
        if start < self.end {
            self.start = start;
            Ok(())
        } else {
            Err(Box::from(format!(
                "Could not modify start value. Provided {} start value is greater than end {}.",
                start, self.end
            )))
        }
    }
    pub fn set_end(&mut self, end: N) -> Result<(), Box<dyn Error>> {
        if end < self.start {
            self.end = end;
            Ok(())
        } else {
            Err(Box::from(format!(
                "Could not modify end value. Provided {} end value is less than start {}.",
                end, self.start
            )))
        }
    }
    pub fn set_chromosome(&mut self, chr: String) {
        self.chr = chr;
    }

    /// Expands limits of the regions by specified value. If length of
    /// the resulting region is more than max_length, [Err] is returned.
    /// Also, if start is less than value [Err] is returned.
    /// Otherwise bounds of the region will be modified themselves
    pub fn expand(mut self, value: N, max_length: Option<N>) -> Self {
        let max_length = max_length.unwrap_or(N::max_value());

        if self.end + value > max_length {
            self.end = max_length;
        } else {
            self.end = self.end + value
        }
        if self.start > value {
            self.start = self.start - value;
        }
        self
    }
}

impl<N> Display for RegionCoordinates<N> where
    N: PrimInt + Unsigned + Clone + Display
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chr, self.start, self.end)
    }
}

#[cfg(test)]
mod tests {
    // TODO
}
