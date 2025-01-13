use log::warn;
use std::cmp::Ordering;
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::ops::{Add, Range, Shr, Sub};

#[derive(Clone, Hash, PartialEq, Eq, Debug)]
pub struct RegionCoordinates {
    pub(crate) chr: String,
    pub(crate) start: u32,
    pub(crate) end: u32,
}

#[derive(Debug, Clone)]
pub struct GenomicPosition {
    chr: String,
    position: u32,
}

impl Display for GenomicPosition {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}", self.chr, self.position)
    }
}

impl GenomicPosition {
    pub fn new(chr: String, position: u32) -> GenomicPosition {
        GenomicPosition { chr, position }
    }

    pub fn chr(&self) -> &str {
        &self.chr
    }

    pub fn position(&self) -> u32 {
        self.position
    }
}

impl Add<u32> for GenomicPosition {
    type Output = GenomicPosition;
    fn add(self, rhs: u32) -> Self::Output {
        GenomicPosition::new(self.chr, self.position + rhs)
    }
}

impl Sub<u32> for GenomicPosition {
    type Output = GenomicPosition;
    fn sub(self, rhs: u32) -> Self::Output {
        GenomicPosition::new(self.chr, self.position - rhs)
    }
}

impl PartialEq for GenomicPosition {
    fn eq(&self, other: &Self) -> bool {
        self.chr == other.chr && self.position == other.position
    }
}

impl PartialOrd for GenomicPosition {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.chr == other.chr {
            Some(self.position.cmp(&other.position))
        } else {
            None
        }
    }
}

impl Shr for GenomicPosition {
    type Output = Option<RegionCoordinates>;

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

impl RegionCoordinates {
    pub fn new(chr: String, start: u32, end: u32) -> Self {
        assert!(
            start <= end,
            "End position can not be less than start! {}:{}-{}",
            chr,
            start,
            end
        );
        Self { chr, start, end }
    }
    pub fn try_from_position(start: GenomicPosition, end: GenomicPosition) -> Option<Self> {
        start >> end
    }

    pub fn into_positions(self) -> (GenomicPosition, GenomicPosition) {
        (
            GenomicPosition::new(self.chr.clone(), self.start),
            GenomicPosition::new(self.chr, self.end),
        )
    }
    pub fn chr(&self) -> &str {
        self.chr.as_str()
    }
    pub fn start(&self) -> u32 {
        self.start
    }
    pub fn end(&self) -> u32 {
        self.end
    }
    pub fn start_gpos(&self) -> GenomicPosition {
        GenomicPosition::new(self.chr.clone(), self.start)
    }
    pub fn end_gpos(&self) -> GenomicPosition {
        GenomicPosition::new(self.chr.clone(), self.end)
    }
    pub fn length(&self) -> u32 {
        self.end - self.start
    }

    pub fn intersect(&self, other: &Self) -> Range<u32> {
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
            rhs.start..rhs.end + 1
        }
        // rhs |      |###|------>
        // lhs | <----|---|
        // OR
        // rhs |      #---------->
        // lhs | <----|
        else if rhs.start > lhs.start && rhs.end >= lhs.start {
            rhs.start..lhs.end + 1
        }
        // rhs | <-------|##|
        // lhs |         |------->
        // OR
        // rhs | <----------#
        // lhs |            |---->
        else if rhs.end >= lhs.start && rhs.end > lhs.end {
            lhs.start..rhs.end + 1
        } else {
            0..0
        }
    }

    pub fn set_start(&mut self, start: u32) -> Result<(), Box<dyn Error>> {
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
    pub fn set_end(&mut self, end: u32) -> Result<(), Box<dyn Error>> {
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
    pub fn expand(mut self, value: u32, max_length: Option<u32>) -> Self {
        let max_length = max_length.unwrap_or(u32::MAX);

        if self.end + value > max_length {
            self.end = max_length;
        } else {
            self.end += value
        }
        if self.start > value {
            self.start -= value;
        }
        self
    }
}

impl Display for RegionCoordinates {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chr, self.start, self.end)
    }
}

#[cfg(test)]
mod tests {
    // TODO
}
