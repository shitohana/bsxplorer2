use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::io::{BufRead, Read, Seek};

use bio::io::fasta::*;
use hashbrown::HashMap;
use log::{debug, info};
use num::{PrimInt, ToPrimitive, Unsigned};

use crate::data_structs::region::{GenomicPosition, RegionCoordinates};
use crate::utils::types::{PosNum, Strand};

/// A wrapper around bio's IndexedReader to handle FASTA file operations
pub struct FastaReader<R>
where
    R: Read + BufRead + Seek, {
    reader: IndexedReader<R>,
}

impl<R> FastaReader<R>
where
    R: Read + BufRead + Seek,
{
    /// Create a new FastaReader from an existing IndexedReader
    pub fn new(reader: IndexedReader<R>) -> Self {
        info!("Fasta reader initialized");
        Self { reader }
    }

    /// Try to create a FastaReader from file handles for FASTA and FAI files
    pub fn try_from_handle(
        fasta_handle: R,
        fai_handle: R,
    ) -> Result<Self, Box<dyn Error>> {
        let reader = IndexedReader::new(fasta_handle, fai_handle)?;
        Ok(Self::new(reader))
    }

    /// Get the underlying FASTA index
    pub fn index(&self) -> &Index { &self.reader.index }

    /// Fetch a genomic region specified by coordinates
    pub fn fetch_region<N>(
        &mut self,
        region_coordinates: RegionCoordinates<N>,
    ) -> Result<Vec<u8>, Box<dyn Error>>
    where
        N: PosNum, {
        let chr = region_coordinates.chr();
        let end_pos = region_coordinates
            .end()
            .to_u64()
            .unwrap_or_default();

        let seq_data = self
            .get_record_metadata(chr)
            .ok_or_else(|| {
                CoverageError::ChromosomeNotFound(chr.to_string())
            })?;

        if seq_data.len < end_pos {
            return Err(CoverageError::OutOfBounds(
                chr.to_string(),
                region_coordinates
                    .end()
                    .to_usize()
                    .unwrap_or_default(),
            )
            .into());
        }

        // Determine fetch end position, handling edge case when end is at
        // sequence length
        let fetch_end =
            if region_coordinates.end() != N::from(seq_data.len).unwrap() {
                end_pos + 1
            }
            else {
                seq_data.len
            };

        self.reader.fetch(
            chr,
            region_coordinates
                .start()
                .to_u64()
                .unwrap_or_default(),
            fetch_end,
        )?;

        let capacity = region_coordinates
            .length()
            .to_usize()
            .unwrap_or_default();
        let mut buffer = Vec::with_capacity(capacity);
        self.reader.read(&mut buffer)?;

        debug!("Fetched sequence for {}", region_coordinates);
        Ok(buffer)
    }

    /// Fetch an entire record by name
    pub fn fetch_record(
        &mut self,
        name: &str,
    ) -> Result<Vec<u8>, Box<dyn Error>> {
        let seq_data = self
            .get_record_metadata(name)
            .ok_or_else(|| {
                CoverageError::ChromosomeNotFound(name.to_string())
            })?;

        let mut buffer = Vec::with_capacity(seq_data.len as usize);
        self.reader.fetch_all(name)?;
        self.reader.read(&mut buffer)?;

        Ok(buffer)
    }

    /// Get metadata for a specific record
    pub fn get_record_metadata(
        &self,
        name: &str,
    ) -> Option<Sequence> {
        self.reader
            .index
            .sequences()
            .into_iter()
            .find(|s| s.name == name)
    }

    /// Get all record names in the FASTA file
    pub fn get_record_names(&self) -> Vec<String> {
        self.reader
            .index
            .sequences()
            .into_iter()
            .map(|s| s.name)
            .collect()
    }
}

/// Extends FastaReader with tracking of how much of each sequence has been read
pub struct FastaCoverageReader<R, V>
where
    R: Read + BufRead + Seek,
    V: PrimInt + Unsigned, {
    inner:    FastaReader<R>,
    coverage: Coverage<V>,
}

impl<R, V> From<FastaReader<R>> for FastaCoverageReader<R, V>
where
    R: Read + BufRead + Seek,
    V: PrimInt + Unsigned + ToPrimitive,
{
    fn from(fasta_reader: FastaReader<R>) -> Self {
        let coverage = Coverage::new(fasta_reader.reader.index.sequences());
        Self {
            inner: fasta_reader,
            coverage,
        }
    }
}

impl<R, V> FastaCoverageReader<R, V>
where
    R: Read + BufRead + Seek,
    V: PosNum,
{
    /// Create a new FastaCoverageReader from a FastaReader
    pub fn new(fasta_reader: FastaReader<R>) -> Self {
        Self::from(fasta_reader)
    }

    /// Get the current position for a given sequence name
    pub fn tell_pos(
        &self,
        name: &str,
    ) -> Option<V> {
        self.coverage.get(name).map(|c| c.read)
    }

    /// Get a reference to the inner FastaReader
    pub fn inner(&self) -> &FastaReader<R> { &self.inner }

    /// Get a mutable reference to the inner FastaReader
    pub(crate) fn inner_mut(&mut self) -> &mut FastaReader<R> {
        &mut self.inner
    }

    /// Get a reference to the coverage tracker
    pub fn coverage(&self) -> &Coverage<V> { &self.coverage }

    /// Get a mutable reference to the coverage tracker
    pub(crate) fn coverage_mut(&mut self) -> &mut Coverage<V> {
        &mut self.coverage
    }

    /// Get sequence up to a specified genomic position
    pub fn get_seq_until<N: PosNum>(
        &mut self,
        position: &GenomicPosition<N>,
    ) -> Result<Vec<u8>, Box<dyn Error>> {
        let chr = position.chr();
        let pos = position.position();

        // Convert position to the right numeric type and update coverage
        let pos_v = V::from(pos).ok_or_else(|| {
            CoverageError::ValueError(format!(
                "Failed to convert position {} to required numeric type",
                pos
            ))
        })?;

        let (fetch_start, fetch_end) = self.coverage.shift_to(chr, pos_v)?;

        // Create region coordinates for fetching
        let read_region = RegionCoordinates::new(
            chr.to_string(),
            fetch_start.to_u32().unwrap_or_default(),
            fetch_end.to_u32().unwrap_or_default(),
            Strand::None,
        );

        self.inner.fetch_region(read_region)
    }

    /// Get the remaining sequence for a chromosome
    pub fn get_seq_leftover(
        &mut self,
        chr: &str,
    ) -> Result<Vec<u8>, Box<dyn Error>> {
        let coverage_row = self.coverage.get(chr).ok_or_else(|| {
            CoverageError::ChromosomeNotFound(chr.to_string())
        })?;

        let end_pos = coverage_row
            .total
            .to_u32()
            .unwrap_or_default();
        let end_gpos = GenomicPosition::new(chr.to_string(), end_pos);

        self.get_seq_until(&end_gpos)
    }
}

/// Tracks position in a single sequence
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct CoverageRow<V>
where
    V: Unsigned + PrimInt, {
    read:  V, // Position up to which the sequence has been read
    total: V, // Total length of the sequence
}

impl<V> CoverageRow<V>
where
    V: Unsigned + PrimInt,
{
    /// Create a new CoverageRow with specified read and total positions
    pub fn new(
        read: V,
        total: V,
    ) -> Self {
        Self { read, total }
    }

    /// Create an empty CoverageRow with only the total length
    pub fn empty(total: V) -> Self {
        Self {
            read: V::zero(),
            total,
        }
    }

    /// Update the read position
    pub fn shift_to(
        &mut self,
        pos: V,
    ) -> Result<(V, V), CoverageError> {
        if self.read > pos {
            return Err(CoverageError::ValueError(format!(
                "New position {} is less than previous position {}",
                pos.to_u64().unwrap_or_default(),
                self.read.to_u64().unwrap_or_default()
            )));
        }

        // Previous position is one after the last read position
        let prev = self.read + V::one();

        // Don't exceed total length
        self.read = std::cmp::min(pos, self.total);

        Ok((prev, self.read))
    }

    /// Check if the entire sequence has been read
    pub fn is_finished(&self) -> bool { self.total == self.read }

    /// Get the current read position
    pub fn read(&self) -> V { self.read }

    /// Get the total sequence length
    pub fn total(&self) -> V { self.total }
}

/// Tracks coverage across multiple sequences
pub struct Coverage<V>
where
    V: Unsigned + PrimInt, {
    mapping: HashMap<String, CoverageRow<V>>,
    order:   Vec<String>, // Preserves the order of sequences
}

/// Errors that can occur during coverage operations
#[derive(Debug)]
pub enum CoverageError {
    ChromosomeNotFound(String),
    PreviousNotFinished(String),
    OutOfBounds(String, usize),
    ValueError(String),
}

impl Display for CoverageError {
    fn fmt(
        &self,
        f: &mut Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            CoverageError::ChromosomeNotFound(name) => {
                write!(f, "Chromosome '{}' not found", name)
            },
            CoverageError::PreviousNotFinished(name) => {
                write!(f, "Previous sequence '{}' not finished", name)
            },
            CoverageError::OutOfBounds(name, size) => {
                write!(f, "Position {} for '{}' is out of bounds", size, name)
            },
            CoverageError::ValueError(desc) => {
                write!(f, "Value error: {}", desc)
            },
        }
    }
}

impl Error for CoverageError {}

impl<V> Coverage<V>
where
    V: Unsigned + PrimInt,
{
    /// Create a new Coverage tracker from sequence information
    fn new(sequences: Vec<Sequence>) -> Self {
        let mapping = sequences
            .iter()
            .map(|seq| {
                let total = V::from(seq.len).unwrap_or_else(|| V::zero());
                (seq.name.clone(), CoverageRow::empty(total))
            })
            .collect();

        let order = sequences
            .into_iter()
            .map(|seq| seq.name)
            .collect();

        Self { mapping, order }
    }

    /// Get coverage information for a sequence by name
    pub fn get(
        &self,
        name: &str,
    ) -> Option<&CoverageRow<V>> {
        self.mapping.get(name)
    }

    /// Find a sequence's position in the order
    pub fn position(
        &self,
        name: &str,
    ) -> Option<usize> {
        self.order
            .iter()
            .position(|n| n == name)
    }

    /// Update the coverage position for a sequence
    pub fn shift_to(
        &mut self,
        name: &str,
        pos: V,
    ) -> Result<(V, V), CoverageError> {
        // Check if all previous sequences are finished
        if let Some(idx) = self.position(name) {
            if idx > 0 {
                if let Some(unfinished) = self.order[..idx].iter().find(|&n| {
                    !self
                        .mapping
                        .get(n)
                        .unwrap()
                        .is_finished()
                }) {
                    return Err(CoverageError::PreviousNotFinished(
                        unfinished.to_string(),
                    ));
                }
            }
        }

        // Update the coverage for this sequence
        self.mapping
            .get_mut(name)
            .ok_or_else(|| CoverageError::ChromosomeNotFound(name.to_string()))
            .and_then(|coverage| coverage.shift_to(pos))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coverage_row() {
        // Test finished
        assert!(CoverageRow::new(1u32, 1u32).is_finished());

        // Test shift without overflow
        let mut first = CoverageRow::new(0u32, 10u32);
        let second = CoverageRow::new(5u32, 10u32);
        first.shift_to(5).unwrap();
        assert_eq!(first, second);

        // Test shift with overflow
        let mut first = CoverageRow::new(0u32, 10u32);
        let (start, end) = first.shift_to(12).unwrap();
        assert_eq!(start, 1u32);
        assert_eq!(end, 10u32);
    }

    #[test]
    fn coverage() {
        let mut coverage: Coverage<u64> = Coverage::new(vec![
            Sequence {
                name: "chr1".to_string(),
                len:  100,
            },
            Sequence {
                name: "chr2".to_string(),
                len:  200,
            },
        ]);
        assert_eq!(coverage.get("chr1").map(|v| v.read), Some(0u64));
        assert_eq!(coverage.get("chr2").map(|v| v.read), Some(0u64));

        let (start, end) = coverage.shift_to("chr1", 50).unwrap();
        assert_eq!(start, 1u64);
        assert_eq!(end, 50u64);

        let (start, end) = coverage.shift_to("chr1", 110).unwrap();
        assert_eq!(start, 51u64);
        assert_eq!(end, 100u64);
        assert!(coverage
            .get("chr1")
            .unwrap()
            .is_finished());

        let (start, end) = coverage.shift_to("chr2", 100).unwrap();
        assert_eq!(start, 1u64);
        assert_eq!(end, 100u64);
    }

    #[test]
    #[should_panic]
    fn coverage_without_finish() {
        let mut coverage: Coverage<u64> = Coverage::new(vec![
            Sequence {
                name: "chr1".to_string(),
                len:  100,
            },
            Sequence {
                name: "chr2".to_string(),
                len:  200,
            },
        ]);
        let _result = coverage.shift_to("chr2", 100).unwrap();
    }

    #[test]
    #[should_panic]
    fn coverage_unknown_chr() {
        let mut coverage: Coverage<u64> = Coverage::new(vec![
            Sequence {
                name: "chr1".to_string(),
                len:  100,
            },
            Sequence {
                name: "chr2".to_string(),
                len:  200,
            },
        ]);
        let _result = coverage.shift_to("chr3", 100).unwrap();
    }
}
