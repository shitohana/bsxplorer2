use crate::data_structs::region::{GenomicPosition, RegionCoordinates};
use crate::utils::types::Strand;
use bio::io::fasta::*;
use hashbrown::HashMap;
use log::{debug, info};
use num::{PrimInt, ToPrimitive, Unsigned};
use serde::Serialize;
use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::io::{BufRead, Read, Seek};

pub struct FastaReader<R>
where
    R: Read + BufRead + Seek,
{
    reader: IndexedReader<R>,
}

impl<R> FastaReader<R>
where
    R: Read + BufRead + Seek,
{
    pub fn new(reader: IndexedReader<R>) -> FastaReader<R> {
        info!("Fasta reader initialized");
        Self { reader }
    }

    pub fn try_from_handle(fasta_handle: R, fai_handle: R) -> Result<Self, Box<dyn Error>> {
        let reader = IndexedReader::new(fasta_handle, fai_handle)?;
        Ok(Self::new(reader))
    }

    pub fn index(&self) -> &Index {
        &self.reader.index
    }

    pub fn fetch_region<N>(
        &mut self,
        region_coordinates: RegionCoordinates<N>,
    ) -> Result<Vec<u8>, Box<dyn Error>>
    where
        N: PrimInt + Unsigned + Display + Serialize,
    {
        if let Some(seq_data) = self.get_record_metadata(region_coordinates.chr()) {
            if seq_data.len < region_coordinates.end().to_u64().unwrap() {
                Err(CoverageError::OutOfBounds(
                    region_coordinates.chr().into(),
                    region_coordinates.end().to_usize().unwrap(),
                )
                .into())
            } else {
                self.reader.fetch(
                    region_coordinates.chr(),
                    region_coordinates.start().to_u64().unwrap(),
                    if region_coordinates.end != N::from(seq_data.len).unwrap() {
                        region_coordinates.end().to_u64().unwrap() + 1
                    } else {
                        seq_data.len
                    },
                )?;
                let mut buffer: Vec<u8> =
                    Vec::with_capacity(region_coordinates.length().to_usize().unwrap());
                self.reader.read(&mut buffer)?;

                debug!("Fetched sequence for {}", region_coordinates);
                Ok(buffer)
            }
        } else {
            Err(CoverageError::ChromosomeNotFound(region_coordinates.chr().into()).into())
        }
    }

    pub fn fetch_record(&mut self, name: &str) -> Result<Vec<u8>, Box<dyn Error>> {
        if let Some(seq_data) = self.get_record_metadata(name) {
            let mut buffer: Vec<u8> = Vec::with_capacity(seq_data.len as usize);
            self.reader.fetch_all(name)?;
            self.reader.read(&mut buffer)?;
            Ok(buffer)
        } else {
            Err(CoverageError::ChromosomeNotFound(name.into()).into())
        }
    }

    pub fn get_record_metadata(&self, name: &str) -> Option<Sequence> {
        self.reader
            .index
            .sequences()
            .into_iter()
            .find(|s| s.name == name)
    }

    pub fn get_record_names(&self) -> Vec<String> {
        self.reader
            .index
            .sequences()
            .into_iter()
            .map(|s| s.name.clone())
            .collect()
    }
}

pub struct FastaCoverageReader<R, V>
where
    R: Read + BufRead + Seek,
    V: PrimInt + Unsigned,
{
    inner: FastaReader<R>,
    coverage: Coverage<V>,
}

impl<R, V> From<FastaReader<R>> for FastaCoverageReader<R, V>
where
    R: Read + BufRead + Seek,
    V: PrimInt + Unsigned + ToPrimitive,
{
    fn from(fasta_reader: FastaReader<R>) -> Self {
        let coverage: Coverage<V> = Coverage::new(fasta_reader.reader.index.sequences());
        Self {
            inner: fasta_reader,
            coverage,
        }
    }
}

impl<R, V> FastaCoverageReader<R, V>
where
    R: Read + BufRead + Seek,
    V: PrimInt + Unsigned + ToPrimitive,
{
    pub fn new(fasta_reader: FastaReader<R>) -> Self {
        Self::from(fasta_reader)
    }

    pub fn tell_pos(&self, name: &str) -> Option<V> {
        self.coverage.get(name.to_string()).map(|c| c.read)
    }

    pub fn inner(&self) -> &FastaReader<R> {
        &self.inner
    }

    pub(crate) fn inner_mut(&mut self) -> &mut FastaReader<R> {
        &mut self.inner
    }

    pub fn coverage(&self) -> &Coverage<V> {
        &self.coverage
    }

    pub(crate) fn coverage_mut(&mut self) -> &mut Coverage<V> {
        &mut self.coverage
    }

    pub fn get_seq_until<N: PrimInt + Unsigned + Display + Serialize>(
        &mut self,
        position: &GenomicPosition<N>,
    ) -> Result<Vec<u8>, Box<dyn Error>> {
        let (chr, pos) = (position.chr(), position.position());

        let (fetch_start, fetch_end) = self.coverage.shift_to(chr, V::from(pos).unwrap())?;

        let read_region = RegionCoordinates::new(
            chr.to_string(),
            fetch_start.to_u32().unwrap(),
            fetch_end.to_u32().unwrap(),
            Strand::None,
        );
        let sequence = self.inner.fetch_region(read_region)?;
        Ok(sequence)
    }

    pub fn get_seq_leftover(&mut self, chr: &str) -> Result<Vec<u8>, Box<dyn Error>> {
        let end_gpos = GenomicPosition::new(
            chr.to_string(),
            self.coverage
                .get(chr.to_string())
                .unwrap()
                .total
                .to_u32()
                .unwrap(),
        );
        self.get_seq_until(&end_gpos)
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct CoverageRow<V>
where
    V: Unsigned + PrimInt,
{
    read: V,
    total: V,
}

impl<V> CoverageRow<V>
where
    V: Unsigned + PrimInt,
{
    pub fn new(read: V, total: V) -> Self {
        Self { read, total }
    }
    pub fn empty(total: V) -> Self {
        Self {
            read: V::from(0).unwrap(),
            total,
        }
    }
    pub fn shift_to(&mut self, pos: V) -> Result<(V, V), CoverageError> {
        if self.read > pos {
            return Err(CoverageError::ValueError(format!(
                "New pos {} is less than previous {}",
                pos.to_u64().unwrap(),
                self.read.to_u64().unwrap()
            )));
        }

        let prev = self.read + V::from(1).unwrap();
        self.read = if pos <= self.total { pos } else { self.total };
        Ok((prev, self.read))
    }
    pub fn is_finished(&self) -> bool {
        self.total == self.read
    }

    pub fn read(&self) -> V {
        self.read
    }
    pub fn total(&self) -> V {
        self.total
    }
}

pub struct Coverage<V>
where
    V: Unsigned + PrimInt,
{
    mapping: HashMap<String, CoverageRow<V>>,
    order: Vec<String>,
}

pub enum CoverageError {
    ChromosomeNotFound(String),
    PreviousNotFinished(String),
    OutOfBounds(String, usize),
    ValueError(String),
}

impl Debug for CoverageError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            CoverageError::ChromosomeNotFound(name) => write!(f, "Chromosome '{}' not found", name),
            CoverageError::PreviousNotFinished(name) => {
                write!(f, "Previous sequence '{}' not finished", name)
            }
            CoverageError::OutOfBounds(name, size) => {
                write!(f, "Position {} for '{}' is out of bounds", name, size)
            }
            CoverageError::ValueError(descr) => {
                write!(f, "{}", descr)
            }
        }
    }
}

impl Display for CoverageError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(&self, f)
    }
}

impl Error for CoverageError {}

impl<V> Coverage<V>
where
    V: Unsigned + PrimInt,
{
    fn new(sequences: Vec<Sequence>) -> Self {
        let mapping = HashMap::from_iter(sequences.iter().map(|seq| {
            (
                seq.name.to_string(),
                CoverageRow::empty(V::from(seq.len).unwrap()),
            )
        }));
        let order = sequences
            .iter()
            .cloned()
            .map(|seq| seq.name.to_string())
            .collect();
        Self { mapping, order }
    }

    pub fn get(&self, name: String) -> Option<&CoverageRow<V>> {
        self.mapping.get(&name)
    }

    pub fn position(&self, name: String) -> Option<usize> {
        self.order.iter().position(|n| *n == name)
    }

    pub fn shift_to(&mut self, name: &str, pos: V) -> Result<(V, V), CoverageError> {
        if let Some(idx) = self.position(name.to_string()) {
            if idx > 0 {
                if let Some(name) = self.order[..idx]
                    .iter()
                    .find(|n| !self.mapping.get(*n).unwrap().is_finished())
                {
                    return Err(CoverageError::PreviousNotFinished(name.to_string()));
                }
            }
        }
        if let Some(coverage) = self.mapping.get_mut(name) {
            coverage.shift_to(pos)
        } else {
            Err(CoverageError::ChromosomeNotFound(name.to_owned()))
        }
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
                len: 100,
            },
            Sequence {
                name: "chr2".to_string(),
                len: 200,
            },
        ]);
        assert_eq!(coverage.get("chr1".to_string()).map(|v| v.read), Some(0u64));
        assert_eq!(coverage.get("chr2".to_string()).map(|v| v.read), Some(0u64));

        let (start, end) = coverage.shift_to("chr1", 50).unwrap();
        assert_eq!(start, 1u64);
        assert_eq!(end, 50u64);
        let (start, end) = coverage.shift_to("chr1", 110).unwrap();
        assert_eq!(start, 51u64);
        assert_eq!(end, 100u64);
        assert!(coverage.get("chr1".to_string()).unwrap().is_finished());
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
                len: 100,
            },
            Sequence {
                name: "chr2".to_string(),
                len: 200,
            },
        ]);
        let (_start, _end) = coverage.shift_to("chr2", 100).unwrap();
    }

    #[test]
    #[should_panic]
    fn coverage_unknown_chr() {
        let mut coverage: Coverage<u64> = Coverage::new(vec![
            Sequence {
                name: "chr1".to_string(),
                len: 100,
            },
            Sequence {
                name: "chr2".to_string(),
                len: 200,
            },
        ]);
        let (start, end) = coverage.shift_to("chr3", 100).unwrap();
    }
}
