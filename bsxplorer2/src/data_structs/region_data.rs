use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use crate::data_structs::batch::BsxBatchMethods;
use crate::data_structs::batch::EncodedBsxBatch;
use crate::data_structs::region::RegionCoordinates;
use crate::utils::types::{Data, PosNum, RefId, Strand};

/// Represents data_structs associated with a genomic region.
#[derive(Clone, Eq, PartialEq)]
pub struct RegionData<R, N, T>
where
    R: RefId,
    N: PosNum,
    T: Data,
{
    /// The chromosome ID.
    chr: R,
    /// The start position of the region.
    start: N,
    /// The end position of the region.
    end: N,
    /// The strand of the region.
    strand: Strand,
    /// The associated data_structs.
    data: T,
    /// Optional attributes associated with the region.
    attributes: HashMap<String, String>,
}

impl<R: RefId, N: PosNum, T: Data> Hash for RegionData<R, N, T> {
    /// Hashes the region data_structs based on chromosome, start, end, and
    /// strand.
    fn hash<H: Hasher>(
        &self,
        state: &mut H,
    ) {
        (self.chr.clone(), self.start, self.end, self.strand).hash(state)
    }
}

impl<R, N> From<RegionCoordinates<N>> for RegionData<R, N, ()>
where
    R: RefId,
    N: PosNum,
{
    /// Creates a new `RegionData` from region coordinates with no data_structs.
    fn from(region: RegionCoordinates<N>) -> Self {
        RegionData {
            chr: region.chr.into(),
            start: region.start,
            end: region.end,
            strand: region.strand,
            data: (),
            attributes: HashMap::default(),
        }
    }
}

impl<R, N, T> RegionData<R, N, T>
where
    R: RefId,
    N: PosNum,
    T: Data,
{
    /// Returns a reference to the chromosome ID.
    pub fn chr(&self) -> &R {
        &self.chr
    }

    /// Returns the start position of the region.
    pub fn start(&self) -> N {
        self.start
    }

    /// Returns the end position of the region.
    pub fn end(&self) -> N {
        self.end
    }

    /// Returns the strand of the region.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns a reference to the associated data_structs.
    pub fn data(&self) -> &T {
        &self.data
    }

    /// Sets the chromosome ID.
    pub fn set_chr(
        &mut self,
        chr: R,
    ) {
        self.chr = chr;
    }

    /// Sets the start position.
    pub fn set_start(
        &mut self,
        start: N,
    ) {
        self.start = start;
    }

    /// Sets the end position.
    pub fn set_end(
        &mut self,
        end: N,
    ) {
        self.end = end;
    }

    /// Sets the strand.
    pub fn set_strand(
        &mut self,
        strand: Strand,
    ) {
        self.strand = strand;
    }

    /// Sets the associated data_structs.
    pub fn set_data(
        &mut self,
        data: T,
    ) {
        self.data = data;
    }

    /// Returns the length of the region.
    pub fn length(&self) -> N {
        self.end.clone() - self.start.clone()
    }

    /// Creates a new `RegionData` with the same region information but
    /// different data_structs.
    pub fn with_data<T2: Data>(
        self,
        data: T2,
    ) -> RegionData<R, N, T2> {
        RegionData {
            chr: self.chr,
            start: self.start,
            end: self.end,
            strand: self.strand,
            data,
            attributes: self.attributes,
        }
    }

    /// Creates a new `RegionData` instance.
    pub fn new(
        chr: R,
        start: N,
        end: N,
        strand: Strand,
        data: T,
        attributes: HashMap<String, String>,
    ) -> Self {
        RegionData {
            chr,
            start,
            end,
            strand,
            data,
            attributes,
        }
    }

    /// Creates a new `RegionData` from region coordinates and a value.
    pub fn from_region_and_value(
        region: RegionCoordinates<N>,
        value: T,
        attributes: HashMap<String, String>,
    ) -> RegionData<String, N, T> {
        RegionData {
            chr: region.chr,
            start: region.start,
            end: region.end,
            strand: region.strand,
            data: value,
            attributes,
        }
    }

    /// Tries to create a `RegionData` instance from an `EncodedBsxBatch`.
    pub fn try_from_batch(
        batch: &EncodedBsxBatch
    ) -> anyhow::Result<RegionData<String, u64, &EncodedBsxBatch>> {
        Ok(RegionData {
            chr: <EncodedBsxBatch as BsxBatchMethods>::chr_val(batch)?
                .to_string(),
            start: batch
                .start_pos()
                .ok_or(anyhow::anyhow!("empty data"))?
                as u64,
            end: batch
                .end_pos()
                .ok_or(anyhow::anyhow!("empty data"))? as u64,
            strand: Strand::None,
            data: batch,
            attributes: HashMap::default(),
        })
    }

    /// Consumes the `RegionData` and returns the associated data_structs.
    pub fn into_data(self) -> T {
        self.data
    }

    /// Separates the region data_structs from the associated data_structs,
    /// returning both.
    pub fn drain_data(self) -> (RegionData<R, N, ()>, T) {
        let data = self.data.clone();
        (self.with_data(()), data)
    }
}

impl<R, N, T> RegionData<R, N, T>
where
    R: RefId,
    N: PosNum,
    T: Data,
    String: From<R>,
{
    /// Converts the `RegionData` into a `RegionCoordinates` instance.
    pub fn into_region(self) -> RegionCoordinates<N> {
        RegionCoordinates {
            chr: self.chr.into(),
            start: self.start,
            end: self.end,
            strand: self.strand,
        }
    }

    /// Returns the region as a `RegionCoordinates` instance.
    pub fn as_region(&self) -> RegionCoordinates<N> {
        RegionCoordinates::new(
            self.chr.clone().into(),
            self.start,
            self.end,
            self.strand,
        )
    }
}
