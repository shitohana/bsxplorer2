use crate::data_structs::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use crate::data_structs::region::RegionCoordinates;
use crate::utils::types::Strand;
use num::{PrimInt, Unsigned};
use serde::Serialize;
use std::fmt::Display;
use std::hash::Hash;

struct RegionData<R, N, T>
where
    R: Eq + Hash,
    N: PrimInt + Unsigned + Clone + Serialize + Display,
    T: Sized + Clone,
{
    chr: R,
    start: N,
    end: N,
    strand: Strand,
    data: T,
}

impl<R, N> From<RegionCoordinates<N>> for RegionData<R, N, ()>
where
    R: Eq + Hash + From<String>,
    N: PrimInt + Unsigned + Clone + Serialize + Display,
{
    fn from(region: RegionCoordinates<N>) -> Self {
        RegionData {
            chr: region.chr.into(),
            start: region.start,
            end: region.end,
            strand: region.strand,
            data: (),
        }
    }
}

impl<R, N, T> RegionData<R, N, T>
where
    R: Eq + Hash,
    N: PrimInt + Unsigned + Clone + Serialize + Display,
    T: Sized + Clone,
    String: From<R>,
{
    pub fn new(chr: R, start: N, end: N, strand: Strand, data: T) -> Self {
        RegionData {
            chr,
            start,
            end,
            strand,
            data,
        }
    }
    pub fn from_region_and_value(
        region: RegionCoordinates<N>,
        value: T,
    ) -> RegionData<String, N, T> {
        RegionData {
            chr: region.chr,
            start: region.start,
            end: region.end,
            strand: region.strand,
            data: value,
        }
    }

    pub fn try_from_batch(
        batch: &EncodedBsxBatch,
    ) -> anyhow::Result<RegionData<String, u64, &EncodedBsxBatch>> {
        let first_position = batch.first_position()?;

        Ok(RegionData {
            chr: first_position.chr().to_string(),
            start: first_position.position(),
            end: batch.last_position()?.position(),
            strand: Strand::None,
            data: batch,
        })
    }

    pub fn into_data(self) -> T {
        self.data
    }

    pub fn into_region(self) -> RegionCoordinates<N> {
        RegionCoordinates {
            chr: self.chr.into(),
            start: self.start,
            end: self.end,
            strand: self.strand,
        }
    }

    // Getters
    pub fn chr(&self) -> &R {
        &self.chr
    }

    pub fn start(&self) -> N {
        self.start
    }

    pub fn end(&self) -> N {
        self.end
    }

    pub fn strand(&self) -> Strand {
        self.strand
    }

    pub fn data(&self) -> &T {
        &self.data
    }

    // Setters
    pub fn set_chr(&mut self, chr: R) {
        self.chr = chr;
    }

    pub fn set_start(&mut self, start: N) {
        self.start = start;
    }

    pub fn set_end(&mut self, end: N) {
        self.end = end;
    }

    pub fn set_strand(&mut self, strand: Strand) {
        self.strand = strand;
    }

    pub fn set_data(&mut self, data: T) {
        self.data = data;
    }
}
