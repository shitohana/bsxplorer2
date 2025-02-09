use crate::data_structs::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use crate::data_structs::region::RegionCoordinates;
use crate::utils::types::{Context, Data, PosNum, RefId, Strand};
use itertools::Itertools;
use polars::prelude::PolarsResult;
use rayon::prelude::*;
use std::alloc::GlobalAlloc;
use std::collections::HashMap;
use std::fmt::Display;
use std::hash::Hash;

struct RegionData<R, N, T>
where
    R: RefId,
    N: PosNum,
    T: Data,
{
    chr: R,
    start: N,
    end: N,
    strand: Strand,
    data: T,
    attributes: HashMap<String, String>,
}

impl<R, N> From<RegionCoordinates<N>> for RegionData<R, N, ()>
where
    R: RefId,
    N: PosNum,
{
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

    pub fn length(&self) -> N {
        self.end.clone() - self.start.clone()
    }

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
            attributes: HashMap::default(),
        })
    }

    pub fn into_data(self) -> T {
        self.data
    }
}

impl<R, N, T> RegionData<R, N, T>
where
    R: RefId,
    N: PosNum,
    T: Data,
    String: From<R>,
{
    pub fn into_region(self) -> RegionCoordinates<N> {
        RegionCoordinates {
            chr: self.chr.into(),
            start: self.start,
            end: self.end,
            strand: self.strand,
        }
    }
}

impl<R, N> RegionData<R, N, EncodedBsxBatch>
where
    R: RefId,
    N: PosNum,
{
    pub fn meth_density(&self) -> PolarsResult<Vec<f32>> {
        self.data.get_density_vals()
    }
    pub fn get_positions(&self) -> PolarsResult<Vec<u32>> {
        self.data.get_position_vals()
    }
    pub fn size(&self) -> usize {
        self.data.height()
    }
    pub fn discretize(&self, resolution: usize) -> PolarsResult<Vec<usize>> {
        let length = self.length();
        let positions_norm = self
            .get_positions()?
            .into_par_iter()
            .map(|x| (x as f64) / length.to_f64().unwrap() * resolution as f64)
            .collect::<Vec<f64>>();
        let mut resolution_iter = (0..resolution).into_iter().peekable();
        let sort_indices = positions_norm
            .into_iter()
            .enumerate()
            .filter_map(|(idx, x)| match resolution_iter.peek() {
                Some(&y) if x > y as f64 => {
                    resolution_iter.next().unwrap();
                    Some(idx)
                }
                _ => None,
            })
            .collect_vec();
        Ok(sort_indices)
    }
    pub fn discrete_density(&self, resolution: usize) -> PolarsResult<Vec<f32>> {
        let sort_indices = self.discretize(resolution)?;
        let meth_density = self.meth_density()?;
        let mut density = vec![0.0; resolution];
        for (idx, val) in sort_indices.into_iter().zip(meth_density.into_iter()) {
            density[idx] = val;
        }
        Ok(density)
    }
    pub fn discrete_counts(&self, resolution: usize) -> PolarsResult<Vec<(i16, i16)>> {
        let sort_indices = self.discretize(resolution)?;
        let counts_meth = self.counts_meth()?;
        let counts_total = self.counts_total()?;
        let mut counts = vec![(0, 0); resolution];
        for (idx, (meth, total)) in sort_indices
            .into_iter()
            .zip(counts_meth.into_iter().zip(counts_total.into_iter()))
        {
            counts[idx] = (meth, total);
        }
        Ok(counts)
    }
    pub fn counts_meth(&self) -> PolarsResult<Vec<i16>> {
        self.data.get_counts_m()
    }
    pub fn counts_total(&self) -> PolarsResult<Vec<i16>> {
        self.data.get_counts_total()
    }
    pub fn filter(mut self, context: Context) -> Self {
        self.data = self.data.filter(Some(context), None);
        self
    }
}
