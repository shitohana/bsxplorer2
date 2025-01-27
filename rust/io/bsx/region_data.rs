use crate::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use crate::region::RegionCoordinates;
use itertools::Itertools;
use polars::prelude::*;
use std::fmt::Display;

pub struct RegionData<T: num::PrimInt + num::Unsigned> {
    batch: EncodedBsxBatch,
    coordinates: RegionCoordinates<T>,
}

#[derive(Default)]
pub enum RegionSamplingStat {
    Max,
    Min,
    Median,
    Mode,
    #[default]
    Mean,
    Std,
    Var,
}

impl<T: num::PrimInt + num::Unsigned + Display> RegionData<T> {
    fn density(&self) -> Vec<Option<f32>> {
        self.batch
            .data()
            .column("density")
            .unwrap()
            .f32()
            .unwrap()
            .into_iter()
            .collect_vec()
    }

    fn density_discrete(
        &self,
        sampling_rate: usize,
        stat: RegionSamplingStat,
    ) -> PolarsResult<Vec<Option<f32>>> {
        let pos_col = self.batch.data().column("pos")?.u32()?;
        let step = self.coordinates.length() / T::from(sampling_rate).unwrap();
        todo!()
    }
}
