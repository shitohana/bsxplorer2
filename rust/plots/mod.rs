#[cfg(feature = "plots")]
mod lineplot;
#[cfg(feature = "plots")]
pub use lineplot::{add_line, LinePlot, LinePlotData};

#[cfg(feature = "plots")]
mod inner {
    use crate::{
        data_structs::{bsx_batch::EncodedBsxBatch, region_data::RegionData},
        utils::types::{PosNum, RefId},
    };

    /// Basic plot trait
    pub trait BsxPlot {
        fn add_region_data<R: RefId, N: PosNum>(
            &mut self,
            region_data: &RegionData<R, N, EncodedBsxBatch>,
        ) -> anyhow::Result<()>;
    }

    /// Modular plot trait
    pub trait ModularBsxPlot<P: BsxPlot> {
        fn modules(&self) -> Vec<P>;
        fn modules_mut(&mut self) -> Vec<&mut P>;

        fn add_module(&mut self, mut module: P) {
            self.modules_mut().push(&mut module);
        }
    }
}

pub(crate) use inner::*;
