pub mod config;
mod data_structs;
mod iterator;
mod segmentation;
#[cfg_attr(coverage_nightly, coverage(off))]
pub(crate) mod tv1d_clone;

pub use data_structs::DMRegion;
pub use iterator::*;
