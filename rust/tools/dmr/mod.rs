pub mod config;
mod data_structs;
mod iterator;
mod segmentation;
pub(crate) mod tv1d_clone;

pub use iterator::*;
pub use data_structs::{DMRegion};