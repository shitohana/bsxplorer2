#![cfg_attr(coverage_nightly, feature(coverage_attribute))]
mod dbscan;
mod segmentation;
mod merge;

pub use segmentation::{
    SegmentAlgorithm, SegmentationData, MethDataBinom, pelt
};
pub use merge::merge_breakpoints;
