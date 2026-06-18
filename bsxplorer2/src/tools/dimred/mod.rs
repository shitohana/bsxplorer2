#![cfg_attr(coverage_nightly, feature(coverage_attribute))]
mod dbscan;
mod merge;
mod segmentation;

pub use merge::{merge_breakpoints, EqFloat, MergeType};
pub use segmentation::{pelt, MethDataBinom, SegmentAlgorithm, SegmentationData};
