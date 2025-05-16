mod index;
#[cfg_attr(coverage_nightly, coverage(off))]
mod ipc;
mod read;
mod region;
mod write;

pub use index::BatchIndex;
pub use read::{BsxFileIterator, BsxFileReader};
pub use region::RegionReader;
pub use write::BsxFileWriter;
