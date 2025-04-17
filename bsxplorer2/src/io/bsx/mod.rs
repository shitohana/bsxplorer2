mod ipc;
mod read;
mod region_read;
mod write;

pub use {
    read::{BsxFileReader, BSXIndex},
    region_read::{RegionReader},
    write::BsxIpcWriter,
};