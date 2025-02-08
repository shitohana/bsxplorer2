#![feature(btree_cursors)]
#![feature(ascii_char)]
#![feature(vec_pop_if)]
#![feature(mpmc_channel)]
#![feature(assert_matches)]
#![feature(path_add_extension)]
#![allow(unused_assignments)]

pub mod data_structs;
pub mod io;
pub mod tools;
pub mod utils;

use crate::io::bsx;
use crate::io::report::schema;
#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg(feature = "python")]
#[pymodule]
fn _lib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Region
    use data_structs::region::{PyGenomicPosition, PyRegionCoordinates};
    m.add_class::<PyGenomicPosition>()?;
    m.add_class::<PyRegionCoordinates>()?;

    // BsxBatch
    use data_structs::bsx_batch::{BsxBatch, EncodedBsxBatch};
    m.add_class::<BsxBatch>()?;
    m.add_class::<EncodedBsxBatch>()?;

    // Types
    use crate::utils::types::{Context, Strand};
    m.add_class::<Context>()?;
    m.add_class::<Strand>()?;

    // IO
    use crate::io::report::*;
    // --- Report
    m.add_class::<schema::ReportTypeSchema>()?;
    m.add_class::<read::ReportReader>()?;
    m.add_class::<write::PyReportWriter>()?;
    // --- Bsx
    m.add_class::<bsx::region_read::RegionReader>()?;
    m.add_class::<bsx::region_read::PyBsxFileReader>()?;
    m.add_class::<bsx::write::PyIpcCompression>()?;
    m.add_class::<bsx::write::PyBsxIpcWriter>()?;

    Ok(())
}
