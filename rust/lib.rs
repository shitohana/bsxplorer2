#![feature(btree_cursors)]
#![feature(ascii_char)]
#![feature(vec_pop_if)]
#![feature(mpmc_channel)]
#![feature(assert_matches)]

pub mod bsx_batch;
pub mod genome;
pub mod io;
pub mod region;
pub mod utils;

use crate::io::bsx;
use crate::io::report::schema;
#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg(feature = "python")]
#[pymodule]
fn _lib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Region
    use crate::region::{PyGenomicPosition, PyRegionCoordinates};
    m.add_class::<PyGenomicPosition>()?;
    m.add_class::<PyRegionCoordinates>()?;

    // BsxBatch
    use crate::bsx_batch::{BsxBatch, EncodedBsxBatch};
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
