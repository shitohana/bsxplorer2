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
#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg(feature = "python")]
#[pymodule]
fn _lib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // ========== Top Module
    use crate::region::{PyGenomicPosition, PyRegionCoordinates};
    m.add_class::<PyGenomicPosition>()?;
    m.add_class::<PyRegionCoordinates>()?;

    use crate::bsx_batch::{BsxBatch, EncodedBsxBatch};
    m.add_class::<BsxBatch>()?;
    m.add_class::<EncodedBsxBatch>()?;

    // ========== Types
    fn register_types_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
        use crate::utils::types::{Context, Strand};
        let types_module = PyModule::new_bound(parent_module.py(), "types")?;
        types_module.add_class::<Context>()?;
        types_module.add_class::<Strand>()?;

        parent_module.add_submodule(&types_module)
    }
    register_types_module(m)?;

    // ========== IO
    fn register_io_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
        // ========== Report
        fn register_report_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
            use crate::io::report::*;
            let child_module = PyModule::new_bound(parent_module.py(), "report")?;

            child_module.add_class::<schema::ReportTypeSchema>()?;
            child_module.add_class::<read::ReportReader>()?;

            parent_module.add_submodule(&child_module)?;
            Ok(())
        }

        // ========== BSX
        fn register_bsx_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
            use crate::io::report::*;
            let child_module = PyModule::new_bound(parent_module.py(), "bsx")?;

            child_module.add_class::<bsx::region_read::RegionReader>()?;
            child_module.add_class::<bsx::region_read::PyBsxFileReader>()?;
            child_module.add_class::<bsx::write::PyIpcCompression>()?;
            child_module.add_class::<bsx::write::PyBsxIpcWriter>()?;

            parent_module.add_submodule(&child_module)?;
            Ok(())
        }

        let io_module = PyModule::new_bound(parent_module.py(), "io")?;
        register_report_module(&io_module)?;
        register_bsx_module(&io_module)?;

        Ok(())
    }
    register_io_module(m)?;

    Ok(())
}
