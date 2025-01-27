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

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
#[pymodule]
fn _lib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    use crate::region::{PyGenomicPosition, PyRegionCoordinates};
    m.add_class::<PyGenomicPosition>()?;
    m.add_class::<PyRegionCoordinates>()?;

    Ok(())
}
