#![allow(unsafe_op_in_unsafe_fn, unused)]
#![warn(unused_imports, unused_braces)]

mod io;
mod types;
mod utils;

use io::register_io_module;
use pyo3::prelude::*;
use types::register_data_structs_module;

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_data_structs_module(m)?;
    register_io_module(m)?;
    Ok(())
}
