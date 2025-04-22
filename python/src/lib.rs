#![allow(unsafe_op_in_unsafe_fn, unused)]
#![warn(unused_imports, unused_braces)]
use pyo3::prelude::*;
mod types;
use types::register_data_structs_module;


#[pymodule]
fn bsx2(m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_data_structs_module(m)?;
    Ok(())
}