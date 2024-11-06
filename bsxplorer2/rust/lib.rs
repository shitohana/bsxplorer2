#![feature(btree_cursors)]
#![feature(ascii_char)]

mod python;
use pyo3::prelude::*;
use python::python::{
    ReportReader, ReportType
};

#[allow(dead_code)]
#[pymodule]
fn _lib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ReportReader>()?;
    m.add_class::<ReportType>()?;
    Ok(())
}
