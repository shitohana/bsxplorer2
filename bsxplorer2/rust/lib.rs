#![feature(btree_cursors)]
#![feature(ascii_char)]
#![feature(variant_count)]

mod reader;
mod writer;
mod bsx_reader;
mod region;
mod utils;

use pyo3::prelude::*;
use reader::{ReportReader, ReportType};

#[allow(dead_code)]
#[pymodule]
fn _lib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ReportReader>()?;
    m.add_class::<ReportType>()?;
    m.add_function(wrap_pyfunction!(writer::convert_report, m)?);
    Ok(())
}
