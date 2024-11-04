mod python;
use pyo3::prelude::*;
use python::python::{
    BSXFilePy as BSXFile,
    BSXIteratorPy as BSXIterator,
    ReportType
};

#[cfg_attr(feature = "python", pymodule)]
fn bsx_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<BSXFile>()?;
    m.add_class::<BSXIterator>()?;
    Ok(())
}
