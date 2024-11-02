use pyo3::prelude::*;
mod bsxipc;
mod genome;
mod report;
mod sequence;
mod utils;

pub use report::ReportType;

#[pymodule]
fn bsx_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ReportType>()?;
    Ok(())
}
