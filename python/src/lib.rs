use pyo3::prelude::*;
mod data_structs;


#[pymodule]
fn bsx2(m: &Bound<'_, PyModule>) -> PyResult<()> {
    Ok(())
}