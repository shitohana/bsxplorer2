use pyo3::prelude::*;

mod bsx;
mod compression;
mod region;
mod report;

pub fn register_io_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let module = PyModule::new_bound(parent_module.py(), "io")?;
    module.add_class::<bsx::PyBsxFileReader>()?;
    module.add_class::<bsx::PyBsxFileWriter>()?;
    module.add_class::<bsx::PyIpcCompression>()?;
    module.add_class::<compression::PyCompression>()?;
    module.add_class::<report::PyReportReader>()?;
    module.add_class::<report::PyReportWriter>()?;
    module.add_class::<region::PyRegionReader>()?;

    parent_module.add_submodule(&module)?;
    Ok(())
}
