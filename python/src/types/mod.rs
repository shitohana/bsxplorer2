use pyo3::prelude::*;

pub mod context_data;
pub mod coords;
pub mod decoded;
pub mod encoded;
pub mod report_schema;
pub mod utils;

pub fn register_data_structs_module(
    parent_module: &Bound<'_, PyModule>
) -> PyResult<()> {
    let module = PyModule::new_bound(parent_module.py(), "types")?;
    module.add_class::<utils::Strand>()?;
    module.add_class::<utils::Context>()?;
    module.add_class::<decoded::PyBsxBatch>()?;
    module.add_class::<encoded::PyEncodedBsxBatch>()?;
    module.add_class::<context_data::PyContextData>()?;
    module.add_class::<report_schema::PyReportTypeSchema>()?;

    parent_module.add_submodule(&module)?;
    Ok(())
}
