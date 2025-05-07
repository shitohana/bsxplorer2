use pyo3::prelude::*;

pub mod annot;
pub mod batch;
pub mod context_data;
pub mod coords;
pub mod index;
pub mod lazy;
pub mod report_schema;
pub mod stats;
pub mod utils;

pub fn register_data_structs_module(
    parent_module: &Bound<'_, PyModule>
) -> PyResult<()> {
    let module = PyModule::new_bound(parent_module.py(), "types")?;

    module.add_class::<utils::PyStrand>()?;
    module.add_class::<utils::PyContext>()?;

    module.add_class::<batch::PyBsxBatch>()?;
    module.add_class::<batch::PyEncodedBsxBatch>()?;

    module.add_class::<context_data::PyContextData>()?;

    module.add_class::<report_schema::PyReportTypeSchema>()?;

    module.add_class::<lazy::PyLazyBsxBatch>()?;
    module.add_class::<lazy::PyLazyEncodedBsxBatch>()?;

    module.add_class::<coords::PyContig>()?;
    module.add_class::<coords::PyGenomicPosition>()?;

    module.add_class::<stats::PyMethylationStats>()?;

    module.add_class::<annot::PyAnnotStore>()?;
    module.add_class::<index::PyBatchIndex>()?;

    parent_module.add_submodule(&module)?;
    Ok(())
}
