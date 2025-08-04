#![allow(unsafe_op_in_unsafe_fn, unused)]
#![warn(unused_imports, unused_braces)]

mod io;
mod data_structs;
mod utils;

use pyo3::prelude::*;

#[pymodule]
fn _bsx2(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<data_structs::utils::PyStrand>()?;
    m.add_class::<data_structs::utils::PyContext>()?;

    m.add_class::<data_structs::batch::PyBsxColumns>()?;
    m.add_class::<data_structs::batch::PyBsxBatch>()?;
    m.add_class::<data_structs::batch::PyAggMethod>()?;

    m.add_class::<data_structs::context_data::PyContextData>()?;

    m.add_class::<data_structs::report_schema::PyReportTypeSchema>()?;

    m.add_class::<data_structs::lazy::PyLazyBsxBatch>()?;

    m.add_class::<data_structs::coords::PyContig>()?;
    m.add_class::<data_structs::coords::PyGenomicPosition>()?;

    m.add_class::<data_structs::annot::PyAnnotStore>()?;
    m.add_class::<data_structs::annot::PyGffEntry>()?;
    m.add_class::<data_structs::annot::PyAnnotStoreIterator>()?;
    m.add_class::<data_structs::annot::PyGffEntryAttributes>()?;
    m.add_class::<data_structs::index::PyBatchIndex>()?;

    m.add_class::<io::bsx::PyBsxFileReader>()?;
    m.add_class::<io::bsx::PyBsxFileWriter>()?;
    m.add_class::<io::bsx::PyIpcCompression>()?;
    m.add_class::<io::compression::PyCompression>()?;
    m.add_class::<io::report::PyReportReader>()?;
    m.add_class::<io::report::PyReportWriter>()?;
    m.add_class::<io::region::PyRegionReader>()?;
    m.add_class::<io::region::PyFilterOperation>()?;
    m.add_class::<io::region::PyRegionReaderIterator>()?;

    Ok(())
}
