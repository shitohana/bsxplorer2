// use crate::reader::ReportType;
// use bsx_rs::io::bsx::write::ConvertReportOptions;
// use pyo3::prelude::*;
// use pyo3_polars::error::PyPolarsErr;
//
// #[pyfunction]
// pub fn convert_report(
//     report_type: ReportType,
//     report_path: &str,
//     output_path: &str,
//     fa_path: &str,
//     fai_path: &str,
//     batch_per_read: Option<usize>,
//     batch_size: Option<usize>,
//     chunk_size: Option<usize>,
//     low_memory: Option<bool>,
//     n_threads: Option<usize>,
//     rechunk: Option<bool>,
// ) -> PyResult<()> {
//     let report_type = report_type.into_rust();
//     let convert_options = ConvertReportOptions::default()
//         .with_batch_size(batch_size.unwrap_or_default())
//         .with_low_memory(low_memory.unwrap_or_default())
//         .with_n_threads(n_threads.unwrap_or_default())
//         .with_rechunk(rechunk.unwrap_or_default())
//         .with_batch_per_read(batch_per_read.unwrap_or_default())
//         .with_chunk_size(chunk_size.unwrap_or_default());
//     match convert_options.finish(report_type, report_path, output_path, fa_path, fai_path) {
//         Ok(_) => Ok(()),
//         Err(e) => Err(PyErr::from(PyPolarsErr::from(e))),
//     }
// }
