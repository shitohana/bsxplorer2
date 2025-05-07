use bsxplorer2::data_structs::batch::{BsxBatch,
                                      BsxBatchMethods,
                                      EncodedBsxBatch,
                                      LazyBsxBatch};
use paste::paste;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

use super::batch::{PyBsxBatch, PyEncodedBsxBatch};
use super::context_data::PyContextData;
use super::report_schema::PyReportTypeSchema;
use super::stats::PyMethylationStats;
use super::utils::{PyContext, PyStrand};

macro_rules! lazy_bsx_batch_wrapper {
    ($wrap_type: ident, $py_name: expr) => {
        paste! {
            #[pyclass(name = $py_name)]
            pub struct [<PyLazy $wrap_type>] {
                inner: LazyBsxBatch<$wrap_type>
            }

            impl From<[<PyLazy $wrap_type>]> for LazyBsxBatch<$wrap_type> {
                fn from(value: [<PyLazy $wrap_type>]) -> Self {
                    value.inner
                }
            }

            impl From<LazyBsxBatch<$wrap_type>> for [<PyLazy $wrap_type>] {
                fn from(value: LazyBsxBatch<$wrap_type>) -> Self {
                    Self { inner: value }
                }
            }

            #[pymethods]
            impl [<PyLazy $wrap_type>] {
                #[new]
                pub fn new(batch: [<Py $wrap_type>]) -> Self {
                    Self { inner: <$wrap_type>::from(batch).lazy() }
                }
                pub fn into_report(&self, report_type: PyReportTypeSchema) -> PyResult<PyDataFrame> {
                    self.inner.clone().into_report(&report_type.into())
                        .map(|df| PyDataFrame(df))
                        .map_err(PyErr::from)
                }
                pub fn collect(&self) -> PyResult<[<Py $wrap_type>]> {
                    self.inner.clone().collect().map(|val: $wrap_type| -> [<Py $wrap_type>] { val.into() }).map_err(PyErr::from)
                }
                pub fn filter_pos_lt(&self, pos: u32) -> Self {
                    Self { inner: self.inner.clone().filter_pos_lt(pos) }
                }
                pub fn filter_pos_gt(&self, pos: u32) -> Self {
                    Self { inner: self.inner.clone().filter_pos_gt(pos) }
                }
                pub fn filter_coverage_lt(&self, coverage: u32) -> Self {
                    Self { inner: self.inner.clone().filter_coverage_lt(coverage) }
                }
                pub fn filter_strand(&self, strand: PyStrand) -> Self {
                    Self { inner: self.inner.clone().filter_strand(strand.into()) }
                }
                pub fn filter_context(&self, context: PyContext) -> Self {
                    Self { inner: self.inner.clone().filter_context(context.into()) }
                }
                pub fn mark_low_coverage(&self, threshold: u32) -> Self {
                    Self { inner: self.inner.clone().mark_low_coverage(threshold) }
                }
                pub fn align_with_contexts(&self, context_data: PyContextData, chr_val: &str) -> Self {
                    Self { inner: self.inner.clone().align_with_contexts(context_data.into(), chr_val) }
                }
                pub fn stats(&self) -> PyResult<PyMethylationStats> {
                    Ok(self.inner.stats()?.into())
                }
            }
        }
    };
}

lazy_bsx_batch_wrapper!(BsxBatch, "LazyBsxBatch");
lazy_bsx_batch_wrapper!(EncodedBsxBatch, "LazyEncodedBsxBatch");
