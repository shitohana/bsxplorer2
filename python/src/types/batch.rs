use std::sync::Arc;

use bsxplorer2::data_structs::batch::{BsxBatch,
                                      BsxBatchBuilder,
                                      BsxBatchMethods,
                                      EncodedBsxBatch};
use paste::paste;
use polars::prelude::{BooleanChunked,
                      DataFrame,
                      IntoSeries,
                      PolarsError,
                      Series};
use pyo3::prelude::*;
use pyo3_polars::{PyDataFrame, PySchema, PySeries};

use super::context_data::PyContextData;
use super::coords::{PyContig, PyGenomicPosition};
use super::report_schema::PyReportTypeSchema;

macro_rules! create_pyseries {
    ($data: expr) => {
        Ok(PySeries($data.clone().into_series()))
    };
}

macro_rules! lazy_bsx_batch_wrapper {
    ($wrap_type: ident, $py_name: expr) => {
        paste! {
            #[pyclass(name = $py_name)]
            #[derive(Debug, Clone)]
            pub struct [<Py $wrap_type>] {
                inner: $wrap_type
            }

            impl From<$wrap_type> for [<Py $wrap_type>] {
                fn from(inner: $wrap_type) -> Self { [<Py $wrap_type>] { inner } }
            }
            impl From<[<Py $wrap_type>]> for $wrap_type {
                fn from(py_batch: [<Py $wrap_type>]) -> Self { py_batch.inner }
            }

            #[pymethods]
            impl [<Py $wrap_type>] {
                #[new]
                #[pyo3(signature = (
                    data,
                    check_nulls=true,
                    check_sorted=true,
                    check_duplicates=true,
                    rechunk=true,
                    check_single_chr=true,
                    context_data=None,
                    report_schema=None
                ))]
                pub fn new(
                    data: PyDataFrame,
                    check_nulls: bool,
                    check_sorted: bool,
                    check_duplicates: bool,
                    rechunk: bool,
                    check_single_chr: bool,
                    context_data: Option<PyContextData>,
                    report_schema: Option<PyReportTypeSchema>,
                ) -> PyResult<Self> {
                    let df: DataFrame = data.into();
                    let builder = {
                        if let Some(report_schema) = report_schema {
                            BsxBatchBuilder::default()
                                .with_report_type(report_schema.to_rust())
                        }
                        else {
                            BsxBatchBuilder::default()
                        }
                    }
                    .with_check_nulls(check_nulls)
                    .with_check_sorted(check_sorted)
                    .with_check_duplicates(check_duplicates)
                    .with_rechunk(rechunk)
                    .with_check_single_chr(check_single_chr)
                    .with_context_data(context_data.map(|v| v.into()));

                    let inner = builder
                        .build::<$wrap_type>(df)
                        .map_err(|e| {
                            pyo3::exceptions::PyValueError::new_err(e.to_string())
                        })?;
                    Ok([<Py $wrap_type>] { inner })
                }

                #[staticmethod]
                #[allow(unsafe_code)]
                pub fn from_dataframe_unchecked(
                    data: PyDataFrame,
                ) -> PyResult<Self> {
                    let df: DataFrame = data.into();
                    let inner = unsafe { <$wrap_type>::new_unchecked(df) };
                    Ok([<Py $wrap_type>] { inner })
                }

                #[staticmethod]
                pub fn schema() -> PyResult<PySchema> {
                    Ok(PySchema(Arc::new(<$wrap_type>::schema())))
                }

                #[staticmethod]
                pub fn empty() -> PyResult<Self> {
                    Ok(Self {
                        inner: <$wrap_type>::empty(),
                    })
                }

                pub fn chr(&self) -> PyResult<PySeries> {
                    create_pyseries!(self.inner.chr())
                }

                pub fn strand(&self) -> PyResult<PySeries> {
                    create_pyseries!(self.inner.strand())
                }

                pub fn context(&self) -> PyResult<PySeries> {
                    create_pyseries!(self.inner.context())
                }

                pub fn count_m(&self) -> PyResult<PySeries> {
                    create_pyseries!(self.inner.count_m())
                }

                pub fn count_total(&self) -> PyResult<PySeries> {
                    create_pyseries!(self.inner.count_total())
                }

                pub fn density(&self) -> PyResult<PySeries> {
                    create_pyseries!(self.inner.density())
                }

                pub fn is_empty(&self) -> bool { self.inner.is_empty() }

                pub fn split_at(
                    &self,
                    index: usize,
                ) -> PyResult<(Self, Self)> {
                    let (batch1, batch2) = self.inner.clone().split_at(index);
                    Ok(([<Py $wrap_type>] { inner: batch1 }, [<Py $wrap_type>] { inner: batch2 }))
                }

                pub fn data(&self) -> PyResult<PyDataFrame> {
                    Ok(PyDataFrame(self.inner.data().clone()))
                }

                pub fn take(&self) -> PyResult<PyDataFrame> {
                    Ok(PyDataFrame(self.inner.clone().take()))
                }

                pub fn chr_val(&self) -> PyResult<String> {
                    Ok(self.inner.chr_val().to_string())
                }

                pub fn start_pos(&self) -> u32 { self.inner.start_pos() }
                pub fn end_pos(&self) -> u32 { self.inner.end_pos() }

                pub fn as_contig(&self) -> PyResult<Option<PyContig>> {
                    Ok(self.inner.as_contig()?.map(|val| val.into()))
                }

                pub fn start_gpos(&self) -> PyGenomicPosition {
                    self.inner.start_gpos().into()
                }
                pub fn end_gpos(&self) -> PyGenomicPosition {
                    self.inner.end_gpos().into()
                }

                pub fn vstack(
                    &self,
                    other: &[<Py $wrap_type>],
                ) -> PyResult<Self> {
                    let new_inner = self
                        .inner
                        .vstack(&other.inner)
                        .map_err(|e| {
                            pyo3::exceptions::PyValueError::new_err(e.to_string())
                        })?;
                    Ok([<Py $wrap_type>] { inner: new_inner })
                }

                pub fn extend(
                    &mut self,
                    other: &[<Py $wrap_type>],
                ) -> PyResult<()> {
                    self.inner
                        .extend(&other.inner)
                        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
                }

                pub fn filter_mask(
                    &self,
                    mask: PySeries,
                ) -> PyResult<Self> {
                    let mask_series: Series = mask.into();
                    let bool_mask: &BooleanChunked = mask_series.bool().map_err(|_| {
                        pyo3::exceptions::PyTypeError::new_err(
                            "Mask must be a boolean Series",
                        )
                    })?;

                    let new_inner = self
                        .inner
                        .filter_mask(bool_mask)
                        .map_err(|e: PolarsError| {
                            pyo3::exceptions::PyValueError::new_err(e.to_string())
                        })?;
                    Ok([<Py $wrap_type>] { inner: new_inner })
                }

                pub fn height(&self) -> usize { self.inner.height() }
                pub fn __len__(&self) -> usize { self.inner.height() }

                pub fn __repr__(&self) -> String {
                    format!(
                        "BsxBatch(shape={}, chr='{}')",
                        self.inner.data().shape().1,
                        self.inner.chr_val()
                    )
                }

                pub fn __richcmp__(
                    &self,
                    other: &Self,
                    op: pyo3::basic::CompareOp,
                ) -> PyResult<bool> {
                    match op {
                        pyo3::basic::CompareOp::Eq => Ok(self.inner == other.inner),
                        pyo3::basic::CompareOp::Ne => Ok(self.inner != other.inner),
                        _ => {
                            Err(pyo3::exceptions::PyNotImplementedError::new_err(
                                "Only == and != are supported for BsxBatch",
                            ))
                        },
                    }
                }
            }
        }
    };
}

lazy_bsx_batch_wrapper!(BsxBatch, "BsxBatch");
lazy_bsx_batch_wrapper!(EncodedBsxBatch, "EncodedBsxBatch");
