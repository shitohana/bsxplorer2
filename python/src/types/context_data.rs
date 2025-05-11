use bsxplorer2::data_structs::batch::{BsxBatch as RsBsxBatch,
                                      EncodedBsxBatch as RsEncodedBsxBatch};
use bsxplorer2::data_structs::context_data::ContextData as RustContextData;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

use super::utils::{PyContext, PyStrand};

#[pyclass(name = "ContextData")]
#[derive(Clone)]
pub struct PyContextData {
    inner: RustContextData,
}

impl From<RustContextData> for PyContextData {
    fn from(data: RustContextData) -> Self {
        Self { inner: data }
    }
}

impl From<PyContextData> for RustContextData {
    fn from(py: PyContextData) -> Self {
        py.inner
    }
}

#[pymethods]
impl PyContextData {
    #[new]
    fn new(sequence: Vec<u8>) -> Self {
        Self {
            inner: RustContextData::from_sequence(&sequence),
        }
    }

    #[staticmethod]
    fn empty() -> Self {
        Self {
            inner: RustContextData::empty(),
        }
    }

    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn take(&self) -> (Vec<u32>, Vec<PyStrand>, Vec<PyContext>) {
        let (ids, strands, contexts) = self.inner.clone().take();
        (
            ids,
            strands
                .into_iter()
                .map(PyStrand::from)
                .collect(),
            contexts
                .into_iter()
                .map(PyContext::from)
                .collect(),
        )
    }

    fn to_decoded_df(&self) -> PyResult<PyDataFrame> {
        let df = self.inner.clone().to_df::<RsBsxBatch>();
        Ok(PyDataFrame(df))
    }

    pub fn to_encoded_df(&self) -> PyResult<PyDataFrame> {
        let df = self
            .inner
            .clone()
            .to_df::<RsEncodedBsxBatch>();
        Ok(PyDataFrame(df))
    }

    pub fn __len__(&self) -> usize {
        self.inner.len()
    }
}
