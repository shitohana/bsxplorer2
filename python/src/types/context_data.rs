use bsxplorer2::data_structs::batch::{BsxBatch as RsBsxBatch,
                                      EncodedBsxBatch as RsEncodedBsxBatch};
use bsxplorer2::data_structs::context_data::ContextData as RustContextData;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

#[pyclass(name = "ContextData")]
#[derive(Clone)]
pub struct PyContextData {
    data: RustContextData,
}

impl From<RustContextData> for PyContextData {
    fn from(data: RustContextData) -> Self { Self { data } }
}

impl From<PyContextData> for RustContextData {
    fn from(py: PyContextData) -> Self { py.data }
}

#[pymethods]
impl PyContextData {
    #[new]
    fn new(sequence: Vec<u8>) -> Self {
        Self {
            data: RustContextData::from_sequence(&sequence),
        }
    }

    fn to_decoded_df(&self) -> PyResult<PyDataFrame> {
        let df = self.data.clone().to_df::<RsBsxBatch>();
        Ok(PyDataFrame(df))
    }

    pub fn to_encoded_df(&self) -> PyResult<PyDataFrame> {
        let df = self
            .data
            .clone()
            .to_df::<RsEncodedBsxBatch>();
        Ok(PyDataFrame(df))
    }
}
