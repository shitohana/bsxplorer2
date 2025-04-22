use pyo3::prelude::*;
use bsxplorer2::data_structs::context_data::ContextData as RustContextData;
use bsxplorer2::data_structs::batch::{BsxBatch as RsBsxBatch, EncodedBsxBatch as RsEncodedBsxBatch};
use pyo3_polars::PyDataFrame;

/// Stores methylation context information derived from a DNA sequence.
///
/// This class pre-calculates the positions, strands, and contexts (CG, CHG, CHH)
/// for cytosines within a given sequence, allowing for efficient alignment
/// with methylation data.
#[pyclass(name="ContextData")]
#[derive(Clone)]
pub struct PyContextData {
    data: RustContextData,
}   

impl PyContextData {
    pub fn to_rust(self) -> RustContextData {
        self.data
    }
}

#[pymethods]
impl PyContextData {
    /// Creates a new ContextData instance from a DNA sequence.
    ///
    /// Parameters
    /// ----------
    /// sequence : bytes
    ///     The DNA sequence as bytes (e.g., b"ACGTACGT").
    ///
    /// Returns
    /// -------
    /// ContextData
    ///     A new instance containing context information.
    #[new]
    fn new(sequence: Vec<u8>) -> Self {
        Self {
            data: RustContextData::from_sequence(&sequence),
        }
    }

    /// Converts the context information into a Polars DataFrame compatible with BsxBatch (decoded format).
    ///
    /// The resulting DataFrame contains 'chr', 'position', 'strand', and 'context' columns
    /// suitable for joining or aligning with decoded methylation data.
    ///
    /// Returns
    /// -------
    /// DataFrame
    ///     A Polars DataFrame with context information.
    fn to_decoded_df(&self) -> PyResult<PyDataFrame> {
        let df = self.data.clone().to_df::<RsBsxBatch>();
        Ok(PyDataFrame(df))
    }

    /// Converts the context information into a Polars DataFrame compatible with EncodedBsxBatch.
    ///
    /// The resulting DataFrame contains 'chr', 'position', 'strand', and 'context' columns
    /// with types suitable for joining or aligning with encoded methylation data (e.g., boolean strand/context).
    ///
    /// Returns
    /// -------
    /// DataFrame
    ///     A Polars DataFrame with encoded context information.
    pub fn to_encoded_df(&self) -> PyResult<PyDataFrame> {
        let df = self.data.clone().to_df::<RsEncodedBsxBatch>();
        Ok(PyDataFrame(df))
    }
}