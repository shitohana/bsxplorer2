use std::sync::Arc;

use bsxplorer2::exports::polars::error::PolarsError;
use bsxplorer2::exports::polars::frame::DataFrame;
use bsxplorer2::exports::polars::prelude::BooleanChunked;
use bsxplorer2::exports::polars::series::{IntoSeries, Series};
use pyo3::prelude::*;
use pyo3::types::{PyType};
use pyo3_polars::{PyDataFrame, PySchema, PySeries};
use bsxplorer2::data_structs::batch::{BsxBatchMethods, EncodedBsxBatch as RsEncodedBsxBatch};
use bsxplorer2::data_structs::batch::BsxBatchBuilder as RsBsxBatchBuilder;
use bsxplorer2::exports::polars::datatypes::CategoricalChunked;

use super::context_data::PyContextData;
use super::report_schema::PyReportTypeSchema;

#[pyclass(name="EncodedBsxBatch")]
#[derive(Clone)]
pub struct PyEncodedBsxBatch {
    pub inner: RsEncodedBsxBatch
}

impl From<RsEncodedBsxBatch> for PyEncodedBsxBatch {
    fn from(batch: RsEncodedBsxBatch) -> Self {
        PyEncodedBsxBatch { inner: batch }
    }
}

impl From<PyEncodedBsxBatch> for RsEncodedBsxBatch {
    fn from(batch: PyEncodedBsxBatch) -> Self {
        batch.inner
    }
}

#[pymethods]
impl PyEncodedBsxBatch {
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
    /// Create a new EncodedBsxBatch from a DataFrame.
    ///
    /// Performs various checks to ensure data integrity unless disabled.
    ///
    /// Parameters
    /// ----------
    /// data : DataFrame
    ///     Input Polars DataFrame containing BS-Seq data.
    /// check_nulls : bool, optional
    ///     If True (default), check for null values in essential columns.
    /// check_sorted : bool, optional
    ///     If True (default), check if the data is sorted by position.
    /// check_duplicates : bool, optional
    ///     If True (default), check for duplicate entries (chr, pos, strand).
    /// rechunk : bool, optional
    ///     If True (default), rechunk the DataFrame for optimal performance.
    /// check_single_chr : bool, optional
    ///     If True (default), ensure the batch contains data for only one chromosome.
    /// context_data : ContextData, optional
    ///     Context data associated with the batch.
    /// report_schema : ReportTypeSchema, optional
    ///     Schema defining the report type.
    ///
    /// Returns
    /// -------
    /// EncodedBsxBatch
    ///     A new instance of EncodedBsxBatch.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If any of the enabled checks fail.
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
                RsBsxBatchBuilder::default().with_report_type(report_schema.to_rust())
            } else {
                RsBsxBatchBuilder::default()
            }
        }
            .with_check_nulls(check_nulls)
            .with_check_sorted(check_sorted)
            .with_check_duplicates(check_duplicates)
            .with_rechunk(rechunk)
            .with_check_single_chr(check_single_chr)
            .with_context_data(context_data.map(|v| v.into()));

        let inner = builder.build::<RsEncodedBsxBatch>(df)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(PyEncodedBsxBatch { inner })
    }

    #[classmethod]
    /// Create a new EncodedBsxBatch from a DataFrame without performing any checks.
    ///
    /// Warning: Use with caution, assumes the input DataFrame is valid and correctly formatted.
    ///
    /// Parameters
    /// ----------
    /// _cls : type
    ///     The class object.
    /// data : DataFrame
    ///     Input Polars DataFrame.
    ///
    /// Returns
    /// -------
    /// EncodedBsxBatch
    ///     A new instance of EncodedBsxBatch.
    pub fn from_dataframe_unchecked(_cls: &Bound<'_, PyType>, data: PyDataFrame) -> PyResult<Self> {
        let df: DataFrame = data.into();
        let inner = unsafe { RsEncodedBsxBatch::new_unchecked(df) };
        Ok(PyEncodedBsxBatch { inner })
    }

    #[staticmethod]
    /// Get the Polars schema for an EncodedBsxBatch.
    ///
    /// Returns
    /// -------
    /// Schema
    ///     The Polars schema definition.
    pub fn schema() -> PyResult<PySchema> {
        Ok(PySchema(Arc::new(RsEncodedBsxBatch::schema())))
    }

    #[staticmethod]
    /// Create an empty EncodedBsxBatch with the correct schema.
    ///
    /// Returns
    /// -------
    /// EncodedBsxBatch
    ///     An empty batch instance.
    pub fn empty() -> PyResult<Self> {
        Ok(PyEncodedBsxBatch { inner: RsEncodedBsxBatch::empty() })
    }

    /// Access the chromosome column as a Polars Series.
    ///
    /// Returns
    /// -------
    /// Series
    ///     The chromosome data (Categorical).
    #[getter]
    pub fn chr(&self) -> PyResult<PySeries> {
        Ok(PySeries(self.inner.chr().clone().into_series()))
    }

    /// Access the chromosome column as a Categorical Series.
    ///
    /// Returns
    /// -------
    /// Series
    ///     The chromosome data as a Polars Categorical Series.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the 'chr' column is not found or is not categorical.
    #[getter]
    pub fn chr_categorical(&self) -> PyResult<PySeries> {
        let categorical_chunked: &CategoricalChunked = self.inner.data()
            .column(bsxplorer2::data_structs::batch::colnames::CHR_NAME)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?
            .categorical()
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(PySeries(categorical_chunked.clone().into_series()))
    }

    /// Access the position column as a Polars Series.
    ///
    /// Returns
    /// -------
    /// Series
    ///     The position data (UInt32).
    #[getter]
    pub fn position(&self) -> PyResult<PySeries> {
        Ok(PySeries(self.inner.position().clone().into_series()))
    }

    /// Access the strand column as a Polars Series.
    ///
    /// Returns
    /// -------
    /// Series
    ///     The strand data (Categorical).
    #[getter]
    pub fn strand(&self) -> PyResult<PySeries> {
        Ok(PySeries(self.inner.strand().clone().into_series()))
    }

    /// Access the context column as a Polars Series.
    ///
    /// Returns
    /// -------
    /// Series
    ///     The context data (Categorical).
    #[getter]
    pub fn context(&self) -> PyResult<PySeries> {
        Ok(PySeries(self.inner.context().clone().into_series()))
    }

    /// Access the methylated counts column as a Polars Series.
    ///
    /// Returns
    /// -------
    /// Series
    ///     The methylated count data (UInt32).
    #[getter]
    pub fn count_m(&self) -> PyResult<PySeries> {
        Ok(PySeries(self.inner.count_m().clone().into_series()))
    }

    /// Access the total counts column as a Polars Series.
    ///
    /// Returns
    /// -------
    /// Series
    ///     The total count data (UInt32).
    #[getter]
    pub fn count_total(&self) -> PyResult<PySeries> {
        Ok(PySeries(self.inner.count_total().clone().into_series()))
    }

    /// Access the density column as a Polars Series.
    ///
    /// Returns
    /// -------
    /// Series
    ///     The density data (Float32).
    #[getter]
    pub fn density(&self) -> PyResult<PySeries> {
        Ok(PySeries(self.inner.density().clone().into_series()))
    }

    /// Check if the batch is empty.
    ///
    /// Returns
    /// -------
    /// bool
    ///     True if the batch contains no rows, False otherwise.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Split the batch into two at the given index.
    ///
    /// Parameters
    /// ----------
    /// index : int
    ///     The row index at which to split the batch.
    ///
    /// Returns
    /// -------
    /// tuple[EncodedBsxBatch, EncodedBsxBatch]
    ///     A tuple containing two new batches, the first with rows up to `index`,
    ///     the second with rows from `index` onwards.
    pub fn split_at(&self, index: usize) -> PyResult<(Self, Self)> {
        let (batch1, batch2) = self.inner.clone().split_at(index);
        Ok((PyEncodedBsxBatch { inner: batch1 }, PyEncodedBsxBatch { inner: batch2 }))
    }

    /// Get a reference to the underlying Polars DataFrame.
    ///
    /// Returns
    /// -------
    /// DataFrame
    ///     A clone of the internal DataFrame.
    #[getter]
    pub fn data(&self) -> PyResult<PyDataFrame> {
        Ok(PyDataFrame(self.inner.data().clone()))
    }

    /// Consume the batch and return the underlying Polars DataFrame.
    ///
    /// Returns
    /// -------
    /// DataFrame
    ///     The internal DataFrame.
    #[getter]
    pub fn take(&self) -> PyResult<PyDataFrame> {
        Ok(PyDataFrame(self.inner.clone().take()))
    }

    /// Get the unique chromosome value present in the batch.
    ///
    /// Assumes the batch contains data for only a single chromosome.
    ///
    /// Returns
    /// -------
    /// str
    ///     The chromosome name.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the batch is empty or contains multiple chromosomes.
    pub fn chr_val(&self) -> PyResult<String> {
        self.inner.chr_val()
            .map(|s| s.to_string())
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Get the starting genomic position in the batch.
    ///
    /// Returns
    /// -------
    /// int or None
    ///     The minimum position value, or None if the batch is empty.
    pub fn start_pos(&self) -> Option<u32> {
        self.inner.start_pos()
    }

    /// Get the ending genomic position in the batch.
    ///
    /// Returns
    /// -------
    /// int or None
    ///     The maximum position value, or None if the batch is empty.
    pub fn end_pos(&self) -> Option<u32> {
        self.inner.end_pos()
    }

    /// Vertically stack this batch with another batch.
    ///
    /// Creates a new batch containing rows from both batches. Performs checks
    /// on the combined data.
    ///
    /// Parameters
    /// ----------
    /// other : EncodedBsxBatch
    ///     The batch to stack underneath this one.
    ///
    /// Returns
    /// -------
    /// EncodedBsxBatch
    ///     A new batch containing the combined data.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the stacking operation results in invalid data (e.g., unsorted, duplicates).
    pub fn vstack(&self, other: &PyEncodedBsxBatch) -> PyResult<Self> {
        let new_inner = self.inner.vstack(&other.inner)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(PyEncodedBsxBatch { inner: new_inner })
    }

    /// Extend this batch by appending rows from another batch in-place.
    ///
    /// Modifies the current batch. Performs checks after extending.
    ///
    /// Parameters
    /// ----------
    /// other : EncodedBsxBatch
    ///     The batch whose rows will be appended.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the extend operation results in invalid data (e.g., unsorted, duplicates).
    pub fn extend(&mut self, other: &PyEncodedBsxBatch) -> PyResult<()> {
        self.inner.extend(&other.inner)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Filter the batch using a boolean mask Series.
    ///
    /// Creates a new batch containing only the rows where the mask is True.
    ///
    /// Parameters
    /// ----------
    /// mask : Series
    ///     A Polars Series of boolean values with the same length as the batch.
    ///
    /// Returns
    /// -------
    /// EncodedBsxBatch
    ///     A new batch containing the filtered rows.
    ///
    /// Raises
    /// ------
    /// TypeError
    ///     If the mask is not a boolean Series.
    /// ValueError
    ///     If the mask length does not match the batch height or other Polars errors occur.
    pub fn filter_mask(&self, mask: PySeries) -> PyResult<Self> {
        let mask_series: Series = mask.into();
        let bool_mask: &BooleanChunked = mask_series.bool()
            .map_err(|_| pyo3::exceptions::PyTypeError::new_err("Mask must be a boolean Series"))?;

        let new_inner = self.inner.filter_mask(bool_mask)
             .map_err(|e: PolarsError| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(PyEncodedBsxBatch { inner: new_inner })
    }

    /// Get the number of rows in the batch.
    ///
    /// Returns
    /// -------
    /// int
    ///     The height (number of rows) of the batch.
    pub fn height(&self) -> usize {
        self.inner.height()
    }

    /// Get the number of rows in the batch (implements Python's `len()`).
    ///
    /// Returns
    /// -------
    /// int
    ///     The height (number of rows) of the batch.
    pub fn __len__(&self) -> usize {
        self.inner.height()
    }

    /// Get a string representation of the batch.
    ///
    /// Returns
    /// -------
    /// str
    ///     A string showing shape and chromosome.
    pub fn __repr__(&self) -> String {
        format!("EncodedBsxBatch(shape={}, chr='{}')", self.inner.data().shape().1, self.inner.chr_val().unwrap_or("N/A"))
    }

    /// Compare two batches for equality (implements `==` and `!=`).
    ///
    /// Parameters
    /// ----------
    /// other : EncodedBsxBatch
    ///     The batch to compare against.
    /// op : CompareOp
    ///     The comparison operation (Eq or Ne).
    ///
    /// Returns
    /// -------
    /// bool
    ///     True if the batches are equal/unequal based on `op`, False otherwise.
    ///
    /// Raises
    /// ------
    /// NotImplementedError
    ///     If the comparison operation is not `==` or `!=`.
    pub fn __richcmp__(&self, other: &Self, op: pyo3::basic::CompareOp) -> PyResult<bool> {
        match op {
            pyo3::basic::CompareOp::Eq => Ok(self.inner == other.inner),
            pyo3::basic::CompareOp::Ne => Ok(self.inner != other.inner),
            _ => Err(pyo3::exceptions::PyNotImplementedError::new_err("Only == and != are supported for EncodedBsxBatch")),
        }
    }
}