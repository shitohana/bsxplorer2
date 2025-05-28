use std::collections::HashMap;
use std::sync::Arc;

use bsxplorer2::data_structs::batch::{
    AggMethod,
    BsxBatch,
    BsxBatchBuilder,
    BsxColumns,
};
use bsxplorer2::data_structs::typedef::DensityType;
use bsxplorer2::data_structs::ContextData;
use bsxplorer2::tools::dimred::SegmentAlgorithm;
use bsxplorer2::utils::get_categorical_dtype;
use polars::prelude::{
    DataFrame,
    IntoSeries,
};
use pyo3::prelude::*;
use pyo3_polars::error::PyPolarsErr;
use pyo3_polars::{
    PyDataFrame,
    PyDataType,
    PySchema,
    PySeries,
};

use super::context_data::PyContextData;
use super::coords::{
    PyContig,
    PyGenomicPosition,
};
use super::lazy::PyLazyBsxBatch;
use super::report_schema::PyReportTypeSchema;
use super::stats::PyMethylationStats;
use super::utils::{
    PyContext,
    PyStrand,
};

#[pyclass(name = "BsxColumns", eq, eq_int)]
#[derive(Clone, PartialEq, Eq)]
pub enum PyBsxColumns {
    Chr,
    Position,
    Strand,
    Context,
    CountM,
    CountTotal,
    Density,
}

#[pymethods]
impl PyBsxColumns {
    #[staticmethod]
    pub fn schema() -> PyResult<PySchema> {
        Ok(PySchema(Arc::new(BsxColumns::schema())))
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            PyBsxColumns::Chr => BsxColumns::Chr.as_str(),
            PyBsxColumns::Position => BsxColumns::Position.as_str(),
            PyBsxColumns::Strand => BsxColumns::Strand.as_str(),
            PyBsxColumns::Context => BsxColumns::Context.as_str(),
            PyBsxColumns::CountM => BsxColumns::CountM.as_str(),
            PyBsxColumns::CountTotal => BsxColumns::CountTotal.as_str(),
            PyBsxColumns::Density => BsxColumns::Density.as_str(),
        }
    }

    #[staticmethod]
    pub fn colnames() -> [&'static str; 7] {
        BsxColumns::colnames()
    }

    #[staticmethod]
    pub fn has_name(name: &str) -> bool {
        BsxColumns::has_name(name)
    }

    pub fn __str__(&self) -> &'static str {
        self.as_str()
    }

    pub fn __repr__(&self) -> String {
        format!("BsxColumns.{}", self.as_str())
    }
}

impl From<PyBsxColumns> for BsxColumns {
    fn from(py_col: PyBsxColumns) -> Self {
        match py_col {
            PyBsxColumns::Chr => BsxColumns::Chr,
            PyBsxColumns::Position => BsxColumns::Position,
            PyBsxColumns::Strand => BsxColumns::Strand,
            PyBsxColumns::Context => BsxColumns::Context,
            PyBsxColumns::CountM => BsxColumns::CountM,
            PyBsxColumns::CountTotal => BsxColumns::CountTotal,
            PyBsxColumns::Density => BsxColumns::Density,
        }
    }
}

impl From<BsxColumns> for PyBsxColumns {
    fn from(col: BsxColumns) -> Self {
        match col {
            BsxColumns::Chr => PyBsxColumns::Chr,
            BsxColumns::Position => PyBsxColumns::Position,
            BsxColumns::Strand => PyBsxColumns::Strand,
            BsxColumns::Context => PyBsxColumns::Context,
            BsxColumns::CountM => PyBsxColumns::CountM,
            BsxColumns::CountTotal => PyBsxColumns::CountTotal,
            BsxColumns::Density => PyBsxColumns::Density,
        }
    }
}

#[pyclass(name = "AggMethod", eq, eq_int)]
#[derive(Clone, PartialEq, Eq)]
pub enum PyAggMethod {
    Mean,
    GeometricMean,
    Median,
    Max,
    Min,
}

impl From<PyAggMethod> for AggMethod {
    fn from(py_agg: PyAggMethod) -> Self {
        match py_agg {
            PyAggMethod::Mean => AggMethod::Mean,
            PyAggMethod::GeometricMean => AggMethod::GeometricMean,
            PyAggMethod::Median => AggMethod::Median,
            PyAggMethod::Max => AggMethod::Max,
            PyAggMethod::Min => AggMethod::Min,
        }
    }
}

#[pyclass(name = "BsxBatch")]
#[derive(Debug, Clone)]
pub struct PyBsxBatch {
    inner: BsxBatch,
}

impl From<BsxBatch> for PyBsxBatch {
    fn from(inner: BsxBatch) -> Self {
        PyBsxBatch { inner }
    }
}
impl From<PyBsxBatch> for BsxBatch {
    fn from(py_batch: PyBsxBatch) -> Self {
        py_batch.inner
    }
}

#[pymethods]
impl PyBsxBatch {
    /// Creates a new BsxBatch from columnar data.
    /// All vectors must have the same length.
    #[new]
    #[pyo3(signature = (chr, chr_dtype, positions, strand, context, count_m, count_total))]
    pub fn new(
        chr: String,
        chr_dtype: Option<PyDataType>, /* Can optionally provide a specific
                                        * categorical dtype */
        positions: Vec<u32>,
        strand: Vec<bool>,
        context: Vec<Option<bool>>,
        count_m: Vec<u16>,
        count_total: Vec<u16>,
    ) -> PyResult<Self> {
        let chr_dtype = chr_dtype.map(|s| s.0);

        let inner = BsxBatch::try_from_columns(
            &chr,
            chr_dtype,
            positions,
            strand,
            context,
            count_m,
            count_total,
        )
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(PyBsxBatch { inner })
    }

    /// Creates a BsxBatch from a Polars DataFrame after validation and casting.
    #[staticmethod]
    #[pyo3(signature = (
        data,
        check_nulls=true,
        check_sorted=true,
        check_duplicates=true,
        rechunk=true,
        check_single_chr=true,
        chr_values=None // Provides values for the categorical dtype if needed
    ))]
    pub fn from_dataframe(
        data: PyDataFrame,
        check_nulls: bool,
        check_sorted: bool,
        check_duplicates: bool,
        rechunk: bool,
        check_single_chr: bool,
        chr_values: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let mut df: DataFrame = data.into();

        let builder = BsxBatchBuilder::no_checks()
            .with_check_nulls(check_nulls)
            .with_check_sorted(check_sorted)
            .with_check_duplicates(check_duplicates)
            .with_rechunk(rechunk)
            .with_check_single_chr(check_single_chr)
            .with_chr_dtype(chr_values.map(get_categorical_dtype));

        // First cast to ensure columns exist and have correct types before checks
        let mut casted_df = builder
            .cast_only(df)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        // Perform checks on the casted DataFrame
        builder
            .checks_only(&casted_df)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        if rechunk {
            casted_df.rechunk_mut()
        }

        // Safety: The checks should guarantee the DataFrame conforms to the BsxBatch
        // schema
        Ok(unsafe { BsxBatch::new_unchecked(casted_df) }.into())
    }

    #[staticmethod]
    pub fn schema() -> PyResult<PySchema> {
        Ok(PySchema(Arc::new(BsxColumns::schema())))
    }

    #[staticmethod]
    #[pyo3(signature = (chr_dtype=None))]
    pub fn empty(chr_dtype: Option<PyDataType>) -> PyResult<Self> {
        let chr_dtype_option = chr_dtype.map(|s| s.0);
        Ok(Self {
            inner: BsxBatch::empty(chr_dtype_option.as_ref()),
        })
    }

    #[staticmethod]
    pub fn concat(batches: Vec<PyBsxBatch>) -> PyResult<Self> {
        let rust_batches: Vec<BsxBatch> =
            batches.into_iter().map(|b| b.inner).collect();
        let result = BsxBatchBuilder::concat(rust_batches)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(result.into())
    }

    // COLUMN GETTERS (return PySeries)
    pub fn chr(&self) -> PySeries {
        PySeries(self.inner.chr().clone().into_series())
    }

    pub fn position(&self) -> PySeries {
        PySeries(self.inner.position().clone().into_series())
    }

    pub fn strand(&self) -> PySeries {
        PySeries(self.inner.strand().clone().into_series())
    }

    pub fn context(&self) -> PySeries {
        PySeries(self.inner.context().clone().into_series())
    }

    pub fn count_m(&self) -> PySeries {
        PySeries(self.inner.count_m().clone().into_series())
    }

    pub fn count_total(&self) -> PySeries {
        PySeries(self.inner.count_total().clone().into_series())
    }

    pub fn density(&self) -> PySeries {
        PySeries(self.inner.density().clone().into_series())
    }

    pub fn column(
        &self,
        name: PyBsxColumns,
    ) -> PySeries {
        PySeries(self.inner.column(name.into()).clone().into_series())
    }

    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    pub fn split_at(
        &self,
        index: usize,
    ) -> (Self, Self) {
        // clone needed because split_at takes Self, not &Self
        let (batch1, batch2) = self.inner.clone().split_at(index);
        (PyBsxBatch { inner: batch1 }, PyBsxBatch { inner: batch2 })
    }

    pub fn rechunk(&mut self) -> Self {
        self.inner.rechunk();
        self.clone()
    }

    pub fn data(&self) -> PyResult<PyDataFrame> {
        Ok(PyDataFrame(self.inner.data().clone()))
    }

    /// Consumes the batch and returns the inner DataFrame.
    pub fn into_dataframe(&self) -> PyResult<PyDataFrame> {
        Ok(PyDataFrame(self.inner.clone().into_inner()))
    }

    pub fn slice(
        &self,
        start: i64,
        length: usize,
    ) -> Self {
        self.inner.slice(start, length).into()
    }

    pub fn add_context_data(
        &self,
        context_data: PyContextData,
    ) -> PyResult<Self> {
        // clone needed because add_context_data takes Self, not &Self
        let context_data: ContextData = context_data.into();
        let result_batch = self
            .inner
            .clone()
            .add_context_data(context_data)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(result_batch.into())
    }

    pub fn extend(
        &mut self,
        other: &Self,
    ) -> PyResult<()> {
        self.inner
            .extend(&other.inner)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    pub fn extend_unchecked(
        &mut self,
        other: &Self,
    ) -> PyResult<()> {
        unsafe {
            self.inner
                .extend_unchecked(&other.inner)
                .map_err(|e| PyPolarsErr::Polars(e).into())
        }
    }

    #[pyo3(signature = (min_size, beta=None))]
    pub fn shrink(
        &self,
        min_size: usize,
        beta: Option<f64>,
    ) -> PyResult<(Vec<f64>, Vec<f64>)> {
        let result = self.inner.shrink(SegmentAlgorithm::Pelt(beta, min_size))?;
        Ok(result)
    }

    pub fn discretise(
        &self,
        n_fragments: usize,
        method: PyAggMethod,
    ) -> PyResult<(Vec<f64>, Vec<f64>)> {
        let result = self.inner.discretise(n_fragments, method.into())?;
        Ok(result)
    }

    pub fn partition(
        &self,
        breakpoints: Vec<usize>,
        method: PyAggMethod,
    ) -> PyResult<(Vec<f64>, Vec<f64>)> {
        let result = self
            .inner
            .partition(breakpoints, AggMethod::from(method).get_fn())?;
        Ok(result)
    }

    // Position / Contig methods
    pub fn seqname(&self) -> Option<String> {
        self.inner.seqname().map(String::from)
    }

    pub fn first_pos(&self) -> Option<u32> {
        self.inner.first_pos()
    }

    pub fn last_pos(&self) -> Option<u32> {
        self.inner.last_pos()
    }

    pub fn first_genomic_pos(&self) -> Option<PyGenomicPosition> {
        self.inner.first_genomic_pos().map(|pos| pos.into())
    }

    pub fn last_genomic_pos(&self) -> Option<PyGenomicPosition> {
        self.inner.last_genomic_pos().map(|pos| pos.into())
    }

    pub fn as_contig(&self) -> Option<PyContig> {
        self.inner.as_contig().map(|contig| contig.into())
    }

    pub fn get_methylation_stats(&self) -> PyResult<PyMethylationStats> {
        let rust_stats = self.inner.get_methylation_stats();
        let py_stats = rust_stats.into();
        Ok(py_stats)
    }

    pub fn get_coverage_dist(&self) -> PyResult<HashMap<u16, u32>> {
        Ok(self.inner.get_coverage_dist().into_iter().collect())
    }

    pub fn get_context_stats(
        &self
    ) -> PyResult<HashMap<PyContext, (DensityType, u32)>> {
        let rust_stats = self.inner.get_context_stats();
        let py_stats = rust_stats.into_iter().map(|(k, v)| (k.into(), v)).collect();
        Ok(py_stats)
    }

    pub fn get_strand_stats(&self) -> PyResult<HashMap<PyStrand, (DensityType, u32)>> {
        let rust_stats = self.inner.get_strand_stats();
        let py_stats = rust_stats.into_iter().map(|(k, v)| (k.into(), v)).collect();
        Ok(py_stats)
    }

    pub fn as_binom(
        &self,
        mean: f64,
        pvalue: f64,
    ) -> PyResult<Self> {
        // clone needed because as_binom takes Self, not &Self
        let binom_batch = self
            .inner
            .clone()
            .as_binom(mean, pvalue)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(binom_batch.into())
    }

    // Report conversion
    pub fn into_report(
        &self,
        report_type: PyReportTypeSchema,
    ) -> PyResult<PyDataFrame> {
        let batch = self.inner.clone(); // Clone the batch to avoid consuming it
        let report_type = report_type.to_rust();
        batch
            .into_report(report_type)
            .map(|df| PyDataFrame(df))
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    // Lazy conversion
    pub fn lazy(&self) -> PyResult<PyLazyBsxBatch> {
        Ok(self.inner.clone().lazy().into())
    }

    // Python special methods
    pub fn height(&self) -> usize {
        self.inner.len()
    }

    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    pub fn __repr__(&self) -> String {
        // Handle Option<&str> for seqname
        let seqname_str = self.inner.seqname().unwrap_or("<empty>");
        format!(
            "BsxBatch(shape={}, chr='{}')",
            self.inner.data().shape().1,
            seqname_str
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


// Removed PyEncodedBsxBatch, encode, and decode functions as they have no
// corresponding Rust types/methods in the provided context.
