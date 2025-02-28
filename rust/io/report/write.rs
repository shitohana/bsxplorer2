use crate::data_structs::bsx_batch::BsxBatch;
use crate::io::report::schema::ReportTypeSchema;
use log::{debug, info, warn};
use polars::io::csv::write::{BatchedWriter as BatchedCsvWriter, CsvWriter};
use polars::prelude::*;
use std::io::Write;

/// `ReportWriter` provides functionality to write report data to a sink
/// in CSV format based on a specified schema.
///
/// # Type Parameters
///
/// * `W`: Any type that implements the `Write` trait
pub struct ReportWriter<W: Write> {
    /// The schema defining the structure of the report
    schema: ReportTypeSchema,
    /// Batched CSV writer that handles the actual writing
    writer: BatchedCsvWriter<W>,
}

impl<W: Write> ReportWriter<W> {
    /// Creates a new `ReportWriter` with the specified sink, schema, and thread count.
    ///
    /// # Arguments
    ///
    /// * `sink` - The destination where the report will be written
    /// * `schema` - The schema defining the report structure and conversion rules
    /// * `n_threads` - Number of threads to use for writing
    ///
    /// # Returns
    ///
    /// A `PolarsResult` containing the new `ReportWriter` or an error
    pub fn try_new(sink: W, schema: ReportTypeSchema, n_threads: usize) -> PolarsResult<Self> {
        debug!("Creating new ReportWriter with {} threads", n_threads);
        let report_options = schema.read_options();

        let writer = CsvWriter::new(sink)
            .include_header(report_options.has_header)
            .with_separator(report_options.parse_options.separator)
            .with_quote_char(report_options.parse_options.quote_char.unwrap_or_default())
            .n_threads(n_threads)
            .batched(&schema.schema())
            .map_err(|e| {
                warn!("Failed to create batched CSV writer: {}", e);
                e
            })?;

        info!("ReportWriter successfully created");
        Ok(Self { schema, writer })
    }

    /// Writes a batch of data to the destination.
    ///
    /// This method converts the provided BSX batch to the report format
    /// according to the schema and writes it.
    ///
    /// # Arguments
    ///
    /// * `batch` - The BSX batch to write
    ///
    /// # Returns
    ///
    /// A `PolarsResult` indicating success or containing an error
    pub fn write_batch(&mut self, batch: BsxBatch) -> PolarsResult<()> {
        debug!("Converting BSX batch to report format");
        let mut converted = self
            .schema
            .report_mutate_from_bsx(batch.into())
            .map_err(|e| {
                warn!("Failed to convert BSX batch: {}", e);
                e
            })?;

        debug!("Rechunking converted data for better performance");
        converted.rechunk_mut();

        debug!("Writing batch to destination");
        self.writer.write_batch(&converted).map_err(|e| {
            warn!("Failed to write batch: {}", e);
            e
        })
    }

    /// Writes a DataFrame directly to the destination.
    ///
    /// # Arguments
    ///
    /// * `df` - The DataFrame to write
    ///
    /// # Returns
    ///
    /// A `PolarsResult` indicating success or containing an error
    pub fn write_df(&mut self, df: &DataFrame) -> PolarsResult<()> {
        debug!(
            "Writing DataFrame directly to destination, size: {}x{}",
            df.height(),
            df.width()
        );
        self.writer.write_batch(df).map_err(|e| {
            warn!("Failed to write DataFrame: {}", e);
            e
        })
    }
}

#[cfg(feature = "python")]
pub use python::PyReportWriter;
#[cfg(feature = "python")]
mod python {
    use super::*;
    use crate::utils::wrap_polars_result;
    use pyo3::prelude::*;
    use pyo3_polars::error::PyPolarsErr;
    use pyo3_polars::PyDataFrame;
    use std::fs::File;
    use std::io::BufWriter;

    /// Python-compatible wrapper for the Rust `ReportWriter`.
    ///
    /// This class provides a Python interface to the underlying Rust implementation,
    /// enabling writing report data to files from Python code.
    #[pyclass(name = "ReportWriter")]
    pub struct PyReportWriter {
        /// The inner Rust ReportWriter implementation
        inner: ReportWriter<BufWriter<File>>,
    }

    #[pymethods]
    impl PyReportWriter {
        /// Creates a new Python-compatible ReportWriter.
        ///
        /// # Arguments
        ///
        /// * `sink` - File path where the report will be written
        /// * `schema` - The schema defining the report structure
        /// * `n_threads` - Number of threads to use for writing
        ///
        /// # Returns
        ///
        /// A PyResult containing the new PyReportWriter or an error
        #[new]
        pub fn new(sink: String, schema: ReportTypeSchema, n_threads: usize) -> PyResult<Self> {
            info!("Creating new PyReportWriter for file: {}", sink);
            let file = BufWriter::new(File::create(&sink).map_err(|e| {
                let err_msg = format!("Failed to create file '{}': {}", sink, e);
                warn!("{}", err_msg);
                PyErr::new::<pyo3::exceptions::PyIOError, _>(err_msg)
            })?);

            wrap_polars_result!(
                ReportWriter::<BufWriter<File>>::try_new(file, schema, n_threads).map(|v| {
                    info!("PyReportWriter successfully created");
                    Self { inner: v }
                })
            )
        }

        /// Writes a batch of data to the file.
        ///
        /// # Arguments
        ///
        /// * `batch` - The BSX batch to write
        ///
        /// # Returns
        ///
        /// A PyResult indicating success or containing an error
        pub fn write_batch(&mut self, batch: BsxBatch) -> PyResult<()> {
            debug!("PyReportWriter: Writing BSX batch");
            wrap_polars_result!(self.inner.write_batch(batch))
        }

        /// Writes a DataFrame directly to the file.
        ///
        /// # Arguments
        ///
        /// * `df` - The DataFrame to write
        ///
        /// # Returns
        ///
        /// A PyResult indicating success or containing an error
        pub fn write_df(&mut self, df: PyDataFrame) -> PyResult<()> {
            debug!("PyReportWriter: Writing DataFrame");
            let df: DataFrame = df.into();
            wrap_polars_result!(self.inner.write_df(&df))
        }
    }
}
