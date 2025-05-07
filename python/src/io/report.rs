use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use bsxplorer2::data_structs::batch::BsxBatch;
use bsxplorer2::exports::polars::frame::DataFrame;
use bsxplorer2::io::compression::Compression;
use bsxplorer2::io::report::{ReportReader as RustReportReader,
                             ReportReaderBuilder,
                             ReportWriter as RustReportWriter};
use pyo3::exceptions::{PyFileNotFoundError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

use super::compression::PyCompression;
use crate::types::batch::PyBsxBatch;
use crate::types::report_schema::PyReportTypeSchema;

#[pyclass(name = "ReportReader", unsendable)]
pub struct PyReportReader {
    reader: Option<RustReportReader<BsxBatch>>,
}

#[pymethods]
impl PyReportReader {
    #[new]
    #[pyo3(signature = (
        path,
        report_type,
        chunk_size=10000,
        fasta_path=None,
        fai_path=None,
        batch_size=100000,
        n_threads=None,
        low_memory=false,
        queue_len=1000,
        compression=None,
    ))]
    fn new(
        path: PathBuf,
        report_type: PyReportTypeSchema,
        chunk_size: usize,
        fasta_path: Option<PathBuf>,
        fai_path: Option<PathBuf>,
        batch_size: usize,
        n_threads: Option<usize>,
        low_memory: bool,
        queue_len: usize,
        compression: Option<PyCompression>,
    ) -> PyResult<Self> {
        let mut builder = ReportReaderBuilder::<BsxBatch>::default()
            .with_report_type(report_type.to_rust())
            .with_chunk_size(chunk_size)
            .with_batch_size(batch_size)
            .with_low_memory(low_memory)
            .with_queue_len(queue_len);

        if let Some(fp) = fasta_path {
            builder = builder.with_fasta_path(fp);
        }
        if let Some(fai) = fai_path {
            builder = builder.with_fai_path(fai);
        }
        if let Some(n) = n_threads {
            builder = builder.with_n_threads(n);
        }

        {
            if let Some(comp) = compression {
                let comp_enum = Compression::from(comp);
                builder = builder.with_compression(comp_enum);
            }
        }

        let reader = builder.build(path).map_err(|e| {
            // Attempt to classify error
            if e.to_string()
                .contains("No such file or directory")
            {
                PyFileNotFoundError::new_err(e.to_string())
            }
            else if e
                .to_string()
                .contains("must be specified")
            {
                PyValueError::new_err(e.to_string())
            }
            else {
                PyRuntimeError::new_err(format!(
                    "Failed to build ReportReader: {}",
                    e
                ))
            }
        })?;

        Ok(Self {
            reader: Some(reader),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> { slf }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<PyBsxBatch>> {
        match slf.reader.as_mut() {
            Some(reader) => {
                match reader.next() {
                    Some(Ok(batch)) => Ok(Some(PyBsxBatch::from(batch))),
                    Some(Err(e)) => {
                        Err(PyRuntimeError::new_err(format!(
                            "Error reading next batch: {}",
                            e
                        )))
                    },
                    None => {
                        slf.reader = None; // Consume the reader
                        Ok(None) // Signals Python StopIteration
                    },
                }
            },
            None => Ok(None), // Already consumed or failed previously
        }
    }
}

#[pyclass(name = "ReportWriter", unsendable)]
pub struct PyReportWriter {
    writer: Option<RustReportWriter>,
}

#[pymethods]
impl PyReportWriter {
    #[new]
    #[pyo3(signature = (
        path,
        schema,
        n_threads = 1,
        compression = None,
        compression_level = None
    ))]
    fn new(
        path: PathBuf,
        schema: PyReportTypeSchema,
        n_threads: usize,
        compression: Option<PyCompression>,
        compression_level: Option<u32>,
    ) -> PyResult<Self> {
        let file = File::create(&path).map_err(|e| {
            PyRuntimeError::new_err(format!(
                "Failed to create output file '{}': {}",
                path.display(),
                e
            ))
        })?;
        let sink = BufWriter::new(file);

        let comp_enum = Compression::from(
            compression.unwrap_or_else(|| PyCompression(Compression::None)),
        );

        let writer = RustReportWriter::try_new(
            sink,
            schema.to_rust(),
            n_threads,
            comp_enum,
            compression_level,
        )
        .map_err(|e| {
            PyRuntimeError::new_err(format!(
                "Failed to create ReportWriter: {}",
                e
            ))
        })?;

        Ok(Self {
            writer: Some(writer),
        })
    }

    pub fn write_batch(
        &mut self,
        batch: PyBsxBatch,
    ) -> PyResult<()> {
        if let Some(writer) = self.writer.as_mut() {
            writer
                .write_batch(batch.into()) // Clone might be necessary depending on ownership
                .map_err(|e| PyRuntimeError::new_err(format!("Failed to write batch: {}", e)))
        }
        else {
            Err(PyRuntimeError::new_err(
                "Writer is closed or uninitialized.",
            ))
        }
    }

    pub fn write_df(
        &mut self,
        df: PyDataFrame,
    ) -> PyResult<()> {
        if let Some(writer) = self.writer.as_mut() {
            let rust_df: DataFrame = df.into();
            writer.write_df(&rust_df).map_err(|e| {
                PyRuntimeError::new_err(format!(
                    "Failed to write DataFrame: {}",
                    e
                ))
            })
        }
        else {
            Err(PyRuntimeError::new_err(
                "Writer is closed or uninitialized.",
            ))
        }
    }

    pub fn close(&mut self) -> PyResult<()> {
        if let Some(writer) = self.writer.take() {
            // The underlying BatchedCsvWriter flushes on drop, which happens
            // when `writer` goes out of scope here. We might add an
            // explicit finish/flush call if the Rust struct exposes one later.
            drop(writer);
            Ok(())
        }
        else {
            // Already closed, maybe warn or just do nothing? Let's return Ok
            // for idempotency.
            Ok(())
        }
    }
}
