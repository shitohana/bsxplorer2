use crate::data_structs::bsx_batch::{BsxBatch, BsxBatchMethods, EncodedBsxBatch};
use crate::utils::get_categorical_dtype;
#[cfg(feature = "python")]
use crate::utils::wrap_polars_result;
use anyhow::{Context, Result};
use itertools::Itertools;
use log::{debug, info, warn};
use polars::datatypes::DataType;
use polars::error::PolarsResult;
use polars::export::arrow::datatypes::Metadata;
use polars::prelude::{IpcCompression, IpcWriterOptions, Schema};
#[cfg(feature = "python")]
use pyo3_polars::error::PyPolarsErr;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;

pub use polars::prelude::IpcCompression as PolarsIpcCompression;

/// Writer for BSX data in Arrow IPC format.
/// Handles serialization of BSX batches to Arrow IPC format with optional compression.
pub struct BsxIpcWriter<W>
where
    W: Write,
{
    /// Underlying Arrow IPC writer that handles the actual serialization
    writer: polars::io::ipc::BatchedWriter<W>,
    /// Schema defining the structure of the BSX data
    schema: Schema,
}

impl<W> BsxIpcWriter<W>
where
    W: Write,
{
    /// Creates a new BSX IPC writer with the given sink, chromosome names, compression settings and metadata.
    ///
    /// # Arguments
    /// * `sink` - The destination to write the IPC data to
    /// * `chr_names` - List of chromosome names to use for categorical encoding
    /// * `compression` - Optional compression algorithm to use
    /// * `custom_metadata` - Optional custom metadata to include in the IPC file
    ///
    /// # Returns
    /// A new `BsxIpcWriter` or an error if initialization fails
    pub fn try_new(
        sink: W,
        chr_names: Vec<String>,
        compression: Option<IpcCompression>,
        custom_metadata: Option<Arc<Metadata>>,
    ) -> anyhow::Result<BsxIpcWriter<W>> {
        debug!(
            "Creating new BsxIpcWriter with {} chromosome names",
            chr_names.len()
        );
        if chr_names.is_empty() {
            warn!("Empty chromosome name list provided to BsxIpcWriter");
        }

        let opts = IpcWriterOptions {
            compression: compression.clone(),
            maintain_order: true,
        };

        info!("Initializing writer with compression: {:?}", compression);
        let chr_dtype = get_categorical_dtype(chr_names);
        let schema = EncodedBsxBatch::get_schema(&chr_dtype);

        let mut writer = opts.to_writer(sink);
        if let Some(metadata) = custom_metadata {
            debug!("Setting custom schema metadata");
            writer.set_custom_schema_metadata(metadata);
        }

        let batched_writer = writer
            .batched(&schema)
            .with_context(|| "Failed to create batched writer")?;

        debug!("Successfully created batched writer");
        Ok(Self {
            writer: batched_writer,
            schema,
        })
    }

    /// Creates a writer from a sink and FASTA index file path.
    ///
    /// This extracts chromosome names from the FASTA index file.
    ///
    /// # Arguments
    /// * `sink` - The destination to write the IPC data to
    /// * `fai_path` - Path to the FASTA index file (.fai)
    /// * `compression` - Optional compression algorithm to use
    /// * `custom_metadata` - Optional custom metadata to include in the IPC file
    ///
    /// # Returns
    /// A new `BsxIpcWriter` or an error if initialization fails
    pub fn try_from_sink_and_fai(
        sink: W,
        fai_path: PathBuf,
        compression: Option<IpcCompression>,
        custom_metadata: Option<Arc<Metadata>>,
    ) -> Result<Self> {
        info!(
            "Creating BsxIpcWriter from FASTA index file: {:?}",
            fai_path
        );

        let index = bio::io::fasta::Index::from_file(&fai_path)
            .with_context(|| format!("Failed to read FASTA index from {:?}", fai_path))?;

        let chr_names = index
            .sequences()
            .into_iter()
            .map(|seq| seq.name)
            .collect_vec();

        debug!(
            "Extracted {} chromosome names from FASTA index",
            chr_names.len()
        );
        Self::try_new(sink, chr_names, compression, custom_metadata)
            .with_context(|| format!("Failed to create writer from FASTA index at {:?}", fai_path))
    }

    /// Creates a writer from a sink and FASTA file path.
    ///
    /// This creates a FASTA index file if it doesn't exist, then extracts chromosome names.
    ///
    /// # Arguments
    /// * `sink` - The destination to write the IPC data to
    /// * `fasta_path` - Path to the FASTA file
    /// * `compression` - Optional compression algorithm to use
    /// * `custom_metadata` - Optional custom metadata to include in the IPC file
    ///
    /// # Returns
    /// A new `BsxIpcWriter` or an error if initialization fails
    pub fn try_from_sink_and_fasta(
        sink: W,
        fasta_path: PathBuf,
        compression: Option<IpcCompression>,
        custom_metadata: Option<Arc<Metadata>>,
    ) -> Result<Self> {
        info!("Creating BsxIpcWriter from FASTA file: {:?}", fasta_path);

        // Create index if it doesn't exist
        noodles::fasta::fs::index(fasta_path.clone())
            .with_context(|| format!("Failed to index FASTA file {:?}", fasta_path))?;

        debug!("FASTA index created or verified");
        let index_path = fasta_path.as_path().with_added_extension("fai");
        debug!("Using FASTA index at: {:?}", index_path);
        Self::try_from_sink_and_fai(sink, index_path, compression, custom_metadata)
            .with_context(|| format!("Failed to create writer from FASTA file {:?}", fasta_path))
    }

    /// Writes an encoded BSX batch to the output.
    ///
    /// # Arguments
    /// * `batch` - The encoded BSX batch to write
    ///
    /// # Returns
    /// Success or a Polars error if the write fails
    pub fn write_encoded_batch(&mut self, batch: EncodedBsxBatch) -> PolarsResult<()> {
        debug!("Writing encoded batch to IPC file");
        self.writer.write_batch(batch.data())
    }

    /// Encodes and writes a BSX batch to the output.
    ///
    /// # Arguments
    /// * `batch` - The BSX batch to encode and write
    ///
    /// # Returns
    /// Success or a Polars error if encoding or writing fails
    pub fn write_batch(&mut self, batch: BsxBatch) -> PolarsResult<()> {
        debug!("Encoding and writing batch to IPC file");
        let encoded = match EncodedBsxBatch::encode(batch, self.get_chr_dtype()) {
            Ok(encoded) => encoded,
            Err(e) => {
                warn!("Failed to encode BSX batch: {}", e);
                return Err(e);
            }
        };
        self.write_encoded_batch(encoded)
    }

    /// Finalizes the IPC file and closes the writer.
    ///
    /// # Returns
    /// Success or a Polars error if closing fails
    pub fn close(&mut self) -> PolarsResult<()> {
        info!("Closing BsxIpcWriter");
        self.writer.finish()
    }

    /// Returns the chromosome data type used by this writer.
    ///
    /// # Returns
    /// Reference to the chromosome DataType
    pub fn get_chr_dtype(&self) -> &DataType {
        // We know "chr" exists in the schema, as we created it in try_new
        self.schema.get("chr").unwrap()
    }
}

impl<W: Write> Drop for BsxIpcWriter<W> {
    /// Ensures resources are properly cleaned up when the writer is dropped.
    fn drop(&mut self) {
        if let Err(e) = self.close() {
            warn!("Error closing BsxIpcWriter during drop: {}", e);
        }
    }
}

#[cfg(feature = "python")]
mod python {
    use super::*;
    use crate::utils::wrap_box_result;
    use crate::wrap_polars_result;
    use pyo3::exceptions::PyIOError;
    use pyo3::prelude::*;
    use std::fs::File;
    use std::io::BufWriter;

    /// Python-exposed enumeration for IPC compression options.
    #[pyclass(name = "IpcCompression")]
    #[derive(Clone)]
    pub enum PyIpcCompression {
        /// ZSTD compression algorithm
        ZSTD,
        /// LZ4 compression algorithm
        LZ4,
    }

    impl From<PyIpcCompression> for IpcCompression {
        fn from(value: PyIpcCompression) -> Self {
            match value {
                PyIpcCompression::ZSTD => IpcCompression::ZSTD,
                PyIpcCompression::LZ4 => IpcCompression::LZ4,
            }
        }
    }

    /// Python-exposed BSX writer for IPC-formatted sequence data.
    #[pyclass(name = "BsxWriter")]
    pub struct PyBsxIpcWriter {
        /// The wrapped Rust writer implementation
        inner: BsxIpcWriter<BufWriter<File>>,
    }

    #[pymethods]
    impl PyBsxIpcWriter {
        /// Creates a new BSX writer for Python.
        ///
        /// # Arguments
        /// * `path` - Output file path for the IPC file
        /// * `fai_path` - Path to FASTA index file to extract chromosome names
        /// * `compression` - Optional compression algorithm to use
        ///
        /// # Returns
        /// A new PyBsxIpcWriter or raises a PyIOError on failure
        #[new]
        fn new(
            path: String,
            fai_path: String,
            compression: Option<PyIpcCompression>,
        ) -> PyResult<Self> {
            debug!(
                "Creating Python BsxWriter to {} using index {}",
                path, fai_path
            );

            // Create the file
            let file = match File::create(PathBuf::from(&path)) {
                Ok(f) => f,
                Err(e) => {
                    let error_msg = format!("Failed to create output file {}: {}", path, e);
                    warn!("{}", error_msg);
                    return Err(PyIOError::new_err(error_msg));
                }
            };

            let res = BsxIpcWriter::try_from_sink_and_fai(
                BufWriter::new(file),
                fai_path.into(),
                compression.map(|c| c.into()),
                None,
            );

            let inner = match wrap_box_result!(PyIOError, res) {
                Ok(inner) => inner,
                Err(e) => {
                    warn!("Failed to create BsxIpcWriter: {}", e);
                    return Err(e);
                }
            };

            Ok(Self { inner })
        }

        /// Writes a BSX batch to the output file.
        ///
        /// # Arguments
        /// * `batch` - The BSX batch to write
        ///
        /// # Returns
        /// None on success, or raises an exception on error
        fn write_batch(&mut self, batch: BsxBatch) -> PyResult<()> {
            debug!("Python interface: writing batch");
            wrap_polars_result!(self.inner.write_batch(batch))
        }

        /// Writes an encoded BSX batch to the output file.
        ///
        /// # Arguments
        /// * `batch` - The encoded BSX batch to write
        ///
        /// # Returns
        /// None on success, or raises an exception on error
        fn write_encoded_batch(&mut self, batch: EncodedBsxBatch) -> PyResult<()> {
            debug!("Python interface: writing encoded batch");
            wrap_polars_result!(self.inner.write_encoded_batch(batch))
        }

        /// Closes the writer and finalizes the output file.
        ///
        /// # Returns
        /// None on success, or raises an exception on error
        fn close(&mut self) -> PyResult<()> {
            info!("Python interface: closing BsxWriter");
            wrap_polars_result!(self.inner.close())
        }
    }
}

#[cfg(feature = "python")]
pub use python::{PyBsxIpcWriter, PyIpcCompression};
