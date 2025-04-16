use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;

use crate::data_structs::batch::BsxBatch;
use crate::data_structs::batch::EncodedBsxBatch;
use crate::data_structs::batch::{colnames, BsxBatchBuilder, BsxBatchMethods};
use crate::utils::get_categorical_dtype;
use anyhow::{Context, Result};
use itertools::Itertools;
use log::{debug, info, warn};
use polars::datatypes::DataType;
use polars::error::{PolarsError, PolarsResult};
use polars::export::arrow::datatypes::Metadata;
pub use polars::prelude::IpcCompression as PolarsIpcCompression;
use polars::prelude::{IpcCompression, IpcWriterOptions, Schema};

/// Writer for BSX data_structs in Arrow IPC format.
/// Handles serialization of BSX batches to Arrow IPC format with optional
/// compression.
pub struct BsxIpcWriter<W>
where
    W: Write,
{
    /// Underlying Arrow IPC writer that handles the actual serialization
    writer: polars::io::ipc::BatchedWriter<W>,
    /// Schema defining the structure of the BSX data_structs
    schema: Schema,
}

impl<W> BsxIpcWriter<W>
where
    W: Write,
{
    /// Creates a new BSX IPC writer with the given sink, chromosome names,
    /// compression settings and metadata.
    ///
    /// # Arguments
    /// * `sink` - The destination to write the IPC data_structs to
    /// * `chr_names` - List of chromosome names to use for categorical encoding
    /// * `compression` - Optional compression algorithm to use
    /// * `custom_metadata` - Optional custom metadata to include in the IPC
    ///   file
    ///
    /// # Returns
    /// A new `BsxIpcWriter` or an error if initialization fails
    pub fn try_new(
        sink: W,
        chr_names: Vec<String>,
        compression: Option<IpcCompression>,
        custom_metadata: Option<Arc<Metadata>>,
    ) -> Result<BsxIpcWriter<W>> {
        debug!(
            "Creating new BsxIpcWriter with {} chromosome names",
            chr_names.len()
        );
        if chr_names.is_empty() {
            warn!("Empty chromosome name list provided to BsxIpcWriter");
        }

        let opts = IpcWriterOptions {
            compression,
            maintain_order: true,
        };

        info!("Initializing writer with compression: {:?}", compression);
        let chr_dtype = get_categorical_dtype(chr_names);
        let mut schema = EncodedBsxBatch::schema();
        schema.set_dtype(colnames::CHR_NAME, chr_dtype);

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
    /// * `sink` - The destination to write the IPC data_structs to
    /// * `fai_path` - Path to the FASTA index file (.fai)
    /// * `compression` - Optional compression algorithm to use
    /// * `custom_metadata` - Optional custom metadata to include in the IPC
    ///   file
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

        let index =
            bio::io::fasta::Index::from_file(&fai_path).with_context(|| {
                format!("Failed to read FASTA index from {:?}", fai_path)
            })?;

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
            .with_context(|| {
                format!(
                    "Failed to create writer from FASTA index at {:?}",
                    fai_path
                )
            })
    }

    /// Creates a writer from a sink and FASTA file path.
    ///
    /// This creates a FASTA index file if it doesn't exist, then extracts
    /// chromosome names.
    ///
    /// # Arguments
    /// * `sink` - The destination to write the IPC data_structs to
    /// * `fasta_path` - Path to the FASTA file
    /// * `compression` - Optional compression algorithm to use
    /// * `custom_metadata` - Optional custom metadata to include in the IPC
    ///   file
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
        noodles::fasta::fs::index(fasta_path.clone()).with_context(|| {
            format!("Failed to index FASTA file {:?}", fasta_path)
        })?;

        debug!("FASTA index created or verified");
        let index_path = format!("{}.fai", fasta_path.to_str().unwrap());
        debug!("Using FASTA index at: {:?}", index_path);
        Self::try_from_sink_and_fai(
            sink,
            index_path.into(),
            compression,
            custom_metadata,
        )
        .with_context(|| {
            format!("Failed to create writer from FASTA file {:?}", fasta_path)
        })
    }

    /// Writes an encoded BSX batch to the output.
    ///
    /// # Arguments
    /// * `batch` - The encoded BSX batch to write
    ///
    /// # Returns
    /// Success or a Polars error if the write fails
    pub fn write_encoded_batch(
        &mut self,
        batch: EncodedBsxBatch,
    ) -> PolarsResult<()> {
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
    pub fn write_batch(
        &mut self,
        batch: BsxBatch,
    ) -> PolarsResult<()> {
        debug!("Encoding and writing batch to IPC file");
        let encoded = match BsxBatchBuilder::encode_batch(
            batch,
            self.get_chr_dtype().clone(),
        ) {
            Ok(encoded) => encoded,
            Err(_) => {
                return Err(PolarsError::ComputeError(
                    "failed to encode batch".into(),
                ));
            },
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

    /// Returns the chromosome data_structs type used by this writer.
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
