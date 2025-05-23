use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result};
use itertools::Itertools;
use log::{debug, info, warn};
use polars::datatypes::DataType;
use polars::error::PolarsResult;
use polars::export::arrow::datatypes::Metadata;
use polars::prelude::{IpcCompression, IpcWriterOptions, Schema};

use crate::data_structs::batch::{create_caregorical_dtype, BsxBatch, BsxColumns};

/// Writer for BSX data in Arrow IPC format with optional compression.
pub struct BsxFileWriter<W>
where
    W: Write, {
    writer: polars::io::ipc::BatchedWriter<W>,
    schema: Schema,
}

impl<W> BsxFileWriter<W>
where
    W: Write,
{
    /// Creates a new BSX IPC writer.
    pub fn try_new(
        sink: W,
        chr_names: Vec<String>,
        compression: Option<IpcCompression>,
        custom_metadata: Option<Arc<Metadata>>,
    ) -> Result<BsxFileWriter<W>> {
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
        let chr_dtype =
            create_caregorical_dtype(chr_names.into_iter().map(Some).collect_vec());
        let mut schema = BsxColumns::schema();
        schema.set_dtype(BsxColumns::Chr.as_str(), chr_dtype);

        let mut writer = opts.to_writer(sink);
        if let Some(metadata) = custom_metadata {
            writer.set_custom_schema_metadata(metadata);
        }

        let batched_writer = writer
            .batched(&schema)
            .with_context(|| "Failed to create batched writer")?;

        Ok(Self {
            writer: batched_writer,
            schema,
        })
    }

    /// Creates a writer from a sink and FASTA index file.
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

        let index = bio::io::fasta::Index::from_file(&fai_path).with_context(|| {
            format!("Failed to read FASTA index from {:?}", fai_path)
        })?;

        let chr_names = index
            .sequences()
            .into_iter()
            .map(|seq| seq.name)
            .collect_vec();

        Self::try_new(sink, chr_names, compression, custom_metadata).with_context(
            || format!("Failed to create writer from FASTA index at {:?}", fai_path),
        )
    }

    /// Creates a writer from a sink and FASTA file.
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

        let index_path = format!("{}.fai", fasta_path.to_str().unwrap());
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

    /// Writes an encoded BSX batch.
    pub fn write_batch(
        &mut self,
        batch: BsxBatch,
    ) -> PolarsResult<()> {
        self.writer.write_batch(batch.data())
    }

    /// Finalizes the IPC file and closes the writer.
    pub fn close(&mut self) -> PolarsResult<()> {
        info!("Closing BsxIpcWriter");
        self.writer.finish()
    }

    /// Returns the chromosome data type.
    pub fn get_chr_dtype(&self) -> &DataType {
        // We know "chr" exists in the schema, as we created it in try_new
        self.schema.get("chr").unwrap()
    }
}

impl<W: Write> Drop for BsxFileWriter<W> {
    /// Ensures resources are properly cleaned up when dropped.
    fn drop(&mut self) {
        if let Err(e) = self.close() {
            warn!("Error closing BsxIpcWriter during drop: {}", e);
        }
    }
}
