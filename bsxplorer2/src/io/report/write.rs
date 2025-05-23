use std::io::{Seek, Write};

use anyhow::anyhow;
use log::{debug, info, warn};
use polars::io::csv::write::{BatchedWriter as BatchedCsvWriter, CsvWriter};
use polars::prelude::*;

use crate::data_structs::batch::BsxBatch;
#[cfg(feature = "compression")]
use crate::io::compression::Compression;
use crate::io::report::schema::ReportType;

/// Writes report data to a sink in CSV format based on a specified schema.
pub struct ReportWriter {
    /// The schema defining the structure of the report
    schema: ReportType,
    /// Batched CSV writer that handles the actual writing
    writer: BatchedCsvWriter<Box<dyn Write>>,
}

impl ReportWriter {
    /// Creates a new ReportWriter
    pub fn try_new<W: Write + Seek + 'static>(
        sink: W,
        schema: ReportType,
        n_threads: usize,
        #[cfg(feature = "compression")] compression: Compression,
        #[cfg(feature = "compression")] compression_level: Option<u32>,
    ) -> anyhow::Result<Self> {
        debug!("Creating new ReportWriter with {} threads", n_threads);
        let report_options = schema.read_options();

        #[cfg(feature = "compression")]
        let sink = compression.get_encoder(sink, compression_level.unwrap_or(1))?;
        #[cfg(not(feature = "compression"))]
        let sink = Box::new(sink) as Box<dyn Write>;

        let writer = CsvWriter::new(sink)
            .include_header(report_options.has_header)
            .with_separator(report_options.parse_options.separator)
            .with_quote_char(
                report_options.parse_options.quote_char.unwrap_or_default(),
            )
            .n_threads(n_threads)
            .batched(&schema.schema())
            .map_err(|e| {
                warn!("Failed to create batched CSV writer: {}", e);
                e
            })?;

        info!("ReportWriter successfully created");
        Ok(Self { schema, writer })
    }

    /// Writes a batch of data to the destination
    pub fn write_batch(
        &mut self,
        batch: BsxBatch,
    ) -> anyhow::Result<()> {
        let mut converted = batch.into_report(self.schema)?;

        converted.rechunk_mut();

        self.writer
            .write_batch(&converted)
            .map_err(|e| anyhow::anyhow!("Failed to write batch: {}", e))
    }

    /// Writes a DataFrame directly to the destination
    pub fn write_df(
        &mut self,
        df: &DataFrame,
    ) -> PolarsResult<()> {
        self.writer.write_batch(df).map_err(|e| {
            warn!("Failed to write DataFrame: {}", e);
            e
        })
    }

    pub fn finish(mut self) -> anyhow::Result<()> {
        self.writer.finish().map_err(|e| anyhow!(e))
    }
}
