/// ***********************************************************************
/// *****
/// * Copyright (c) 2025
/// The Prosperity Public License 3.0.0
///
/// Contributor: [shitohana](https://github.com/shitohana)
///
/// Source Code: https://github.com/shitohana/BSXplorer
/// ***********************************************************************
/// ****

/// ***********************************************************************
/// *****
/// * Copyright (c) 2025
/// ***********************************************************************
/// ****
use std::io::Write;

use log::{debug, info, warn};
use polars::io::csv::write::{BatchedWriter as BatchedCsvWriter, CsvWriter};
use polars::prelude::*;

use crate::data_structs::bsx_batch::BsxBatch;
use crate::io::report::schema::ReportTypeSchema;

/// `ReportWriter` provides functionality to write report data_structs to a sink
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
    /// Creates a new `ReportWriter` with the specified sink, schema, and thread
    /// count.
    ///
    /// # Arguments
    ///
    /// * `sink` - The destination where the report will be written
    /// * `schema` - The schema defining the report structure and conversion
    ///   rules
    /// * `n_threads` - Number of threads to use for writing
    ///
    /// # Returns
    ///
    /// A `PolarsResult` containing the new `ReportWriter` or an error
    pub fn try_new(
        sink: W,
        schema: ReportTypeSchema,
        n_threads: usize,
    ) -> PolarsResult<Self> {
        debug!("Creating new ReportWriter with {} threads", n_threads);
        let report_options = schema.read_options();

        let writer = CsvWriter::new(sink)
            .include_header(report_options.has_header)
            .with_separator(report_options.parse_options.separator)
            .with_quote_char(
                report_options
                    .parse_options
                    .quote_char
                    .unwrap_or_default(),
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

    /// Writes a batch of data_structs to the destination.
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
    pub fn write_batch(
        &mut self,
        batch: BsxBatch,
    ) -> PolarsResult<()> {
        debug!("Converting BSX batch to report format");
        let mut converted = self
            .schema
            .report_mutate_from_bsx(batch.into())
            .map_err(|e| {
                warn!("Failed to convert BSX batch: {}", e);
                e
            })?;

        debug!("Rechunking converted data_structs for better performance");
        converted.rechunk_mut();

        debug!("Writing batch to destination");
        self.writer
            .write_batch(&converted)
            .map_err(|e| {
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
    pub fn write_df(
        &mut self,
        df: &DataFrame,
    ) -> PolarsResult<()> {
        debug!(
            "Writing DataFrame directly to destination, size: {}x{}",
            df.height(),
            df.width()
        );
        self.writer
            .write_batch(df)
            .map_err(|e| {
                warn!("Failed to write DataFrame: {}", e);
                e
            })
    }
}
