use std::io::{
    Seek,
    Write,
};

use anyhow::anyhow;
use polars::io::csv::write::{
    BatchedWriter as BatchedCsvWriter,
    CsvWriter,
};
use polars::prelude::*;

use crate::data_structs::batch::BsxBatch;
#[cfg(feature = "compression")]
use crate::io::compression::Compression;
use crate::io::report::schema::ReportType;

/// Writes report data to a sink in CSV format based on a specified schema.
pub struct ReportWriter {
    /// The schema defining the structure of the report
    schema:        ReportType,
    /// Batched CSV writer that handles the actual writing
    writer:        BatchedCsvWriter<Box<dyn Write>>,
    finish_output: FinishOutput,
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
        let report_options = schema.read_options();

        #[cfg(feature = "compression")]
        let (sink, finish_output) =
            compression.get_report_encoder(sink, compression_level.unwrap_or(1))?;
        #[cfg(not(feature = "compression"))]
        let (sink, finish_output) = shared_output(sink, |mut w| w.flush());

        let writer = CsvWriter::new(sink)
            .include_header(report_options.has_header)
            .with_separator(report_options.parse_options.separator)
            .with_quote_char(
                report_options.parse_options.quote_char.unwrap_or_default(),
            )
            .n_threads(n_threads)
            .batched(&schema.schema())?;

        Ok(Self {
            schema,
            writer,
            finish_output,
        })
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
        self.writer.write_batch(df)
    }

    pub fn finish(mut self) -> anyhow::Result<()> {
        self.writer.finish().map_err(|e| anyhow!(e))?;
        drop(self.writer);
        (self.finish_output)().map_err(|e| anyhow!(e))
    }
}

pub(crate) type FinishOutput = Box<dyn FnOnce() -> std::io::Result<()>>;

/// Share the encoder with an explicit finalizer so CSV finishing cannot hide
/// compression footer or sink flush errors.
pub(crate) fn shared_output<W: std::io::Write + 'static>(
    writer: W,
    finish: impl FnOnce(W) -> std::io::Result<()> + 'static,
) -> (Box<dyn std::io::Write>, FinishOutput) {
    use std::cell::RefCell;
    use std::rc::Rc;

    struct SharedWriter<W>(Rc<RefCell<Option<W>>>);
    impl<W: std::io::Write> std::io::Write for SharedWriter<W> {
        fn write(
            &mut self,
            bytes: &[u8],
        ) -> std::io::Result<usize> {
            self.0
                .borrow_mut()
                .as_mut()
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::BrokenPipe,
                        "report output is closed",
                    )
                })?
                .write(bytes)
        }

        fn flush(&mut self) -> std::io::Result<()> {
            self.0
                .borrow_mut()
                .as_mut()
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::BrokenPipe,
                        "report output is closed",
                    )
                })?
                .flush()
        }
    }
    let shared = Rc::new(RefCell::new(Some(writer)));
    let handle = Box::new(SharedWriter(shared.clone()));
    let finalizer = Box::new(move || {
        finish(shared.borrow_mut().take().expect("output finalized once"))
    });
    (handle, finalizer)
}
