use crate::data_structs::bsx_batch::BsxBatch;
use crate::io::report::schema::ReportTypeSchema;
use polars::io::csv::write::{BatchedWriter as BatchedCsvWriter, CsvWriter};
use polars::prelude::*;
use std::io::Write;

pub struct ReportWriter<W: Write> {
    schema: ReportTypeSchema,
    writer: BatchedCsvWriter<W>,
}

impl<W: Write> ReportWriter<W> {
    pub fn try_new(sink: W, schema: ReportTypeSchema, n_threads: usize) -> PolarsResult<Self> {
        let report_options = schema.read_options();

        let writer = CsvWriter::new(sink)
            .include_header(report_options.has_header)
            .with_separator(report_options.parse_options.separator)
            .with_quote_char(report_options.parse_options.quote_char.unwrap_or_default())
            .n_threads(n_threads)
            .batched(&schema.schema())?;

        Ok(Self { schema, writer })
    }

    pub fn write_batch(&mut self, batch: BsxBatch) -> PolarsResult<()> {
        let mut converted = self.schema.report_mutate_from_bsx(batch.into())?;
        converted.rechunk_mut();
        self.writer.write_batch(&converted)
    }

    fn write_df(&mut self, df: &DataFrame) -> PolarsResult<()> {
        self.writer.write_batch(df)
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

    #[pyclass(name = "ReportWriter")]
    pub struct PyReportWriter {
        inner: ReportWriter<BufWriter<File>>,
    }

    #[pymethods]
    impl PyReportWriter {
        #[new]
        pub fn new(sink: String, schema: ReportTypeSchema, n_threads: usize) -> PyResult<Self> {
            let file = BufWriter::new(File::create(sink)?);
            wrap_polars_result!(
                ReportWriter::<BufWriter<File>>::try_new(file, schema, n_threads)
                    .map(|v| Self { inner: v })
            )
        }

        pub fn write_batch(&mut self, batch: BsxBatch) -> PyResult<()> {
            wrap_polars_result!(self.inner.write_batch(batch))
        }

        pub fn write_df(&mut self, df: PyDataFrame) -> PyResult<()> {
            let df: DataFrame = df.into();
            wrap_polars_result!(self.inner.write_df(&df))
        }
    }
}
