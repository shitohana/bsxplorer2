use log::info;
use polars::datatypes::{DataType, PlIndexMap, PlSmallStr};
use polars::io::csv::write::BatchedWriter;
use polars::io::mmap::MmapBytesReader;
use polars::prelude::*;
use std::io::Write;
use std::num::NonZero;
use std::ops::Div;

#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub enum ReportType {
    BISMARK,
    CGMAP,
    BEDGRAPH,
    COVERAGE,
}

impl ReportType {
    /// Get polars [Schema] of the report type
    pub fn get_schema(&self) -> Schema {
        match self {
            ReportType::BISMARK => Schema::from(PlIndexMap::from_iter([
                (PlSmallStr::from("chr"), DataType::String),
                (PlSmallStr::from("position"), DataType::UInt64),
                (PlSmallStr::from("strand"), DataType::String),
                (PlSmallStr::from("count_m"), DataType::UInt32),
                (PlSmallStr::from("count_um"), DataType::UInt32),
                (PlSmallStr::from("context"), DataType::String),
                (PlSmallStr::from("trinuc"), DataType::String),
            ])),
            ReportType::CGMAP => Schema::from(PlIndexMap::from_iter([
                (PlSmallStr::from("chr"), DataType::String),
                (PlSmallStr::from("nuc"), DataType::String),
                (PlSmallStr::from("position"), DataType::UInt64),
                (PlSmallStr::from("context"), DataType::String),
                (PlSmallStr::from("dinuc"), DataType::String),
                (PlSmallStr::from("count_m"), DataType::UInt32),
                (PlSmallStr::from("count_total"), DataType::UInt32),
            ])),
            ReportType::BEDGRAPH => Schema::from(PlIndexMap::from_iter([
                (PlSmallStr::from("chr"), DataType::String),
                (PlSmallStr::from("start"), DataType::UInt64),
                (PlSmallStr::from("end"), DataType::UInt64),
                (PlSmallStr::from("density"), DataType::Float64),
            ])),
            ReportType::COVERAGE => Schema::from(PlIndexMap::from_iter([
                (PlSmallStr::from("chr"), DataType::String),
                (PlSmallStr::from("start"), DataType::UInt64),
                (PlSmallStr::from("end"), DataType::UInt64),
                (PlSmallStr::from("density"), DataType::Float64),
                (PlSmallStr::from("count_m"), DataType::UInt32),
                (PlSmallStr::from("count_um"), DataType::UInt32),
            ])),
        }
    }

    /// Get preconfigured CSV reader for specified report type;
    ///
    /// # Defaults
    /// chunk_size = 10_000
    /// low_memory = false
    /// rechunk = true
    pub fn get_reader<R: MmapBytesReader>(
        &self,
        handle: R,
        chunk_size: Option<usize>,
        low_memory: Option<bool>,
        n_threads: Option<usize>,
        rechunk: Option<bool>,
    ) -> CsvReader<R> {
        let mut read_options = CsvReadOptions::default()
            // As no yet supported formats have headers
            .with_has_header(false)
            .with_chunk_size(chunk_size.unwrap_or(10_000))
            .with_low_memory(low_memory.unwrap_or(false))
            .with_n_threads(n_threads)
            .with_schema(Some(SchemaRef::from(self.get_schema())))
            .with_rechunk(rechunk.unwrap_or(true))
            .with_parse_options({
                CsvParseOptions::default()
                    .with_separator(b'\t')
                    .with_try_parse_dates(false)
                    .with_quote_char(Some(b'#'))
            });
        match self {
            ReportType::BEDGRAPH => {
                read_options = read_options.with_skip_rows(1);
            }
            _ => {}
        };
        read_options.into_reader_with_file_handle(handle)
    }

    pub fn get_writer<W: Write>(
        &self,
        sink: W,
        batch_size: Option<usize>,
        n_threads: Option<usize>,
    ) -> BatchedWriter<W> {
        let writer = match self {
            ReportType::COVERAGE
            | ReportType::BEDGRAPH
            | ReportType::BISMARK
            | ReportType::CGMAP => CsvWriter::new(sink)
                .with_separator(b'\t')
                .with_batch_size(NonZero::new(batch_size.unwrap_or_default()).unwrap())
                .include_header(false)
                .n_threads(n_threads.unwrap_or_default())
                .batched(&self.get_schema())
                .expect("could not create the csv writer"),
        };
        info!("Opened writer with {self:?} with batch_size {batch_size:?} with n_threads {n_threads:?}");
        writer
    }

    /// Modify raw [LazyFrame] to Universal BSX schema
    ///
    /// # Returns
    ///
    /// [LazyFrame] with columns renamed/mutated to universal schema
    pub(crate) fn to_universal_mutate(&self, lazy_frame: LazyFrame) -> LazyFrame {
        let context_encoder = when(col("context").eq(lit("CG")))
            .then(lit(true))
            .when(col("context").eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean);

        match self {
            ReportType::BISMARK => lazy_frame
                .with_column((col("count_m") + col("count_um")).alias("count_total"))
                .with_columns([
                    (col("count_m") / col("count_total").cast(DataType::Float64)).cast(DataType::Float64).alias("density"),
                    col("strand").eq(lit("+")).alias("strand"),
                    col("count_m") / col("count_total").alias("context"),
                ]),
            ReportType::CGMAP => lazy_frame.with_columns([
                col("nuc").eq(lit("C")).alias("strand"),
                context_encoder.alias("context"),
            ]),
            ReportType::BEDGRAPH => lazy_frame
                .rename(["start"], ["position"], true)
                .drop(["end"])
                .with_columns([
                    lit(NULL).alias("count_m"),
                    lit(NULL).alias("count_total"),
                    col("density").div(lit(100)).alias("density"),
                ]),
            ReportType::COVERAGE => lazy_frame
                .rename(["start"], ["position"], true)
                .drop(["end"]),
        }
    }
}
