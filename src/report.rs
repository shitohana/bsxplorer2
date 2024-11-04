use crate::sequence::{scan_fasta, SeqData};
use polars::export::arrow::datatypes::IntegerType;
use polars::io::ipc::BatchedWriter;
use polars::prelude::*;
use std::fs::File;

/// Get schema of BSXplorer internal format
pub fn get_universal_schema(chr_names: Vec<&str>) -> Schema {
    let categories = DataType::Enum(
        Some({
            let mut cat_builder =
                CategoricalChunkedBuilder::new(Default::default(), 0, Default::default());
            for chr_name in chr_names {
                let _ = &cat_builder.append(Some(chr_name));
            }
            cat_builder
                .finish()
                .get_rev_map()
                .clone()
        }),
        CategoricalOrdering::Physical,
    );

    assert!(categories.contains_categoricals());

    let mut schema = Schema::default();
    schema.insert(PlSmallStr::from_str("chr"), categories);
    schema.insert(PlSmallStr::from_str("position"), DataType::UInt64);
    schema.insert(PlSmallStr::from_str("strand"), DataType::Boolean);
    schema.insert(PlSmallStr::from_str("context"), DataType::Boolean);
    schema.insert(PlSmallStr::from_str("count_m"), DataType::UInt16);
    schema.insert(PlSmallStr::from_str("count_total"), DataType::UInt16);
    schema.insert(PlSmallStr::from_str("density"), DataType::Float32);

    schema
}

/// Type of methylation report
pub enum ReportType {
    BEDGRAPH,
    COVERAGE,
    CGMAP,
    BISMARK,
}

impl ReportType {
    /// Get corresponding polars schema for report type
    pub fn schema(&self) -> Schema {
        match self {
            ReportType::BISMARK => {
                let mut schema = Schema::default();
                schema.insert(PlSmallStr::from_str("chr"), DataType::String);
                schema.insert(PlSmallStr::from_str("position"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("strand"), DataType::String);
                schema.insert(PlSmallStr::from_str("count_m"), DataType::UInt32);
                schema.insert(PlSmallStr::from_str("count_um"), DataType::UInt32);
                schema.insert(PlSmallStr::from_str("context"), DataType::String);
                schema.insert(PlSmallStr::from_str("trinuc"), DataType::String);
                schema
            }
            ReportType::COVERAGE => {
                let mut schema = Schema::default();
                schema.insert(PlSmallStr::from_str("chr"), DataType::String);
                schema.insert(PlSmallStr::from_str("start"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("end"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("density"), DataType::Float32);
                schema.insert(PlSmallStr::from_str("count_m"), DataType::UInt32);
                schema.insert(PlSmallStr::from_str("count_um"), DataType::UInt32);
                schema
            }
            ReportType::CGMAP => {
                let mut schema = Schema::default();
                schema.insert(PlSmallStr::from_str("chr"), DataType::String);
                schema.insert(PlSmallStr::from_str("nuc"), DataType::String);
                schema.insert(PlSmallStr::from_str("position"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("context"), DataType::String);
                schema.insert(PlSmallStr::from_str("dinuc"), DataType::String);
                schema.insert(PlSmallStr::from_str("count_m"), DataType::UInt32);
                schema.insert(PlSmallStr::from_str("count_total"), DataType::UInt32);
                schema
            }
            ReportType::BEDGRAPH => {
                let mut schema = Schema::default();
                schema.insert(PlSmallStr::from_str("chr"), DataType::String);
                schema.insert(PlSmallStr::from_str("start"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("end"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("density"), DataType::Float32);
                schema
            }
        }
    }
    /// Get read options for report type
    fn read_options(
        &self,
        chunk_size: Option<usize>,
        low_memory: Option<bool>,
        n_threads: Option<usize>,
    ) -> CsvReadOptions {
        match self {
            ReportType::BISMARK
            | ReportType::COVERAGE
            | ReportType::CGMAP
            | ReportType::BEDGRAPH => CsvReadOptions::default()
                .with_has_header(false)
                .with_chunk_size(chunk_size.unwrap_or(2 << 20))
                .with_low_memory(low_memory.unwrap_or(false))
                .with_n_threads(n_threads)
                .with_schema(Some(SchemaRef::new(self.schema())))
                .with_parse_options({
                    CsvParseOptions::default()
                        .with_separator(b'\t')
                        .with_try_parse_dates(false)
                        .with_quote_char(Some(b'#'))
                }),
        }
    }

    /// Get mutation function to convert report format to BSX
    fn mutate_func(&self, data_frame: DataFrame, schema: &Schema) -> PolarsResult<DataFrame> {
        let lazy_frame = data_frame.lazy();
        let context_expr: Expr = when(col("context").eq(lit("CG")))
            .then(lit(true))
            .when(col("context").eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean);
        match self {
            ReportType::BEDGRAPH => todo!(),
            ReportType::COVERAGE => todo!(),
            ReportType::CGMAP => lazy_frame.with_columns([col("nuc").eq(lit("C")).alias("strand")]),
            ReportType::BISMARK => lazy_frame
                .with_column((col("count_m") + col("count_um")).alias("count_total"))
                .with_columns([
                    (col("count_m") / col("count_total")).alias("density"),
                    col("strand")
                        .eq(lit("+"))
                        .alias("strand"),
                ]),
        }
        .select([
            col("chr")
                .cast(schema.get_field("chr").unwrap().dtype)
                .cast(DataType::from_arrow(
                    &ArrowDataType::Dictionary(
                        IntegerType::UInt32,
                        Box::new(ArrowDataType::Utf8),
                        false,
                    ),
                    false,
                )),
            col("position").cast(
                schema
                    .get_field("position")
                    .unwrap()
                    .dtype,
            ),
            col("strand").cast(
                schema
                    .get_field("strand")
                    .unwrap()
                    .dtype,
            ),
            context_expr.alias("context").cast(
                schema
                    .get_field("context")
                    .unwrap()
                    .dtype,
            ),
            col("count_m").cast(
                schema
                    .get_field("count_m")
                    .unwrap()
                    .dtype,
            ),
            col("count_total").cast(
                schema
                    .get_field("count_total")
                    .unwrap()
                    .dtype,
            ),
            col("density").cast(
                schema
                    .get_field("density")
                    .unwrap()
                    .dtype,
            ),
        ])
        .collect()
    }
}

impl ReportType {
    /// Master function to convert report into BSX schema
    pub fn convert_report(
        &self,
        input_path: &str,
        output_path: &str,
        fasta_path: &str,
        chunk_size: Option<usize>,
        low_memory: Option<bool>,
        n_threads: Option<usize>,
    ) {
        // Init chromosome names
        let seq_data: SeqData = scan_fasta(fasta_path);
        let chr_names: Vec<&str> = seq_data
            .chr_names
            .iter()
            .map(|s| &s as &str)
            .collect::<Vec<&str>>();

        // Init reader
        let mut binding: CsvReader<File> = self
            .read_options(chunk_size, low_memory, n_threads)
            .try_into_reader_with_file_path(Some(input_path.parse().unwrap()))
            .expect("could not load CSV");
        let mut reader: BatchedCsvReader = binding
            .batched_borrowed()
            .expect("could not borrow CSV");

        // Init schema
        let schema: Schema = get_universal_schema(chr_names);

        // Init writer with categorical chromosome names
        let mut writer: BatchedWriter<File> = IpcWriterOptions::default()
            .to_writer(File::create(output_path).expect("could not create output file"))
            .batched(&schema)
            .expect("could not create ipc writer");

        // Read
        while let Some(batch_vec) = reader.next_batches(10).unwrap() {
            for batch in batch_vec {
                self.mutate_func(batch, &schema)
                    .expect("could not mutate record function")
                    .partition_by(["chr", "context", "strand"], true)
                    .expect("could not partition batch")
                    .iter_mut()
                    .try_for_each(|x| writer.write_batch(x.align_chunks()))
                    .expect("could not write batch");
            }
        }
        writer
            .finish()
            .expect("could not finish write");
    }
}
