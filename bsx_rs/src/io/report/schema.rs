use crate::bsx_batch::BsxBatch;
use crate::utils::{hashmap_from_arrays, schema_from_arrays};
use itertools::Itertools;
use polars::prelude::*;
use std::ops::Div;

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3_polars::PySchema;

#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "python", pyclass)]
pub enum ReportTypeSchema {
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
}

impl ReportTypeSchema {
    pub const fn col_names(&self) -> &[&'static str] {
        match self {
            Self::Bismark => &[
                "chr", "position", "strand", "count_m", "count_um", "context", "trinuc",
            ],
            Self::Coverage => &["chr", "start", "end", "density", "count_m", "count_um"],
            Self::CgMap => &[
                "chr",
                "nuc",
                "position",
                "context",
                "dinuc",
                "density",
                "count_m",
                "count_total",
            ],
            Self::BedGraph => &["chr", "start", "end", "density"],
        }
    }

    const fn col_types(&self) -> &[DataType] {
        match self {
            Self::Bismark => &[
                DataType::String,
                DataType::UInt64,
                DataType::String,
                DataType::UInt32,
                DataType::UInt32,
                DataType::String,
                DataType::String,
            ],
            Self::CgMap => &[
                DataType::String,
                DataType::String,
                DataType::UInt64,
                DataType::String,
                DataType::String,
                DataType::Float64,
                DataType::UInt32,
                DataType::UInt32,
            ],
            Self::BedGraph => &[
                DataType::String,
                DataType::UInt64,
                DataType::UInt64,
                DataType::Float64,
            ],
            Self::Coverage => &[
                DataType::String,
                DataType::UInt64,
                DataType::UInt64,
                DataType::Float64,
                DataType::UInt32,
                DataType::UInt32,
            ],
        }
    }

    pub const fn chr_col(&self) -> &'static str {
        match self {
            Self::Bismark | Self::BedGraph | Self::Coverage | Self::CgMap => "chr",
        }
    }

    pub const fn position_col(&self) -> &'static str {
        match self {
            Self::Coverage | Self::BedGraph => "start",

            Self::Bismark | Self::CgMap => "position",
        }
    }

    pub const fn context_col(&self) -> Option<&'static str> {
        match self {
            Self::Bismark | Self::CgMap => Some("context"),

            Self::BedGraph | Self::Coverage => None,
        }
    }

    pub const fn strand_col(&self) -> Option<&'static str> {
        match self {
            Self::Bismark | Self::CgMap => Some("strand"),

            Self::BedGraph | Self::Coverage => None,
        }
    }

    pub const fn need_align(&self) -> bool {
        match self {
            Self::Bismark | Self::CgMap => false,

            Self::BedGraph | Self::Coverage => true,
        }
    }

    pub fn schema(&self) -> Schema {
        schema_from_arrays(self.col_names(), self.col_types())
    }

    pub fn hashmap(&self) -> PlHashMap<&str, DataType> {
        hashmap_from_arrays(self.col_names(), self.col_types())
    }

    pub fn read_options(&self) -> CsvReadOptions {
        let mut read_options = CsvReadOptions::default()
            // As no yet supported formats have headers
            .with_has_header(false)
            .with_schema(Some(SchemaRef::from(self.schema())))
            .with_parse_options({
                CsvParseOptions::default()
                    .with_separator(b'\t')
                    .with_try_parse_dates(false)
                    .with_quote_char(Some(b'#'))
            });
        if let Self::BedGraph = self {
            read_options = read_options.with_skip_rows(1);
        };
        read_options
    }

    /// Returns DataFrame with columns and types as in [BsxBatch],
    /// but context, strand and counts data are possibly invalid or non
    /// present.
    pub fn bsx_mutate(&self) -> fn(DataFrame) -> PolarsResult<DataFrame> {
        match self {
            ReportTypeSchema::Bismark => |df: DataFrame| {
                df.lazy()
                    .with_column((col("count_m") + col("count_um")).alias("count_total"))
                    .with_column(
                        (col("count_m")
                            .cast(DataType::Float64)
                            .div(col("count_total")))
                        .alias("density"),
                    )
                    .select(BsxBatch::col_names().iter().map(|s| col(*s)).collect_vec())
                    .cast(BsxBatch::hashmap(), true)
                    .collect()
            },
            ReportTypeSchema::CgMap => |df: DataFrame| {
                df.lazy()
                    .with_columns([when(col("nuc").eq(lit("C")))
                        .then(lit("+"))
                        .when(col("nuc").eq(lit("G")))
                        .then(lit("-"))
                        .otherwise(lit("."))
                        .alias("strand")])
                    .select(BsxBatch::col_names().iter().map(|s| col(*s)).collect_vec())
                    .cast(BsxBatch::hashmap(), true)
                    .collect()
            },
            ReportTypeSchema::BedGraph => |df: DataFrame| {
                df.lazy()
                    .with_columns([
                        col("start").alias("position"),
                        lit(".").alias("strand"),
                        lit(NULL).alias("context"),
                        lit(NULL).alias("count_m"),
                        lit(NULL).alias("count_total"),
                    ])
                    .select(BsxBatch::col_names().iter().map(|s| col(*s)).collect_vec())
                    .cast(BsxBatch::hashmap(), true)
                    .collect()
            },
            ReportTypeSchema::Coverage => |df: DataFrame| {
                df.lazy()
                    .with_column((col("count_m") + col("count_um")).alias("count_total"))
                    .with_columns([
                        col("start").alias("position"),
                        lit(".").alias("strand"),
                        lit(NULL).alias("context"),
                        (col("count_m")
                            .cast(DataType::Float64)
                            .div(col("count_total")))
                        .alias("density"),
                    ])
                    .select(BsxBatch::col_names().iter().map(|s| col(*s)).collect_vec())
                    .cast(BsxBatch::hashmap(), true)
                    .collect()
            },
        }
    }

    pub fn report_mutate_from_bsx(&self, df: DataFrame) -> PolarsResult<DataFrame> {
        match self {
            ReportTypeSchema::Bismark => df
                .lazy()
                .with_column((col("count_total") - col("count_m")).alias("count_um"))
                .with_column(col("context").alias("trinuc"))
                .select(self.col_names().iter().map(|s| col(*s)).collect_vec())
                .cast(self.hashmap(), true)
                .collect(),
            ReportTypeSchema::CgMap => df
                .lazy()
                .with_columns([when(col("strand").eq(lit("+")))
                    .then(lit("C"))
                    .when(col("strand").eq(lit("-")))
                    .then(lit("G"))
                    .otherwise(lit("."))
                    .alias("nuc")])
                .with_column(col("context").alias("dinuc"))
                .select(self.col_names().iter().map(|s| col(*s)).collect_vec())
                .cast(self.hashmap(), true)
                .collect(),
            ReportTypeSchema::BedGraph => df
                .lazy()
                .with_columns([col("position").alias("start"), col("position").alias("end")])
                .drop_nans(Some(vec![col("density")]))
                .select(self.col_names().iter().map(|s| col(*s)).collect_vec())
                .cast(self.hashmap(), true)
                .collect(),
            ReportTypeSchema::Coverage => df
                .lazy()
                .with_columns([col("position").alias("start"), col("position").alias("end")])
                .drop_nans(Some(vec![col("density")]))
                .with_column((col("count_total") - col("count_m")).alias("count_um"))
                .select(self.col_names().iter().map(|s| col(*s)).collect_vec())
                .cast(self.hashmap(), true)
                .collect(),
        }
    }
}

#[cfg(feature = "python")]
#[pymethods]
impl ReportTypeSchema {
    fn get_schema(&self) -> PyResult<PySchema> {
        Ok(PySchema(SchemaRef::new(self.schema())))
    }
}

#[cfg(test)]
mod report_schema_test {
    use crate::bsx_batch::BsxBatch;
    use crate::io::report::schema::ReportTypeSchema;
    use crate::io::report::*;

    #[test]
    fn bismark_conversion() {
        let report_type = ReportTypeSchema::Bismark;

        let input_df = df![
            "chr" => ["1", "2", "3"],
            "position" => [1, 2, 3],
            "strand" => [".", ".", "."],
            "context" => ["CG", "CHG", "CHH"],
            "count_m" => [0, 1, 0],
            "count_um" => [1, 0, 1],
            "trinuc" => ["CG", "CHG", "CHH"]
        ]
        .unwrap()
        .lazy()
        .cast(report_type.hashmap(), true)
        .collect()
        .unwrap()
        .select(report_type.col_names().iter().cloned())
        .unwrap();

        let output_df = df![
            "chr" => ["1", "2", "3"],
            "position" => [1, 2, 3],
            "strand" => [".", ".", "."],
            "context" => ["CG", "CHG", "CHH"],
            "count_m" => [0, 1, 0],
            "count_total" => [1, 1, 1],
            "density" => [0, 1, 0]
        ]
        .unwrap()
        .lazy()
        .cast(BsxBatch::hashmap(), true)
        .collect()
        .unwrap();

        let mutate_func = report_type.bsx_mutate();
        let bsx_batch = mutate_func(input_df.clone()).unwrap();
        assert_eq!(bsx_batch, output_df);
        let reverse_transform = report_type.report_mutate_from_bsx(bsx_batch).unwrap();
        assert_eq!(reverse_transform, input_df);
    }

    #[test]
    fn cgmap_conversion() {
        let report_type = ReportTypeSchema::CgMap;

        let input_df = df![
            "chr" => ["1", "2", "3"],
            "nuc" => ["G", "C", "C"],
            "position" => [1, 2, 3],
            "context" => ["CG", "CHG", "CHH"],
            "dinuc" => ["CG", "CHG", "CHH"],
            "density" => [0, 1, 0],
            "count_m" => [0, 1, 0],
            "count_total" => [1, 1, 1],
        ]
        .unwrap()
        .lazy()
        .cast(report_type.hashmap(), true)
        .collect()
        .unwrap();

        let output_df = df![
            "chr" => ["1", "2", "3"],
            "position" => [1, 2, 3],
            "strand" => ["-", "+", "+"],
            "context" => ["CG", "CHG", "CHH"],
            "count_m" => [0, 1, 0],
            "count_total" => [1, 1, 1],
            "density" => [0, 1, 0]
        ]
        .unwrap()
        .lazy()
        .cast(BsxBatch::hashmap(), true)
        .collect()
        .unwrap();

        let mutate_func = report_type.bsx_mutate();
        let bsx_batch = mutate_func(input_df.clone()).unwrap();
        assert_eq!(bsx_batch, output_df);
        let reverse_transform = report_type.report_mutate_from_bsx(bsx_batch).unwrap();
        assert_eq!(reverse_transform, input_df);
    }

    #[test]
    fn coverage_conversion() {
        let report_type = ReportTypeSchema::Coverage;
        let input_df = df![
            "chr" => ["1", "2", "3"],
            "start" => [1, 2, 3],
            "end" => [1, 2, 3],
            "density" => [0, 1, 0],
            "count_m" => [0, 1, 0],
            "count_um" => [1, 0, 1],
        ]
        .unwrap()
        .lazy()
        .cast(report_type.hashmap(), true)
        .collect()
        .unwrap();

        let output_df = df![
            "chr" => ["1", "2", "3"],
            "position" => [1, 2, 3],
            "strand" => [".", ".", "."],
            "context" => [None::<&str>, None, None],
            "count_m" => [0, 1, 0],
            "count_total" => [1, 1, 1],
            "density" => [0, 1, 0]
        ]
        .unwrap()
        .lazy()
        .cast(BsxBatch::hashmap(), true)
        .collect()
        .unwrap();

        let mutate_func = report_type.bsx_mutate();
        let bsx_batch = mutate_func(input_df.clone()).unwrap();
        assert_eq!(bsx_batch, output_df);
        let reverse_transform = report_type.report_mutate_from_bsx(bsx_batch).unwrap();
        assert_eq!(reverse_transform, input_df);
    }

    #[test]
    fn bedgraph_conversion() {
        let report_type = ReportTypeSchema::BedGraph;

        let input_df = df![
            "chr" => ["1", "2", "3"],
            "start" => [1, 2, 3],
            "end" => [1, 2, 3],
            "density" => [0, 1, 0],
        ]
        .unwrap()
        .lazy()
        .cast(report_type.hashmap(), true)
        .collect()
        .unwrap();

        let output_df = df![
            "chr" => ["1", "2", "3"],
            "position" => [1, 2, 3],
            "strand" => [".", ".", "."],
            "context" => [None::<&str>, None, None],
            "count_m" => [None::<u32>, None, None],
            "count_total" => [None::<u32>, None, None],
            "density" => [0, 1, 0]
        ]
        .unwrap()
        .lazy()
        .cast(BsxBatch::hashmap(), true)
        .collect()
        .unwrap();

        let mutate_func = report_type.bsx_mutate();
        let bsx_batch = mutate_func(input_df.clone()).unwrap();
        assert_eq!(bsx_batch, output_df);
        let reverse_transform = report_type.report_mutate_from_bsx(bsx_batch).unwrap();
        assert_eq!(reverse_transform, input_df);
    }
}
