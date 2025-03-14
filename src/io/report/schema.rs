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
use std::ops::Div;

use itertools::Itertools;
use log::{debug, trace, warn};
use polars::prelude::*;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3_polars::PySchema;

use crate::data_structs::bsx_batch::BsxBatch;
use crate::utils::{hashmap_from_arrays, schema_from_arrays};

/// Represents different methylation report file formats supported by the
/// application.
///
/// Each format has its own column structure and data_structs types that need to
/// be handled differently during import and transformation operations.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "python", pyclass)]
pub enum ReportTypeSchema {
    /// Bismark methylation extractor output format
    Bismark,
    /// CG methylation map format
    CgMap,
    /// BedGraph methylation density format
    BedGraph,
    /// Coverage report format with methylated/unmethylated counts
    Coverage,
}

impl ReportTypeSchema {
    /// Returns the column names for the specific report format.
    ///
    /// Each report type has a predefined set of columns that are expected
    /// in files of that format.
    pub const fn col_names(&self) -> &[&'static str] {
        match self {
            Self::Bismark => {
                &[
                    "chr", "position", "strand", "count_m", "count_um",
                    "context", "trinuc",
                ]
            },
            Self::Coverage => {
                &["chr", "start", "end", "density", "count_m", "count_um"]
            },
            Self::CgMap => {
                &[
                    "chr",
                    "nuc",
                    "position",
                    "context",
                    "dinuc",
                    "density",
                    "count_m",
                    "count_total",
                ]
            },
            Self::BedGraph => &["chr", "start", "end", "density"],
        }
    }

    /// Returns the data_structs types for each column in the report format.
    ///
    /// The order of data_structs types corresponds to the order of column names
    /// returned by `col_names()`.
    const fn col_types(&self) -> &[DataType] {
        match self {
            Self::Bismark => {
                &[
                    DataType::String, // chr
                    DataType::UInt64, // position
                    DataType::String, // strand
                    DataType::UInt32, // count_m
                    DataType::UInt32, // count_um
                    DataType::String, // context
                    DataType::String, // trinuc
                ]
            },
            Self::CgMap => {
                &[
                    DataType::String,  // chr
                    DataType::String,  // nuc
                    DataType::UInt64,  // position
                    DataType::String,  // context
                    DataType::String,  // dinuc
                    DataType::Float64, // density
                    DataType::UInt32,  // count_m
                    DataType::UInt32,  // count_total
                ]
            },
            Self::BedGraph => {
                &[
                    DataType::String,  // chr
                    DataType::UInt64,  // start
                    DataType::UInt64,  // end
                    DataType::Float64, // density
                ]
            },
            Self::Coverage => {
                &[
                    DataType::String,  // chr
                    DataType::UInt64,  // start
                    DataType::UInt64,  // end
                    DataType::Float64, // density
                    DataType::UInt32,  // count_m
                    DataType::UInt32,  // count_um
                ]
            },
        }
    }

    /// Returns the name of the chromosome column for this report type.
    pub const fn chr_col(&self) -> &'static str {
        match self {
            Self::Bismark | Self::BedGraph | Self::Coverage | Self::CgMap => {
                "chr"
            },
        }
    }

    /// Returns the name of the position column for this report type.
    ///
    /// Different report formats may use different column names to represent the
    /// genomic position.
    pub const fn position_col(&self) -> &'static str {
        match self {
            Self::Coverage | Self::BedGraph => "start",
            Self::Bismark | Self::CgMap => "position",
        }
    }

    /// Returns the name of the context column if it exists in this report type.
    ///
    /// Some formats might not include methylation context information.
    pub const fn context_col(&self) -> Option<&'static str> {
        match self {
            Self::Bismark | Self::CgMap => Some("context"),
            Self::BedGraph | Self::Coverage => None,
        }
    }

    /// Returns the name of the strand column if it exists in this report type.
    ///
    /// Some formats might not include strand information.
    pub const fn strand_col(&self) -> Option<&'static str> {
        match self {
            Self::Bismark | Self::CgMap => Some("strand"),
            Self::BedGraph | Self::Coverage => None,
        }
    }

    /// Indicates whether this report type needs alignment when processing.
    ///
    /// Some report formats contain range data_structs that requires special
    /// alignment handling.
    pub const fn need_align(&self) -> bool {
        match self {
            Self::Bismark | Self::CgMap => false,
            Self::BedGraph | Self::Coverage => true,
        }
    }

    /// Creates a Polars Schema for this report format.
    ///
    /// The schema defines the column names and their corresponding data_structs
    /// types.
    pub fn schema(&self) -> Schema {
        debug!("Creating schema for {:?} report format", self);
        schema_from_arrays(self.col_names(), self.col_types())
    }

    /// Creates a HashMap mapping column names to their data_structs types for
    /// this report format.
    pub fn hashmap(&self) -> PlHashMap<&str, DataType> {
        trace!("Creating column type HashMap for {:?} report format", self);
        hashmap_from_arrays(self.col_names(), self.col_types())
    }

    /// Creates CSV read options configured for this report format.
    ///
    /// The returned options are configured with appropriate settings for
    /// parsing files of this specific report type.
    pub fn read_options(&self) -> CsvReadOptions {
        debug!("Configuring CSV read options for {:?} format", self);
        let mut read_options = CsvReadOptions::default()
            .with_has_header(false) // None of the supported formats have headers
            .with_schema(Some(SchemaRef::from(self.schema())))
            .with_parse_options({
                CsvParseOptions::default()
                    .with_separator(b'\t')
                    .with_try_parse_dates(false)
                    .with_quote_char(Some(b'#'))
            });

        // BedGraph format has a header line that needs to be skipped
        if let Self::BedGraph = self {
            debug!("BedGraph format detected - configuring to skip header row");
            read_options = read_options.with_skip_rows(1);
        };

        read_options
    }

    /// Returns a function that converts a DataFrame from the current report
    /// format to BsxBatch format.
    ///
    /// The returned function transforms data_structs from the specific report
    /// format into the standard internal format (BsxBatch), handling column
    /// renames, type conversions, and calculated fields.
    ///
    /// # Returns
    /// A function that takes a DataFrame and returns a PolarsResult<DataFrame>
    /// in BsxBatch format. The returned DataFrame will contain all the
    /// columns required by BsxBatch, but some fields may be NULL or
    /// placeholder values depending on the source format.
    pub fn bsx_mutate(&self) -> fn(DataFrame) -> PolarsResult<DataFrame> {
        debug!("Creating bsx_mutate function for {:?} format", self);
        match self {
            ReportTypeSchema::Bismark => {
                |df: DataFrame| {
                    debug!("Converting Bismark format to BsxBatch format");
                    trace!(
                        "Input DataFrame shape: {}x{}",
                        df.height(),
                        df.width()
                    );

                    let result = df
                    .lazy()
                    // Calculate total count by adding methylated and unmethylated counts
                    .with_column((col("count_m") + col("count_um")).alias("count_total"))
                    // Calculate methylation density as the ratio of methylated to total count
                    .with_column(
                        col("count_m")
                            .cast(DataType::Float64)
                            .div(col("count_total"))
                            .alias("density"),
                    )
                    // Select only the columns needed for BsxBatch format
                    .select(BsxBatch::col_names().iter().map(|s| col(*s)).collect_vec())
                    // Ensure all columns have the correct data_structs types
                    .cast(BsxBatch::hashmap(), true)
                    .collect();

                    match &result {
                        Ok(df) => {
                            trace!(
                                "Transformed DataFrame shape: {}x{}",
                                df.height(),
                                df.width()
                            )
                        },
                        Err(e) => {
                            warn!("Error converting Bismark format: {}", e)
                        },
                    }

                    result
                }
            },
            ReportTypeSchema::CgMap => {
                |df: DataFrame| {
                    debug!("Converting CgMap format to BsxBatch format");
                    trace!(
                        "Input DataFrame shape: {}x{}",
                        df.height(),
                        df.width()
                    );

                    let result = df
                    .lazy()
                    // Determine strand based on nucleotide (C is +, G is -)
                    .with_columns([when(col("nuc").eq(lit("C")))
                        .then(lit("+"))
                        .when(col("nuc").eq(lit("G")))
                        .then(lit("-"))
                        .otherwise(lit("."))
                        .alias("strand")])
                    // Select only the columns needed for BsxBatch format
                    .select(BsxBatch::col_names().iter().map(|s| col(*s)).collect_vec())
                    // Ensure all columns have the correct data_structs types
                    .cast(BsxBatch::hashmap(), true)
                    .collect();

                    match &result {
                        Ok(df) => {
                            trace!(
                                "Transformed DataFrame shape: {}x{}",
                                df.height(),
                                df.width()
                            )
                        },
                        Err(e) => warn!("Error converting CgMap format: {}", e),
                    }

                    result
                }
            },
            ReportTypeSchema::BedGraph => {
                |df: DataFrame| {
                    debug!("Converting BedGraph format to BsxBatch format");
                    trace!(
                        "Input DataFrame shape: {}x{}",
                        df.height(),
                        df.width()
                    );

                    let result = df
                    .lazy()
                    // Add required columns with placeholder values
                    .with_columns([
                        col("start").alias("position"),
                        lit(".").alias("strand"),
                        lit(NULL).alias("context"),
                        lit(NULL).alias("count_m"),
                        lit(NULL).alias("count_total"),
                    ])
                    // Select only the columns needed for BsxBatch format
                    .select(BsxBatch::col_names().iter().map(|s| col(*s)).collect_vec())
                    // Ensure all columns have the correct data_structs types
                    .cast(BsxBatch::hashmap(), true)
                    .collect();

                    match &result {
                        Ok(df) => {
                            trace!(
                                "Transformed DataFrame shape: {}x{}",
                                df.height(),
                                df.width()
                            )
                        },
                        Err(e) => {
                            warn!("Error converting BedGraph format: {}", e)
                        },
                    }

                    result
                }
            },
            ReportTypeSchema::Coverage => {
                |df: DataFrame| {
                    debug!("Converting Coverage format to BsxBatch format");
                    trace!(
                        "Input DataFrame shape: {}x{}",
                        df.height(),
                        df.width()
                    );

                    let result = df
                    .lazy()
                    // Calculate total count by adding methylated and unmethylated counts
                    .with_column((col("count_m") + col("count_um")).alias("count_total"))
                    // Add required columns and rename existing ones
                    .with_columns([
                        col("start").alias("position"),
                        lit(".").alias("strand"),
                        lit(NULL).alias("context"),
                        // Recalculate density for consistency
                        col("count_m")
                            .cast(DataType::Float64)
                            .div(col("count_total"))
                            .alias("density"),
                    ])
                    // Select only the columns needed for BsxBatch format
                    .select(BsxBatch::col_names().iter().map(|s| col(*s)).collect_vec())
                    // Ensure all columns have the correct data_structs types
                    .cast(BsxBatch::hashmap(), true)
                    .collect();

                    match &result {
                        Ok(df) => {
                            trace!(
                                "Transformed DataFrame shape: {}x{}",
                                df.height(),
                                df.width()
                            )
                        },
                        Err(e) => {
                            warn!("Error converting Coverage format: {}", e)
                        },
                    }

                    result
                }
            },
        }
    }

    /// Converts a DataFrame from BsxBatch format back to the original report
    /// format.
    ///
    /// This is the inverse operation of `bsx_mutate()`, transforming from the
    /// internal representation back to the specific report format.
    ///
    /// # Parameters
    /// * `df` - A DataFrame in BsxBatch format
    ///
    /// # Returns
    /// A PolarsResult<DataFrame> in the original report format
    pub fn report_mutate_from_bsx(
        &self,
        df: DataFrame,
    ) -> PolarsResult<DataFrame> {
        debug!("Converting from BsxBatch to {:?} format", self);
        trace!("Input DataFrame shape: {}x{}", df.height(), df.width());

        let result = match self {
            ReportTypeSchema::Bismark => {
                debug!("Transforming to Bismark format");
                df.lazy()
                    // Calculate unmethylated count from total and methylated counts
                    .with_column((col("count_total") - col("count_m")).alias("count_um"))
                    // Copy context column to trinuc (they're the same in Bismark)
                    .with_column(col("context").alias("trinuc"))
                    // Select only the columns needed for Bismark format
                    .select(self.col_names().iter().map(|s| col(*s)).collect_vec())
                    // Ensure all columns have the correct data_structs types
                    .cast(self.hashmap(), true)
                    .collect()
            },
            ReportTypeSchema::CgMap => {
                debug!("Transforming to CgMap format");
                df.lazy()
                    // Determine nucleotide based on strand (+ is C, - is G)
                    .with_columns([when(col("strand").eq(lit("+")))
                        .then(lit("C"))
                        .when(col("strand").eq(lit("-")))
                        .then(lit("G"))
                        .otherwise(lit("."))
                        .alias("nuc")])
                    // Copy context column to dinuc (they're the same in CgMap)
                    .with_column(col("context").alias("dinuc"))
                    // Select only the columns needed for CgMap format
                    .select(self.col_names().iter().map(|s| col(*s)).collect_vec())
                    // Ensure all columns have the correct data_structs types
                    .cast(self.hashmap(), true)
                    .collect()
            },
            ReportTypeSchema::BedGraph => {
                debug!("Transforming to BedGraph format");
                df.lazy()
                    // Convert position to start/end columns
                    .with_columns([col("position").alias("start"), col("position").alias("end")])
                    // Remove any rows with NaN density values
                    .drop_nans(Some(vec![col("density")]))
                    // Select only the columns needed for BedGraph format
                    .select(self.col_names().iter().map(|s| col(*s)).collect_vec())
                    // Ensure all columns have the correct data_structs types
                    .cast(self.hashmap(), true)
                    .collect()
            },
            ReportTypeSchema::Coverage => {
                debug!("Transforming to Coverage format");
                df.lazy()
                    // Convert position to start/end columns
                    .with_columns([col("position").alias("start"), col("position").alias("end")])
                    // Remove any rows with NaN density values
                    .drop_nans(Some(vec![col("density")]))
                    // Calculate unmethylated count from total and methylated counts
                    .with_column((col("count_total") - col("count_m")).alias("count_um"))
                    // Select only the columns needed for Coverage format
                    .select(self.col_names().iter().map(|s| col(*s)).collect_vec())
                    // Ensure all columns have the correct data_structs types
                    .cast(self.hashmap(), true)
                    .collect()
            },
        };

        match &result {
            Ok(df) => {
                trace!(
                    "Transformed DataFrame shape: {}x{}",
                    df.height(),
                    df.width()
                )
            },
            Err(e) => warn!("Error converting from BsxBatch: {}", e),
        }

        result
    }
}

#[cfg(feature = "python")]
#[pymethods]
impl ReportTypeSchema {
    /// Gets the schema for this report type in a Python-compatible format.
    ///
    /// # Returns
    /// A PySchema object that can be used in Python code
    fn get_schema(&self) -> PyResult<PySchema> {
        debug!("Creating PySchema for {:?} report format", self);
        Ok(PySchema(SchemaRef::new(self.schema())))
    }
}

#[cfg(test)]
mod report_schema_test {
    use crate::data_structs::bsx_batch::BsxBatch;
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
        let reverse_transform = report_type
            .report_mutate_from_bsx(bsx_batch)
            .unwrap();
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
        let reverse_transform = report_type
            .report_mutate_from_bsx(bsx_batch)
            .unwrap();
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
        let reverse_transform = report_type
            .report_mutate_from_bsx(bsx_batch)
            .unwrap();
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
        let reverse_transform = report_type
            .report_mutate_from_bsx(bsx_batch)
            .unwrap();
        assert_eq!(reverse_transform, input_df);
    }
}
