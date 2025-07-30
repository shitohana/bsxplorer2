use std::convert::Infallible;
use std::fmt::Display;
use std::str::FromStr;

use polars::prelude::*;

use crate::utils::{
    hashmap_from_arrays,
    schema_from_arrays,
};

/// Supported methylation report file formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReportType {
    /// Bismark methylation extractor output format
    Bismark,
    /// CG methylation map format
    CgMap,
    /// BedGraph methylation density format
    BedGraph,
    /// Coverage report format with methylated/unmethylated counts
    Coverage,
}

impl FromStr for ReportType {
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "bismark" => Ok(ReportType::Bismark),
            "cgmap" => Ok(ReportType::CgMap),
            "bedgraph" => Ok(ReportType::BedGraph),
            "coverage" => Ok(ReportType::Coverage),
            _ => unimplemented!(),
        }
    }
}

impl Display for ReportType {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        let str = match self {
            ReportType::Bismark => String::from("bismark"),
            ReportType::CgMap => String::from("cgmap"),
            ReportType::BedGraph => String::from("bedgraph"),
            ReportType::Coverage => String::from("coverage"),
        };
        write!(f, "{}", str)
    }
}

impl ReportType {
    /// Returns column names for this report format.
    pub const fn col_names(&self) -> &[&'static str] {
        match self {
            Self::Bismark => {
                &[
                    "chr", "position", "strand", "count_m", "count_um", "context",
                    "trinuc",
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

    /// Returns data types for each column.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub const fn col_types(&self) -> &[DataType] {
        match self {
            Self::Bismark => {
                &[
                    DataType::String, // chr
                    DataType::UInt32, // position
                    DataType::String, // strand
                    DataType::UInt16, // count_m
                    DataType::UInt16, // count_um
                    DataType::String, // context
                    DataType::String, // trinuc
                ]
            },
            Self::CgMap => {
                &[
                    DataType::String,  // chr
                    DataType::String,  // nuc
                    DataType::UInt32,  // position
                    DataType::String,  // context
                    DataType::String,  // dinuc
                    DataType::Float32, // density
                    DataType::UInt16,  // count_m
                    DataType::UInt16,  // count_total
                ]
            },
            Self::BedGraph => {
                &[
                    DataType::String,  // chr
                    DataType::UInt32,  // start
                    DataType::UInt32,  // end
                    DataType::Float32, // density
                ]
            },
            Self::Coverage => {
                &[
                    DataType::String,  // chr
                    DataType::UInt32,  // start
                    DataType::UInt32,  // end
                    DataType::Float32, // density
                    DataType::UInt16,  // count_m
                    DataType::UInt16,  // count_um
                ]
            },
        }
    }

    /// Returns chromosome column name.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub const fn chr_col(&self) -> &'static str {
        match self {
            Self::Bismark | Self::BedGraph | Self::Coverage | Self::CgMap => "chr",
        }
    }

    /// Returns position column name.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub const fn position_col(&self) -> &'static str {
        match self {
            Self::Coverage | Self::BedGraph => "start",
            Self::Bismark | Self::CgMap => "position",
        }
    }

    /// Returns context column name if available.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub const fn context_col(&self) -> Option<&'static str> {
        match self {
            Self::Bismark | Self::CgMap => Some("context"),
            Self::BedGraph | Self::Coverage => None,
        }
    }

    /// Returns strand column name if available.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub const fn strand_col(&self) -> Option<&'static str> {
        match self {
            Self::Bismark | Self::CgMap => Some("strand"),
            Self::BedGraph | Self::Coverage => None,
        }
    }

    /// Whether this format needs alignment when processing.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub const fn need_align(&self) -> bool {
        match self {
            Self::Bismark | Self::CgMap => false,
            Self::BedGraph | Self::Coverage => true,
        }
    }

    /// Creates a Polars Schema for this format.
    pub fn schema(&self) -> Schema {
        schema_from_arrays(self.col_names(), self.col_types())
    }

    /// Creates a HashMap of column names to data types.
    pub fn hashmap(&self) -> PlHashMap<&str, DataType> {
        hashmap_from_arrays(self.col_names(), self.col_types())
    }

    /// Creates CSV read options for this format.
    pub fn read_options(&self) -> CsvReadOptions {
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
            read_options = read_options.with_skip_rows(1);
        };

        read_options
    }
}
