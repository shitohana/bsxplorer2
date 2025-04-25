use std::error::Error;

use polars::prelude::*;
use serde::Deserialize;
use crate::utils::{hashmap_from_arrays, schema_from_arrays};

/// Supported methylation report file formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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
    /// Returns column names for this report format.
    pub const fn col_names(&self) -> &[&'static str] {
        match self {
            Self::Bismark => &[
                "chr", "position", "strand", "count_m", "count_um", "context",
                "trinuc",
            ],
            Self::Coverage => {
                &["chr", "start", "end", "density", "count_m", "count_um"]
            },
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

    /// Returns data types for each column.
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

    /// Returns chromosome column name.
    pub const fn chr_col(&self) -> &'static str {
        match self {
            Self::Bismark | Self::BedGraph | Self::Coverage | Self::CgMap => {
                "chr"
            },
        }
    }

    /// Returns position column name.
    pub const fn position_col(&self) -> &'static str {
        match self {
            Self::Coverage | Self::BedGraph => "start",
            Self::Bismark | Self::CgMap => "position",
        }
    }

    /// Returns context column name if available.
    pub const fn context_col(&self) -> Option<&'static str> {
        match self {
            Self::Bismark | Self::CgMap => Some("context"),
            Self::BedGraph | Self::Coverage => None,
        }
    }

    /// Returns strand column name if available.
    pub const fn strand_col(&self) -> Option<&'static str> {
        match self {
            Self::Bismark | Self::CgMap => Some("strand"),
            Self::BedGraph | Self::Coverage => None,
        }
    }

    /// Whether this format needs alignment when processing.
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

    /// Returns a `csv::StringRecord` with the column names for this format.
    pub(crate) fn header_record(&self) -> csv::StringRecord {
        let mut record = csv::StringRecord::new();
        for name in self.col_names() {
            record.push_field(name);
        }
        record
    }
}

pub(crate) trait ReportRow<'a>: Deserialize<'a> {
    fn get_chr(&self) -> String;
    fn get_pos(&self) -> usize;
}

#[derive(Deserialize, Debug)]
pub(crate) struct BismarkRow {
    chr: String,
    position: usize,
    strand: String,
    count_m: u32,
    count_um: u32,
    context: String,
    trinuc: String,
}

impl<'a> ReportRow<'a> for BismarkRow {
    fn get_chr(&self) -> String {
        self.chr.clone()
    }
    fn get_pos(&self) -> usize {
        self.position
    }
}

#[derive(Deserialize, Debug)]
pub(crate) struct CgMapRow {
    chr: String,
    nuc: String,
    position: usize,
    context: String,
    dinuc: String,
    density: f64,
    count_m: u32,
    count_total: u32,
}

impl<'a> ReportRow<'a> for CgMapRow {
    fn get_chr(&self) -> String {
        self.chr.clone()
    }
    fn get_pos(&self) -> usize {
        self.position
    }
}

#[derive(Deserialize, Debug)]
pub(crate) struct BedGraphRow {
    chr: String,
    start: usize,
    end: usize,
    density: f64,
}

impl<'a> ReportRow<'a> for BedGraphRow {
    fn get_chr(&self) -> String {
        self.chr.clone()
    }
    fn get_pos(&self) -> usize {
        self.start
    }
}

#[derive(Deserialize, Debug)]
pub(crate) struct CoverageRow {
    chr: String,
    start: usize,
    end: usize,
    density: f64,
    count_m: u32,
    count_um: u32,
}

impl<'a> ReportRow<'a> for CoverageRow {
    fn get_chr(&self) -> String {
        self.chr.clone()
    }
    fn get_pos(&self) -> usize {
        self.start
    }
}