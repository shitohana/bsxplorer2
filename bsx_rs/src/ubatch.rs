use std::ops::Div;
use std::cmp::Ordering;
use crate::io::bsx::region_reader::ReadFilters;
use crate::io::report::types::ReportType;
use crate::region::{RegionCoordinates, RegionData};
use crate::utils::types::{Context, Strand};
use itertools::Itertools;
use polars::prelude::*;
use std::fmt::Display;
use std::ops::Not;
use polars::export::num::WrappingNeg;
use crate::io::report::read::BatchStats;

#[derive(Debug, Clone)]
pub struct UniversalBatch {
    data: DataFrame,
}

impl UniversalBatch {
    fn colnames() -> Vec<&'static str> {
        vec![
            "chr",
            "strand",
            "position",
            "context",
            "count_m",
            "count_total",
            "density",
        ]
    }

    pub fn new(mut data: DataFrame) -> Self {
        let mut data_schema = data.schema();
        
        assert!(
            Self::colnames().into_iter().all(|name| data_schema.get(name).is_some()),
            "DataFrame columns does not match schema"
        );
        
        // Schema consistency check
        assert!(
            data_schema.get("chr").unwrap().is_string(),
            "DataFrame column 'chr' must be categorical"
        );
        assert!(
            data_schema.get("strand").unwrap().is_bool(),
            "DataFrame column 'strand' must be bool"
        );
        assert!(
            data_schema.get("context").unwrap().is_bool(),
            "DataFrame column 'strand' must be bool"
        );

        for (col, dtype) in [
            ("position", DataType::UInt32),
            ("count_m", DataType::UInt16),
            ("count_total", DataType::UInt16),
        ] {
            let orig_dtype = data_schema.get(col).unwrap().clone();
            assert!(
                orig_dtype.is_unsigned_integer(),
                "DataFrame column '{col}' must be unsigned integer"
            );
            if orig_dtype != dtype {
                data_schema.set_dtype(col, dtype);
            }
        }
        {
            let orig_dtype = data_schema.get("density").unwrap().clone();
            assert!(
                orig_dtype.is_float(),
                "DataFrame column 'density' must be float {orig_dtype}"
            );
            if orig_dtype != DataType::Float64 {
                data_schema.set_dtype("density", DataType::Float64);
            }
        }

        // Null check
        for col in ["chr", "position"] {
            assert_eq!(
                data.column(col).unwrap().null_count(),
                0,
                "DataFrame column '{col}' cannot be null"
            );
        }
        let data = data.lazy()
            .cast(
                PlHashMap::from_iter(
                    data_schema.iter()
                        .map(|(x, _dtype)| (x.as_str(), _dtype.clone()))
                ), 
                true
            ).collect().unwrap();

        Self { data }
    }
    
    pub fn from_report_type(data: DataFrame, report_type: ReportType) -> Result<Self, PolarsError> {
        let lazy_frame = data.lazy();

        let context_encoder = when(col("context").eq(lit("CG")))
            .then(lit(true))
            .when(col("context").eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean);
        
        let res = match report_type {
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
        }.collect();
        
        Ok(Self::new(res?))
    }

    pub fn into_report_type(self, report_type: ReportType) -> DataFrame {
        let target_schema = report_type.get_schema();
        let mut data_lazy = self.data.lazy();
        if target_schema.get("strand").is_some() {
            data_lazy = data_lazy.with_column(
                when(col("strand").eq(lit(false)))
                    .then(lit("-"))
                    .when(col("strand").eq(lit(true)))
                    .then(lit("+"))
                    .otherwise(lit(".")),
            )
        }
        if target_schema.get("context").is_some() {
            data_lazy = data_lazy.with_column(
                when(col("context").is_null())
                    .then(lit("CHH"))
                    .when(col("context").eq(lit(false)))
                    .then(lit("CHG"))
                    .otherwise(lit("CG")),
            )
        }
        let new_data = data_lazy
            .cast(
                PlHashMap::from_iter(
                    target_schema
                        .iter_names_and_dtypes()
                        .map(|(name, dtype)| (name.as_str(), dtype.clone())),
                ),
                true,
            )
            .collect()
            .unwrap();
        new_data
    }

    /// Get reference to inner DataFrame
    pub fn get_data(&self) -> &DataFrame {
        &self.data
    }
    /// Get mutable reference to inner DataFrame
    pub fn get_data_mut(&mut self) -> &mut DataFrame {
        &mut self.data
    }
    /// Get first position in batch
    pub fn first_position(&self) -> u32 {
        self.data
            .column("position")
            .unwrap()
            .u32()
            .unwrap()
            .first()
            .unwrap()
    }
    /// Get last position in batch
    pub fn last_position(&self) -> u32 {
        self.data
            .column("position")
            .unwrap()
            .u32()
            .unwrap()
            .last()
            .unwrap()
    }
    /// Check if chromosome categorical indexes are unique
    pub fn unique_chr(&self) -> bool {
        self.data
            .column("chr")
            .unwrap()
            .str()
            .unwrap()
            .n_unique()
            .unwrap()
            == 1
    }
    pub fn get_chr(&self) -> String {
        self.data.column("chr").unwrap().str().unwrap().first().unwrap().to_string()
    }
    
    /// Filter self by [ReadFilters]
    pub fn filter(mut self, filter: &ReadFilters) -> Self {
        if let Some(context) = filter.context {
            self = self.filter_context(context)
        }
        if let Some(strand) = filter.strand {
            self = self.filter_strand(strand)
        }
        self
    }
    /// Create new [UniversalBatch] by [RegionCoordinates] slice.
    pub fn slice(&self, coordinates: &RegionCoordinates) -> RegionData {
        let positions = self.data.column("position").unwrap().u32().unwrap();
        let start = positions
            .iter()
            .position(|x| coordinates.start < x.unwrap())
            .unwrap_or_else(|| 0);
        let slice_length = positions
            .iter()
            .skip(start)
            .position(|x| coordinates.end < x.unwrap())
            .unwrap_or(positions.len() - start - 1);

        let new_data = self.data.slice(start as i64, slice_length);
        RegionData::new(new_data, coordinates.clone())
    }
    /// Filter self by [Context]
    pub fn filter_context(mut self, context: Context) -> Self {
        let condition = match context {
            Context::CG => self
                .data
                .column("context")
                .unwrap()
                .bool()
                .unwrap()
                .to_owned(),
            Context::CHG => self.data.column("context").unwrap().bool().unwrap().not(),
            Context::CHH => self
                .data
                .column("context")
                .unwrap()
                .bool()
                .unwrap()
                .is_null(),
            Context::ALL => return self,
        };
        self.data = self.data.filter(&condition).unwrap();
        self
    }
    /// Filter self by [Strand]
    pub fn filter_strand(mut self, strand: Strand) -> Self {
        let condition = match strand {
            Strand::Forward => self
                .data
                .column("strand")
                .unwrap()
                .bool()
                .unwrap()
                .to_owned(),
            Strand::Reverse => self.data.column("strand").unwrap().bool().unwrap().not(),
            _ => return self,
        };
        self.data = self.data.filter(&condition).unwrap();
        self
    }
}

impl Display for UniversalBatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.data)
    }
}

impl From<DataFrame> for UniversalBatch {
    fn from(data_frame: DataFrame) -> Self {
        UniversalBatch::new(data_frame)
    }
}

impl From<UniversalBatch> for DataFrame {
    fn from(data_frame: UniversalBatch) -> Self {
        data_frame.data
    }
}

impl Eq for UniversalBatch {}

impl PartialEq<Self> for UniversalBatch {
    fn eq(&self, other: &Self) -> bool {
        self.first_position() == other.first_position()
    }
}

impl PartialOrd<Self> for UniversalBatch {
    /// Order implemented as reversed for correct ordering in Max-heap
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        ((self.first_position() as i32).wrapping_neg()).partial_cmp(&(other.first_position() as i32).wrapping_neg())
    }
}

impl Ord for UniversalBatch {
    /// Order implemented as reversed for correct ordering in Max-heap
    fn cmp(&self, other: &Self) -> Ordering {
        ((self.first_position() as i32).wrapping_neg()).cmp(&(other.first_position() as i32).wrapping_neg())
    }
}