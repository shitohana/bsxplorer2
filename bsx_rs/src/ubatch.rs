use crate::io::bsx::read::ReadFilters;
use crate::io::report::types::ReportType;
use crate::region::{RegionCoordinates, RegionData};
use crate::utils::types::{Context, Strand};
use itertools::Itertools;
use polars::prelude::*;
use std::fmt::Display;
use std::ops::Not;

struct UniversalBatch {
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
        assert_eq!(
            Self::colnames()
                .iter()
                .map(|&name| name.to_string())
                .collect_vec(),
            data.schema()
                .iter_names()
                .map(|val| val.to_string())
                .collect_vec(),
            "DataFrame columns does not match schema"
        );
        // Schema consistency check
        assert!(
            data_schema.get("chr").unwrap().is_categorical(),
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
                "DataFrame column 'density' must be float"
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

        data = data
            .select_with_schema(
                data_schema.iter_names().map(|x| x.clone()),
                &SchemaRef::from(data_schema.clone()),
            )
            .unwrap();

        Self { data }
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
            .categorical()
            .unwrap()
            .physical()
            .n_unique()
            .unwrap()
            == 1
    }
    /// Filter self by [ReadFilters]
    pub fn filter(mut self, filter: ReadFilters) -> Self {
        if let Some(context) = filter.context {
            self = self.filter_context(context)
        }
        if let Some(strand) = filter.strand {
            self = self.filter_strand(strand)
        }
        self
    }
    /// Create new [UniversalBatch] by [RegionCoordinates] slice.
    pub fn slice(&mut self, coordinates: RegionCoordinates) -> RegionData {
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
        RegionData::new(new_data, coordinates)
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
