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

/* TODO:    Create batch cache as a continous universal batch, maybe sorted
            maybe split by chromosome, which is extended using binary heap 
            sorting of incoming batches.
            Then .next() method of ReportReader should only split such batch cache 
            table at specified index. No additional split to chunk size is needed. 
*/

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum UBatchCol {
    Chr,
    Strand,
    Position,
    Context,
    CountM,
    CountUm,
    Density
}

impl UBatchCol {
    fn as_str(&self) -> &'static str {
        match self { 
            Self::Chr => "chr",
            Self::Strand => "strand",
            Self::Position => "position",
            Self::Context => "context",
            Self::CountM => "count_m",
            Self::CountUm => "count_um",
            Self::Density => "density"
        }
    }
    
    fn data_type(&self) -> DataType {
        match self { 
            Self::Chr => DataType::Categorical(None, CategoricalOrdering::Physical),
            Self::Strand => DataType::Boolean,
            Self::Position => DataType::UInt64,
            Self::Context => DataType::Boolean,
            Self::CountM => DataType::UInt16,
            Self::CountUm => DataType::UInt16,
            Self::Density => DataType::Float64,
        }
    }
    
    fn all_cols() -> Vec<String> {
        vec![
            Self::Chr.as_str().to_string(),
            Self::Strand.as_str().to_string(),
            Self::Position.as_str().to_string(),
            Self::Context.as_str().to_string(),
            Self::CountM.as_str().to_string(),
            Self::CountUm.as_str().to_string(),
            Self::Density.as_str().to_string(),
        ]
    }
    
    fn schema() -> Schema {
        Schema::from_iter([
            (Self::Chr.as_str().into(), Self::Chr.data_type()),
            (Self::Strand.as_str().into(), Self::Strand.data_type()),
            (Self::Position.as_str().into(), Self::Position.data_type()),
            (Self::Context.as_str().into(), Self::Context.data_type()),
            (Self::CountM.as_str().into(), Self::CountM.data_type()),
            (Self::CountUm.as_str().into(), Self::CountUm.data_type()),
            (Self::Density.as_str().into(), Self::Density.data_type()),
        ])
    }
    
    fn can_cast(&self, data_type: DataType) -> bool {
        match self { 
            Self::Chr => data_type.is_categorical(),
            Self::Strand => data_type.is_bool(),
            Self::Position => data_type.is_unsigned_integer(),
            Self::Context => data_type.is_bool(),
            Self::CountM => data_type.is_unsigned_integer(),
            Self::CountUm => data_type.is_unsigned_integer(),
            Self::Density => data_type.is_float(),
        }
    }
    
    fn from_string(string: String) -> Self {
        match string.as_str() {
            "chr" => Self::Chr,
            "strand" => Self::Strand,
            "position" => Self::Position,
            "context" => Self::Context,
            "count_m" => Self::CountM,
            "count_um" => Self::CountUm,
            "density" => Self::Density,
            _ => panic!("Unknown BSX data type {}", string)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_str() {
        assert_eq!(UBatchCol::Chr.as_str(), "chr");
        assert_eq!(UBatchCol::Strand.as_str(), "strand");
        assert_eq!(UBatchCol::Position.as_str(), "position");
        assert_eq!(UBatchCol::Context.as_str(), "context");
        assert_eq!(UBatchCol::CountM.as_str(), "count_m");
        assert_eq!(UBatchCol::CountUm.as_str(), "count_um");
        assert_eq!(UBatchCol::Density.as_str(), "density");
    }

    #[test]
    fn test_data_type() {
        assert_eq!(UBatchCol::Chr.data_type(), DataType::Categorical(None, CategoricalOrdering::Physical));
        assert_eq!(UBatchCol::Strand.data_type(), DataType::Boolean);
        assert_eq!(UBatchCol::Position.data_type(), DataType::UInt64);
        assert_eq!(UBatchCol::Context.data_type(), DataType::Boolean);
        assert_eq!(UBatchCol::CountM.data_type(), DataType::UInt16);
        assert_eq!(UBatchCol::CountUm.data_type(), DataType::UInt16);
        assert_eq!(UBatchCol::Density.data_type(), DataType::Float64);
    }

    #[test]
    fn test_all_cols() {
        let expected = vec![
            "chr".to_string(),
            "strand".to_string(),
            "position".to_string(),
            "context".to_string(),
            "count_m".to_string(),
            "count_um".to_string(),
            "density".to_string(),
        ];
        assert_eq!(UBatchCol::all_cols(), expected);
    }

    #[test]
    fn test_schema() {
        let schema = UBatchCol::schema();
        let expected_schema = Schema::from_iter([
            ("chr".into(), DataType::Categorical(None, CategoricalOrdering::Physical)),
            ("strand".into(), DataType::Boolean),
            ("position".into(), DataType::UInt64),
            ("context".into(), DataType::Boolean),
            ("count_m".into(), DataType::UInt16),
            ("count_um".into(), DataType::UInt16),
            ("density".into(), DataType::Float64),
        ]);
        assert_eq!(schema, expected_schema);
    }

    #[test]
    fn test_from_string() {
        assert_eq!(UBatchCol::from_string("chr".to_string()), UBatchCol::Chr);
        assert_eq!(UBatchCol::from_string("strand".to_string()), UBatchCol::Strand);
        assert_eq!(UBatchCol::from_string("position".to_string()), UBatchCol::Position);
        assert_eq!(UBatchCol::from_string("context".to_string()), UBatchCol::Context);
        assert_eq!(UBatchCol::from_string("count_m".to_string()), UBatchCol::CountM);
        assert_eq!(UBatchCol::from_string("count_um".to_string()), UBatchCol::CountUm);
        assert_eq!(UBatchCol::from_string("density".to_string()), UBatchCol::Density);
    }

    #[test]
    #[should_panic(expected = "Unknown BSX data type invalid")]
    fn test_from_string_invalid() {
        UBatchCol::from_string("invalid".to_string());
    }
}


#[derive(Debug, Clone)]
pub struct UBatch {
    data: DataFrame,
}

impl UBatch {
    pub fn new(mut data: DataFrame) -> Self {
        let mut data_schema = data.schema();
        assert!(
            Self::colnames().into_iter().all(|name| data_schema.get(&*name).is_some()),
            "DataFrame columns does not match schema ({:?})", data_schema
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
    
    fn colnames() -> Vec<String> {
        UBatchCol::all_cols()
    }
    
    pub fn row_count(&self) -> usize {
        self.data.height()
    }
    
    /// Creates a [RegionCoordinates] objects from [UBatch] positions.
    /// 
    /// # Warning
    /// This method does not perform any checks
    pub fn get_region_coordinates(&self) -> RegionCoordinates {
        RegionCoordinates::new(self.get_chr(), self.first_position(), self.last_position())
    }
    
    pub fn partition_by(mut self, column: UBatchCol, include_key: bool) -> Vec<UBatch> {
        let new_data = self.data.partition_by([column.as_str()], include_key).unwrap();
        new_data.into_iter().map(|df| UBatch::from(df)).collect()
    }
    
    pub fn split(&self, offset: i64) -> (UBatch, UBatch) {
        let (first, second) = self.data.split_at(offset);
        (first.into(), second.into())
    }
    
    /// This method checks if the internal data is sorted by position and performs sorting
    /// if data was not sorted
    /// 
    /// # Raises
    /// If data has different chromosomes
    pub fn check_pos_sorted(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        if !self.unique_chr() { return Err(Box::from("Chromosomes differ!")) }
        
        let was_sorted = self.data.column("position")?.u32()?.iter().is_sorted();
        if !was_sorted {
            self.data.sort(["position"], SortMultipleOptions::default())?;
        }
        Ok(())
    }
    
    fn _strand_expr() -> Expr {
        when(col("strand").eq(lit("+"))).then(true).when(col("strand").eq(lit("-"))).then(false).otherwise(lit(NULL)).cast(DataType::Boolean)
    }
    fn _context_expr() -> Expr {
        when(col("context").eq(lit("CG")))
            .then(lit(true))
            .when(col("context").eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
    }
    fn _nuc_expr() -> Expr {
        when(col("nuc").eq(lit("C"))).then(true).when(col("strand").eq(lit("G"))).then(false).otherwise(lit(NULL)).cast(DataType::Boolean)
    }

    pub fn from_report_type(data: DataFrame, report_type: &ReportType) -> Result<Self, Box<dyn std::error::Error>> {
        let lazy_frame = data.lazy();
        
        let res = match report_type {
            ReportType::BISMARK => lazy_frame
                .with_column((col("count_m") + col("count_um")).alias("count_total"))
                .with_columns([
                    (col("count_m") / col("count_total").cast(DataType::Float64)).cast(DataType::Float64).alias("density"),
                    Self::_strand_expr().alias("strand"),
                    Self::_context_expr().alias("context"),
                ]),
            ReportType::CGMAP => lazy_frame.with_columns([
                Self::_nuc_expr().alias("strand"),
                Self::_context_expr().alias("context"),
            ]),
            ReportType::BEDGRAPH => lazy_frame
                .rename(["start"], ["position"], true)
                .drop(["end"])
                .with_columns([
                    lit(NULL).alias("strand"),
                    lit(NULL).alias("context"),
                    lit(NULL).alias("count_m"),
                    lit(NULL).alias("count_total"),
                    col("density").div(lit(100)).alias("density"),
                ]),
            ReportType::COVERAGE => lazy_frame
                .rename(["start"], ["position"], true)
                .drop(["end"])
                .with_column((col("count_m") + col("count_um")).alias("count_total"))
                .with_columns([
                    lit(NULL).alias("strand"),
                    lit(NULL).alias("context"),
                    col("density").div(lit(100)).alias("density"),
                ]),
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
    pub fn partition_chr(self) -> Vec<UBatch> {
        self.data
            .partition_by(["chr"], true).unwrap()
            .into_iter()
            .map(|d| UBatch::new(d))
            .collect()
    }
    pub fn partition_by_chunks(self, chunk_size: usize) -> Vec<UBatch> {
        self.data
            .lazy()
            .with_row_index("index", None)
            .with_column(col("index").floor_div(lit(chunk_size as u32)))
            .collect().unwrap()
            .partition_by(["index"], false)
            .unwrap()
            .into_iter()
            .map(|d| UBatch::new(d))
            .collect()
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
    
    pub fn vstack(mut self, batch: &UBatch) -> UBatch {
        self.data = self.data.vstack(&batch.data).unwrap();
        self
    }
    
    pub fn extend(&mut self, batch: &UBatch, to_end: bool) {
        if to_end {
            self.data.extend(&batch.data).unwrap()
        } else {
            self.data = batch.get_data().clone().vstack(&self.data).unwrap()
        }
    }
    
    /// Create new [UBatch] by [RegionCoordinates] slice.
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

impl Display for UBatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.data)
    }
}

impl From<DataFrame> for UBatch {
    fn from(data_frame: DataFrame) -> Self {
        UBatch::new(data_frame)
    }
}

impl From<UBatch> for DataFrame {
    fn from(data_frame: UBatch) -> Self {
        data_frame.data
    }
}

