use crate::io::report::types::ReportType;
use crate::region::RegionCoordinates;
use crate::utils::types::BSXResult;
use arrow::datatypes::ArrowNativeType;
use log::debug;
use polars::prelude::*;
use std::error::Error;
use std::fmt::Display;
use std::ops::Div;

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum BSXCol {
    Chr,
    Strand,
    Position,
    Context,
    CountM,
    CountTotal,
    Density,
}

impl BSXCol {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Chr => "chr",
            Self::Strand => "strand",
            Self::Position => "position",
            Self::Context => "context",
            Self::CountM => "count_m",
            Self::CountTotal => "count_total",
            Self::Density => "density",
        }
    }

    pub fn data_type(&self) -> DataType {
        match self {
            Self::Chr => DataType::Categorical(None, CategoricalOrdering::Physical),
            Self::Strand => DataType::Boolean,
            Self::Position => DataType::UInt32,
            Self::Context => DataType::Boolean,
            Self::CountM => DataType::UInt16,
            Self::CountTotal => DataType::UInt16,
            Self::Density => DataType::Float64,
        }
    }

    pub fn all_cols() -> Vec<String> {
        vec![
            Self::Chr.into(),
            Self::Strand.into(),
            Self::Position.into(),
            Self::Context.into(),
            Self::CountM.into(),
            Self::CountTotal.into(),
            Self::Density.into(),
        ]
    }

    fn all_variants() -> Vec<BSXCol> {
        vec![
            BSXCol::Chr,
            BSXCol::Strand,
            BSXCol::Position,
            BSXCol::Context,
            BSXCol::CountM,
            BSXCol::CountTotal,
            BSXCol::Density,
        ]
    }

    fn schema() -> Schema {
        Schema::from_iter([
            (Self::Chr.as_str().into(), Self::Chr.data_type()),
            (Self::Strand.as_str().into(), Self::Strand.data_type()),
            (Self::Position.as_str().into(), Self::Position.data_type()),
            (Self::Context.as_str().into(), Self::Context.data_type()),
            (Self::CountM.as_str().into(), Self::CountM.data_type()),
            (
                Self::CountTotal.as_str().into(),
                Self::CountTotal.data_type(),
            ),
            (Self::Density.as_str().into(), Self::Density.data_type()),
        ])
    }

    fn can_cast(&self, data_type: DataType) -> bool {
        match self {
            Self::Chr => data_type.is_enum(),
            Self::Strand => data_type.is_bool(),
            Self::Position => data_type.is_unsigned_integer(),
            Self::Context => data_type.is_bool(),
            Self::CountM => data_type.is_unsigned_integer(),
            Self::CountTotal => data_type.is_unsigned_integer(),
            Self::Density => data_type.is_float(),
        }
    }

    fn is_non_null(&self) -> bool {
        match self {
            Self::Chr => true,
            Self::Strand => false,
            Self::Position => true,
            Self::Context => false,
            Self::CountM => false,
            Self::CountTotal => false,
            Self::Density => false,
        }
    }

    fn from_string(string: String) -> Self {
        match string.as_str() {
            "chr" => Self::Chr,
            "strand" => Self::Strand,
            "position" => Self::Position,
            "context" => Self::Context,
            "count_m" => Self::CountM,
            "count_total" => Self::CountTotal,
            "density" => Self::Density,
            _ => panic!("Unknown BSX data type {}", string),
        }
    }
}

impl From<BSXCol> for DataType {
    fn from(col: BSXCol) -> Self {
        col.data_type()
    }
}

impl From<BSXCol> for Expr {
    fn from(column: BSXCol) -> Self {
        col(PlSmallStr::from(column))
    }
}

impl From<BSXCol> for PlSmallStr {
    fn from(column: BSXCol) -> Self {
        column.as_str().into()
    }
}

impl From<BSXCol> for String {
    fn from(column: BSXCol) -> Self {
        column.as_str().to_string()
    }
}

/// Trait for all internal batch formats of BSXplorer.
/// All batches will have columns from [BSXCol]
pub trait BSXBatch {
    /* --- Constructors --- */
    /// Unsafely construct [BSXBatch].
    fn from_df(data_frame: DataFrame) -> Self;
    /// Get reference to inner [DataFrame].

    /* --- Getters --- */
    fn get_data(&self) -> &DataFrame;
    /// Get mutable reference to inner [DataFrame].
    fn get_data_mut(&mut self) -> &mut DataFrame;
    /// Get first position in batch
    fn get_first_pos(&self) -> u32 {
        self.get_data()
            .column(BSXCol::Position.as_str())
            .unwrap()
            .u32()
            .unwrap()
            .first()
            .unwrap()
    }
    /// Get last position in batch
    fn get_last_pos(&self) -> u32 {
        self.get_data()
            .column(BSXCol::Position.as_str())
            .unwrap()
            .u32()
            .unwrap()
            .last()
            .unwrap()
    }
    /// Get mapping of chromosomes and their categorical indices
    fn get_chr_revmap(&self) -> &Arc<RevMapping> {
        self.get_data()
            .column(BSXCol::Chr.as_str())
            .unwrap()
            .categorical()
            .unwrap()
            .get_rev_map()
    }
    /// Get index of chromosome from first row of the batch
    fn get_chr_idx(&self) -> usize {
        self.get_data()
            .column(BSXCol::Chr.as_str())
            .unwrap()
            .categorical()
            .unwrap()
            .physical()
            .first()
            .unwrap()
            .as_usize()
    }
    /// Get name of the first chromosome in the batch
    fn get_chr(&self) -> Option<String> {
        if self.check_chr_unique().is_ok() {
            let rev_map = self.get_chr_revmap();
            let first_idx = self.get_chr_idx();
            Some(rev_map.get(first_idx as u32).to_string())
        } else {
            None
        }
    }

    fn get_region_coordinates(&self) -> Option<RegionCoordinates> {
        if let Some(chr) = self.get_chr() {
            RegionCoordinates::new(chr, self.get_first_pos(), self.get_last_pos()).into()
        } else {
            None
        }
    }

    /// Get row count
    fn length(&self) -> usize {
        self.get_data().height()
    }

    fn _set_data(&mut self, data_frame: DataFrame);

    /* --- Checks --- */
    fn check_types(&self) -> Result<(), Box<dyn Error>> {
        let data = self.get_data();
        let data_schema = data.schema();
        let ucol_all = BSXCol::all_variants();
        // Check all columns are coercible
        for c in ucol_all.iter() {
            let dtype = data_schema.get(c.as_str()).unwrap().clone();
            if !c.can_cast(dtype.clone()) {
                return Err(Box::from(format!("Can not cast {c} from {dtype} to {c}")));
            }
        }
        // Check non-null
        for c in ucol_all.iter().filter(|c| c.is_non_null()) {
            let null_count = data.column(c.as_str()).unwrap().null_count();
            if null_count != 0 {
                return Err(Box::from(format!("DataFrame column '{c}' cannot be null")));
            }
        }
        Ok(())
    }
    fn check_chr_unique(&self) -> BSXResult<()> {
        if self
            .get_data()
            .column(BSXCol::Chr.as_str())
            .unwrap()
            .categorical()
            .unwrap()
            .physical()
            .unique()
            .unwrap()
            .len()
            == 1
        {
            Ok(())
        } else {
            Err(Box::from("Chromosome not unique"))
        }
    }
    fn check_cols(&self) -> Result<(), Box<dyn Error>> {
        let data = self.get_data();
        let data_schema = data.schema();
        let ucol_all = BSXCol::all_variants();

        // Check all cols present
        for c in ucol_all.iter() {
            if !data_schema.contains(c.as_str()) {
                return Err(Box::from(format!("Field {c} is missing")));
            }
        }
        Ok(())
    }
    fn check_position_sorted(&self) -> BSXResult<()> {
        if self
            .get_data()
            .column(BSXCol::Position.as_str())
            .unwrap()
            .u32()
            .unwrap()
            .iter()
            .is_sorted()
        {
            Ok(())
        } else {
            Err("Data is not sorted".into())
        }
    }

    fn from_cast(data: LazyFrame, chr_names: &Vec<String>) -> Result<Self, Box<dyn Error>>
    where
        Self: Sized,
    {
        let mut bsx_schema = BSXCol::schema();

        let categories = DataType::Enum(
            Some({
                let mut cat_builder =
                    CategoricalChunkedBuilder::new(Default::default(), 0, Default::default());
                for chr_name in chr_names {
                    let _ = &cat_builder.append(Some(chr_name.as_str()));
                }
                cat_builder.finish().get_rev_map().clone()
            }),
            CategoricalOrdering::Physical,
        );
        assert!(categories.contains_categoricals());

        bsx_schema.with_column(BSXCol::Chr.into(), categories);
        debug!("New schema {bsx_schema:?}");

        // TODO dodelat'

        let new = data
            .clone()
            .lazy()
            .cast(
                PlHashMap::from_iter(
                    bsx_schema
                        .iter_names_and_dtypes()
                        .map(|(n, dtype)| (n.as_str(), dtype.clone())),
                ),
                true,
            )
            .collect();
        match new {
            Ok(df) => Ok(Self::from_df(df)),
            Err(e) => Err(e.into()),
        }
    }

    fn sort_positions(&mut self) {
        let sorted = self
            .get_data()
            .sort(["position"], SortMultipleOptions::default())
            .unwrap();
        self._set_data(sorted);
    }

    /* --- Mutations --- */
    /// Modify [BSXBatch] inplace by inserting rows
    fn extend(&mut self, batch: &UBatch) {
        self.get_data_mut().extend(&batch.get_data()).unwrap();
    }
    /// Create new [BSXBatch] with concatenated data
    fn vstack(&self, batch: &Self) -> Self
    where
        Self: Sized,
    {
        Self::from_df(self.get_data().vstack(&batch.get_data()).unwrap())
    }
    fn filter(&self, expr: Expr) -> PolarsResult<Self>
    where
        Self: Sized,
    {
        match self.get_data().clone().lazy().filter(expr).collect() {
            Ok(df) => Ok(Self::from_df(df)),
            Err(e) => Err(e),
        }
    }

    /* --- Partitions --- */
    /// Partition batch by [BSXCol]
    fn partition_by_bsx(&self, column: BSXCol) -> Vec<UBatch> {
        let new_data = self
            .get_data()
            .partition_by([column.as_str()], true)
            .unwrap();
        new_data.into_iter().map(|df| UBatch::from(df)).collect()
    }
    /// Partition by custom column name. Name is not guarantied
    /// to be present in internal [DataFrame].
    fn partition_by(&mut self, column: String) -> PolarsResult<Vec<UBatch>> {
        match self.get_data().partition_by([column.as_str()], true) {
            Ok(cols) => Ok(cols.into_iter().map(|df| UBatch::from(df)).collect()),
            Err(e) => Err(e),
        }
    }
    fn split(&self, offset: i64) -> (Self, Self)
    where
        Self: Sized,
    {
        let (first, second) = self.get_data().split_at(offset);
        (Self::from_df(first), Self::from_df(second))
    }

    /* --- Conversion --- */
    fn _strand_expr() -> Expr {
        when(col("strand").eq(lit("+")))
            .then(true)
            .when(col("strand").eq(lit("-")))
            .then(false)
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
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
        when(col("nuc").eq(lit("C")))
            .then(true)
            .when(col("strand").eq(lit("G")))
            .then(false)
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
    }

    fn from_report_type(data: DataFrame, report_type: &ReportType) -> Result<Self, Box<dyn Error>>
    where
        Self: Sized,
    {
        let df = report_type_to_bsx(report_type, data)?;
        Ok(Self::from_df(df))
    }

    fn get_report_type(&self, report_type: ReportType) -> DataFrame {
        let target_schema = report_type.get_schema();
        let mut data_lazy = self.get_data().clone().lazy();
        if target_schema.get("strand").is_some() {
            data_lazy = data_lazy.with_column(
                when(Expr::from(BSXCol::Strand).eq(lit(false)))
                    .then(lit("-"))
                    .when(col("strand").eq(lit(true)))
                    .then(lit("+"))
                    .otherwise(lit(".")),
            )
        }
        if target_schema.get("context").is_some() {
            data_lazy = data_lazy.with_column(
                when(Expr::from(BSXCol::Context).is_null())
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
}

impl Display for BSXCol {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

#[derive(Debug, Clone)]
pub struct UBatch {
    data: DataFrame,
}

impl BSXBatch for UBatch {
    fn from_df(batch: DataFrame) -> Self {
        let obj = UBatch { data: batch };
        obj
    }
    fn get_data(&self) -> &DataFrame {
        &self.data
    }
    fn get_data_mut(&mut self) -> &mut DataFrame {
        &mut self.data
    }

    fn _set_data(&mut self, data_frame: DataFrame) {
        self.data = data_frame;
    }
}

impl UBatch {
    fn new(data: DataFrame) -> Self {
        let obj = UBatch { data };
        match obj.check_types() {
            Ok(_) => obj,
            Err(e) => {
                panic!("{}", e)
            }
        }
    }
    
    pub fn into_data(self) -> DataFrame {
        self.data
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

fn _strand_expr() -> Expr {
    when(col("strand").eq(lit("+")))
        .then(true)
        .when(col("strand").eq(lit("-")))
        .then(false)
        .otherwise(lit(NULL))
        .cast(BSXCol::Strand.data_type())
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
    when(col("nuc").eq(lit("C")))
        .then(true)
        .when(col("strand").eq(lit("G")))
        .then(false)
        .otherwise(lit(NULL))
        .cast(DataType::Boolean)
}

/// Calculates obligatory columns from [BSXBatch], performs cast for
/// all columns except [BSXCol::Chr]
pub fn report_type_to_bsx(report_type: &ReportType, df: DataFrame) -> BSXResult<DataFrame> {
    let lazy_frame = df.lazy();

    let res = match report_type {
        ReportType::BISMARK => lazy_frame
            // count_total
            .with_column(
                (col("count_m") + col("count_um"))
                    .cast(BSXCol::CountTotal.data_type())
                    .alias(BSXCol::CountTotal.as_str()),
            )
            // density
            .with_column(
                (col("count_m").div(BSXCol::CountTotal.into()))
                    .cast(BSXCol::Density.data_type())
                    .alias(BSXCol::Density.as_str()),
            )
            .with_columns([
                _strand_expr().alias(BSXCol::Strand.as_str()),
                _context_expr().alias(BSXCol::Context.as_str()),
            ]),
        ReportType::CGMAP => lazy_frame
            .with_columns([
                _nuc_expr().alias("strand"),
                _context_expr().alias("context"),
            ]),
        ReportType::BEDGRAPH => lazy_frame
            .rename(["start"], ["position"], true)
            .drop(["end"])
            .with_columns([
                lit(NULL).alias("strand").cast(BSXCol::Strand.data_type()),
                lit(NULL).alias("context").cast(BSXCol::Context.data_type()),
                lit(NULL).alias("count_m").cast(BSXCol::CountM.data_type()),
                lit(NULL).alias("count_total").cast(BSXCol::CountTotal.data_type()),
                col("density").div(lit(100)).alias("density").cast(BSXCol::Density.data_type()),
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
    }.cast(
        PlHashMap::from_iter(vec![
            (BSXCol::Position.as_str(), BSXCol::Position.data_type()),
            (BSXCol::Strand.as_str(), BSXCol::Strand.data_type()),
            (BSXCol::Context.as_str(), BSXCol::Context.data_type()),
            (BSXCol::Density.as_str(), BSXCol::Density.data_type()),
            (BSXCol::CountTotal.as_str(), BSXCol::CountTotal.data_type()),
            (BSXCol::CountM.as_str(), BSXCol::CountM.data_type()),
        ]), false
    )
    .collect()?;

    Ok(res)
}

#[cfg(test)]
mod tests {
    use std::assert_matches::assert_matches;
    use crate::io::report::readv3::{cast_chr_to_categorical, get_enum_chr};
    use super::*;

    struct TestBatch {
        data: DataFrame,
    }

    impl Default for TestBatch {
        fn default() -> Self {
            let data = df!(
                "chr" => &["chr1", "chr1", "chr2", "chr2", "chr2"],
                "strand" => &[Some(true), Some(false), Some(false), Some(true), Some(true)],
                "position" => &[1, 2, 6, 7, 8],
                "context" => &[Some(true), Some(false), None, Some(true), None],
                "count_m" => &[10, 20, 30, 0, 0],
                "count_total" => &[10, 20, 30, 40, 50],
                "density" => &[1.0, 1.0, 1.0, 0.0, 0.0],
            ).unwrap().lazy().cast(
                PlHashMap::from_iter([
                    (BSXCol::Chr.as_str(), DataType::String),
                    (BSXCol::Strand.as_str(), BSXCol::Strand.data_type()),
                    (BSXCol::Position.as_str(), BSXCol::Position.data_type()),
                    (BSXCol::Context.as_str(), BSXCol::Context.data_type()),
                    (BSXCol::CountM.as_str(), BSXCol::CountM.data_type()),
                    (BSXCol::CountTotal.as_str(), BSXCol::CountTotal.data_type()),
                    (BSXCol::Density.as_str(), BSXCol::Density.data_type()),
                ]), true
            ).collect().unwrap();
            let chr_type = get_enum_chr(&vec![String::from("chr1"), String::from("chr2")]);
            let data = cast_chr_to_categorical(data, &chr_type).unwrap();
            Self {data}
        }
    }

    impl BSXBatch for TestBatch {
        fn from_df(data_frame: DataFrame) -> Self {
            Self {data: data_frame}
        }

        fn get_data(&self) -> &DataFrame {
            &self.data
        }

        fn get_data_mut(&mut self) -> &mut DataFrame {
            &mut self.data
        }

        fn _set_data(&mut self, data_frame: DataFrame) {
            self.data = data_frame;
        }
    }

    #[test]
    fn test_pos_retr() {
        let batch = TestBatch::default();
        assert_eq!(batch.get_first_pos(), 1);
        assert_eq!(batch.get_last_pos(), 8);
    }

    #[test]
    fn test_chr_retr() {
        let batch = TestBatch::default();
        // Splitting to unique chr
        let mut partitioned = batch.partition_by_bsx(BSXCol::Chr);
        partitioned.sort_by_cached_key(|b| b.get_chr().unwrap());
        let batch = partitioned.first().unwrap();
        assert_eq!(batch.get_chr().unwrap().as_str(), "chr1");
        assert_eq!(batch.get_chr_idx(), 0);

    }

    #[test]
    fn test_diff_chr() {
        let batch = TestBatch::default();
        assert_eq!(batch.get_chr(), None);
    }
    
    #[test]
    fn region_coordinates() {
        let batch = TestBatch::default();
        // Splitting to unique chr
        let mut partitioned = batch.partition_by_bsx(BSXCol::Chr);
        partitioned.sort_by_cached_key(|b| b.get_chr().unwrap());
        let batch = partitioned.first().unwrap();
        let correct = RegionCoordinates::new("chr1".into(), 1, 2);
        assert_eq!(batch.get_region_coordinates().unwrap(), correct)
    }
    
    #[test]
    fn length() {
        let batch = TestBatch::default();
        assert_eq!(batch.length(), 5)
    }
    
    #[test]
    fn test_chr_unique() {
        let batch = TestBatch::default();
        assert_matches!(batch.check_chr_unique(), Err(_));
    }
}