use std::any::Any;
use std::ops::Deref;

use hashbrown::HashMap;
use itertools::{izip, Itertools};
use num::Zero;
use polars::frame::column::ScalarColumn;
use polars::prelude::*;
use statrs::distribution::DiscreteCDF;
use statrs::statistics::Statistics;

use crate::data_structs::context_data::ContextData;
use crate::data_structs::coords::{Contig, GenomicPosition};
use crate::data_structs::enums::{Context, IPCEncodedEnum, Strand};
use crate::data_structs::methstats::MethylationStats;
use crate::data_structs::typedef::BsxSmallStr;
use crate::io::report::ReportType;
use crate::{plsmallstr, with_field_fn};

macro_rules! create_empty_series {
    ($col: ident) => {
        Series::new_empty(plsmallstr!(BsxCol::$col.as_str()), &BsxCol::$col.dtype())
    };
}

enum BsxColumns {
    Chr,
    Position,
    Strand,
    Context,
    CountM,
    CountTotal,
    Density,
}

impl BsxColumns {
    pub const fn as_str(&self) -> &'static str {
        match self {
            BsxCol::Chr => "chr",
            BsxCol::Position => "position",
            BsxCol::Strand => "strand",
            BsxCol::Context => "context",
            BsxCol::CountM => "count_m",
            BsxCol::CountTotal => "count_total",
            BsxCol::Density => "density",
        }
    }

    pub const fn dtype(&self) -> DataType {
        match self {
            BsxCol::Chr => create_empty_categorical_dtype(),
            BsxCol::Position => DataType::UInt32,
            BsxCol::Strand => DataType::Boolean,
            BsxCol::Context => DataType::Boolean,
            BsxCol::CountM => DataType::UInt16,
            BsxCol::CountTotal => DataType::UInt16,
            BsxCol::Density => DataType::Float32,
        }
    }

    pub const fn colnames() -> [&'static str; 7] {
        [
            BsxCol::Chr.as_str(),
            BsxCol::Position.as_str(),
            BsxCol::Strand.as_str(),
            BsxCol::Context.as_str(),
            BsxCol::CountM.as_str(),
            BsxCol::CountTotal.as_str(),
            BsxCol::Density.as_str(),
        ]
    }

    pub fn has_name(name: &str) -> bool {
        Self::colnames().contains(&name)
    }

    #[inline(always)]
    pub fn col(&self) -> Expr {
        col(self.as_str())
    }

    pub fn create_anyvalue(
        &self,
        value: Box<dyn Any>,
    ) -> Option<AnyValue> {
        match self {
            BsxCol::Chr => {
                value
                    .downcast_ref::<String>()
                    .map(|v| AnyValue::StringOwned(v.into()))
            },
            BsxCol::Position => {
                value.downcast_ref::<u32>().map(|v| AnyValue::UInt32(*v))
            },
            BsxCol::Strand => {
                value.downcast_ref::<bool>().map(|v| AnyValue::Boolean(*v))
            },
            BsxCol::Context => {
                value.downcast_ref::<Option<bool>>().map(|v| {
                    v.map(|ctx| AnyValue::Boolean(ctx))
                        .unwrap_or(AnyValue::Null)
                })
            },
            BsxCol::CountM => value.downcast_ref::<u16>().map(|v| AnyValue::UInt16(*v)),
            BsxCol::CountTotal => {
                value.downcast_ref::<u16>().map(|v| AnyValue::UInt16(*v))
            },
            BsxCol::Density => {
                value.downcast_ref::<f32>().map(|v| AnyValue::Float32(*v))
            },
        }
    }

    pub fn create_series<T>(
        &self,
        data: Vec<T>,
    ) -> PolarsResult<Series>
    where
        T: Sized + 'static, {
        let any_vec = data
            .into_iter()
            .map(|value| self.create_anyvalue(Box::new(value)))
            .collect::<Option<Vec<_>>>()
            .ok_or(PolarsError::SchemaMismatch(
                format!("Could not downcast type {}", stringify!(T)).into(),
            ))?;
        Series::from_any_values(plsmallstr!(self.as_str()), &any_vec, true)
    }
}

use BsxColumns as BsxCol;

const fn create_empty_categorical_dtype() -> DataType {
    DataType::Categorical(None, CategoricalOrdering::Physical)
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BsxBatchBuilder {
    pub chr_dtype:        Option<DataType>,
    pub check_nulls:      bool,
    pub check_sorted:     bool,
    pub check_duplicates: bool,
    pub rechunk:          bool,
    pub check_single_chr: bool,
}

impl Default for BsxBatchBuilder {
    fn default() -> Self {
        Self::all_checks()
    }
}

macro_rules! name_dtype_tuple {
    ($enum_var: expr) => {
        ($enum_var.as_str().into(), $enum_var.dtype())
    };
}

// Public methods
impl BsxBatchBuilder {
    with_field_fn!(chr_dtype, Option<DataType>);

    with_field_fn!(check_duplicates, bool);

    with_field_fn!(check_sorted, bool);

    with_field_fn!(rechunk, bool);

    with_field_fn!(check_nulls, bool);

    with_field_fn!(check_single_chr, bool);

    /// Creates a builder with all data validation checks enabled.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn all_checks() -> Self {
        Self {
            chr_dtype:        None,
            check_duplicates: true,
            check_sorted:     true,
            rechunk:          true,
            check_nulls:      true,
            check_single_chr: true,
        }
    }

    /// Creates a builder with all validation checks disabled.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn no_checks() -> Self {
        Self {
            chr_dtype:        None,
            check_duplicates: false,
            check_sorted:     false,
            rechunk:          false,
            check_nulls:      false,
            check_single_chr: false,
        }
    }

    pub fn with_chr_values<S, P>(
        mut self,
        chr_values: Vec<String>,
    ) -> Self
    where
        S: AsRef<str>,
        P: AsRef<[S]>, {
        let dtype = Some(create_caregorical_dtype(
            chr_values.into_iter().map(|v| Some(v)).collect_vec(),
        ));
        self.chr_dtype = dtype;
        self
    }

    pub fn checks_only(
        &self,
        df: &DataFrame,
    ) -> PolarsResult<()> {
        if self.check_nulls && !build::check_no_nulls(df)? {
            return Err(PolarsError::NoData(
                "Null values found in 'position' column".into(),
            ));
        }
        if self.check_single_chr && !build::check_single_chr(df)? {
            return Err(PolarsError::InvalidOperation(
                "There should be exactly one chromosome value per batch".into(),
            ));
        }
        if self.check_sorted && !build::check_sorted(df)? {
            return Err(PolarsError::InvalidOperation(
                "The 'position' column is not sorted".into(),
            ));
        }
        if self.check_duplicates && !build::check_no_duplicates(df)? {
            return Err(PolarsError::Duplicate(
                "'position' column contains duplicates".into(),
            ));
        }
        Ok(())
    }

    pub fn cast_only(
        &self,
        df: DataFrame,
    ) -> PolarsResult<DataFrame> {
        let mut df = df;
        build::set_flags(&mut df)?;
        df.lazy()
            .cast(
                HashMap::from_iter([
                    (
                        BsxCol::Chr.as_str().into(),
                        self.chr_dtype
                            .clone()
                            .unwrap_or(create_empty_categorical_dtype()),
                    ),
                    name_dtype_tuple!(BsxCol::Position),
                    name_dtype_tuple!(BsxCol::Strand),
                    name_dtype_tuple!(BsxCol::Context),
                    name_dtype_tuple!(BsxCol::CountM),
                    name_dtype_tuple!(BsxCol::CountTotal),
                    name_dtype_tuple!(BsxCol::Density),
                ]),
                true,
            )
            .collect()
    }

    pub fn build_from_report(
        &self,
        df: DataFrame,
        report_type: ReportType,
    ) -> PolarsResult<BsxBatch> {
        let lf = df.lazy();
        let converted = match report_type {
            ReportType::Bismark => build::bismark_expr(lf),
            ReportType::CgMap => build::cgmap_expr(lf),
            ReportType::BedGraph => build::bedgraph_expr(lf),
            ReportType::Coverage => build::coverage_expr(lf),
        }
        .collect()?;

        self.checks_only(&converted)?;
        let casted = self.cast_only(converted)?;

        Ok(unsafe { BsxBatch::new_unchecked(casted) })
    }
}

pub fn create_caregorical_dtype<S, P>(chr_values: P) -> DataType
where
    S: AsRef<str>,
    P: AsRef<[Option<S>]>, {
    use polars::export::arrow::array::Utf8ViewArray;
    let categories = Utf8ViewArray::from_slice(chr_values);
    let rev_mapping = Arc::new(RevMapping::build_local(categories));
    DataType::Categorical(Some(rev_mapping), CategoricalOrdering::Physical)
}

mod build {
    use polars::series::IsSorted;

    use super::*;

    /// Check duplicates in the position column
    pub fn check_no_duplicates(df: &DataFrame) -> PolarsResult<bool> {
        Ok(df.column(BsxCol::Position.as_str())?.n_unique()? == df.height())
    }

    /// Check sortedness of the position column
    pub fn check_sorted(df: &DataFrame) -> PolarsResult<bool> {
        Ok(df
            .column(BsxCol::Position.as_str())?
            .u32()?
            .iter()
            .is_sorted())
    }

    /// Check nulls in the chromosome and position columns
    pub fn check_no_nulls(df: &DataFrame) -> PolarsResult<bool> {
        let chr_nulls = df.column(BsxCol::Chr.as_str())?.null_count();
        let pos_nulls = df.column(BsxCol::Position.as_str())?.null_count();
        Ok(chr_nulls == 0 && pos_nulls == 0)
    }

    /// Check that the chromosome column contains only a single value
    pub fn check_single_chr(df: &DataFrame) -> PolarsResult<bool> {
        let chr = df.column(BsxCol::Chr.as_str())?;
        Ok(chr.n_unique()? == 1)
    }

    pub fn set_flags(df: &mut DataFrame) -> PolarsResult<()> {
        let mut pos_col = df.column(BsxCol::Position.as_str())?.to_owned();
        pos_col.set_sorted_flag(IsSorted::Ascending);
        df.with_column(pos_col)?;
        Ok(())
    }

    /// Creates an expression to convert nucleotide values to strand symbols.
    pub fn nuc_to_bool_expr() -> Expr {
        when(col("nuc").eq(lit("C")))
            .then(lit(true))
            .when(col("nuc").eq(lit("G")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
    }

    pub fn strand_to_bool_expr() -> Expr {
        when(col("strand").eq(lit("+")))
            .then(lit(true))
            .when(col("strand").eq(lit("-")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
    }

    pub fn context_to_bool_expr() -> Expr {
        when(col("context").eq(lit("CG")))
            .then(lit(true))
            .when(col("context").eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
    }

    pub fn count_total_col_expr() -> Expr {
        col("count_m") + col("count_um")
    }

    pub fn density_col_expr() -> Expr {
        col("count_m").cast(DataType::Float64)
            / col("count_total").cast(DataType::Float64)
    }

    /// Processes Bismark format data into standardized BSX format.
    pub fn bismark_expr(lf: LazyFrame) -> LazyFrame {
        lf.with_column(count_total_col_expr().alias(BsxCol::CountTotal.as_str()))
            .with_columns([
                context_to_bool_expr().alias(BsxCol::Context.as_str()),
                strand_to_bool_expr().alias(BsxCol::Strand.as_str()),
                density_col_expr().alias(BsxCol::Density.as_str()),
            ])
            .select(BsxCol::colnames().into_iter().map(col).collect_vec())
    }

    /// Converts CG map format data into standardized BSX format.
    pub fn cgmap_expr(lf: LazyFrame) -> LazyFrame {
        lf.with_columns([
            nuc_to_bool_expr().alias(BsxCol::Strand.as_str()),
            nuc_to_bool_expr().alias(BsxCol::Context.as_str()),
        ])
        .select(BsxCol::colnames().into_iter().map(col).collect_vec())
    }

    /// Processes coverage format data into standardized BSX format.
    pub fn coverage_expr(lf: LazyFrame) -> LazyFrame {
        lf
            .with_columns([
                count_total_col_expr().alias(BsxCol::CountTotal.as_str()),
                col(ReportType::Coverage.position_col()).alias(BsxCol::Position.as_str()),
                lit(NULL).alias("strand"),
                lit(NULL).alias("context"),
            ])
            .select(BsxCol::colnames().into_iter().map(col).collect_vec())
    }

    /// Processes BedGraph format data into standardized BSX format.
    pub fn bedgraph_expr(lf: LazyFrame) -> LazyFrame {
        lf.with_columns([
            col(ReportType::BedGraph.position_col()).alias(BsxCol::Position.as_str()),
            lit(NULL).alias(BsxCol::Strand.as_str()),
            lit(NULL).alias(BsxCol::Context.as_str()),
            lit(NULL).alias(BsxCol::CountM.as_str()),
            lit(NULL).alias(BsxCol::CountTotal.as_str()),
        ])
        .select(BsxCol::colnames().into_iter().map(col).collect_vec())
    }

    #[cfg(test)]
    mod tests {
        use rstest::{fixture, rstest};
        use super::*;
        use crate::data_structs::enums::{Context, Strand};

        #[fixture]
        fn test_df() -> DataFrame {
            df!(
                BsxCol::Chr.as_str() => ["chr1", "chr1", "chr1"],
                BsxCol::Position.as_str() => [100u32, 150, 200],
                BsxCol::Strand.as_str() => [true, false, true],
                BsxCol::Context.as_str() => [true, false, true],
                BsxCol::CountM.as_str() => [5u16, 10, 15],
                BsxCol::CountTotal.as_str() => [10u16, 20, 30],
                BsxCol::Density.as_str() => [0.5f32, 0.5, 0.5]
            )
            .unwrap()
        }

        #[fixture]
        fn no_cols_df() -> DataFrame {
            DataFrame::empty()
        }

        #[rstest]
        #[case::pos_duplicates({
            df!(
                BsxCol::Chr.as_str() => ["chr1", "chr1", "chr1"],
                BsxCol::Position.as_str() => [100u32, 100, 200], // Duplicate at position 100
            ).unwrap()
        }, check_no_duplicates)]
        #[case::unsorted({
            df!(
                BsxCol::Position.as_str() => [200u32, 100, 150], // Unsorted positions
            ).unwrap()
        }, check_sorted)]
        #[case::nulls_chr({
            df!(
                BsxCol::Chr.as_str() => [Some("chr1"), Some("chr1"), None],
                BsxCol::Position.as_str() => [100u32, 150, 200],
            ).unwrap()
        }, check_no_nulls)]
        #[case::nulls_pos({
            df!(
                BsxCol::Chr.as_str() => ["chr1", "chr1", "chr1"],
                BsxCol::Position.as_str() => [Some(100), Some(150), None],
            ).unwrap()
        }, check_no_nulls)]
        #[case::multi_chr({
            df!(
                BsxCol::Chr.as_str() => ["chr1", "chr2", "chr3"],
                BsxCol::Position.as_str() => [100u32, 150, 200],
            ).unwrap()
        }, check_single_chr)]
        fn test_checks(
            test_df: DataFrame,
            no_cols_df: DataFrame,
            #[case] should_not_pass_df: DataFrame,
            #[case] check_fn: fn(&DataFrame) -> PolarsResult<bool>
        ) {
            assert_eq!(check_fn(&test_df).unwrap(), true);
            assert_eq!(check_fn(&should_not_pass_df).unwrap(), false);
            assert!(check_fn(&no_cols_df).is_err());
        }

        #[rstest]
        fn test_set_flags(mut test_df: DataFrame) {
            assert!(set_flags(&mut test_df).is_ok());

            // Verify the sorting flag was set
            let pos_col = test_df.column(BsxCol::Position.as_str()).unwrap();
            assert_eq!(pos_col.is_sorted_flag(), IsSorted::Ascending);
        }

        #[test]
        fn test_conversion_expressions() {
            // Test nuc_to_bool_expr
            let df = df!("nuc" => ["C", "G", "A"]).unwrap();
            let res = df.lazy()
                .with_column(nuc_to_bool_expr().alias("result"))
                .collect()
                .unwrap();
            let results = res.column("result").unwrap().bool().unwrap();
            assert_eq!(results.get(0), Some(true));
            assert_eq!(results.get(1), Some(false));
            assert_eq!(results.get(2), None);

            // Test strand_to_bool_expr
            let df = df!("strand" => ["+", "-", "."]).unwrap();
            let res = df.lazy()
                .with_column(strand_to_bool_expr().alias("result"))
                .collect()
                .unwrap();
            let results = res.column("result").unwrap().bool().unwrap();
            assert_eq!(results.get(0), Some(true));
            assert_eq!(results.get(1), Some(false));
            assert_eq!(results.get(2), None);

            // Test context_to_bool_expr
            let df = df!("context" => ["CG", "CHG", "CHH"]).unwrap();
            let res = df.lazy()
                .with_column(context_to_bool_expr().alias("result"))
                .collect()
                .unwrap();
            let results = res.column("result").unwrap().bool().unwrap();
            assert_eq!(results.get(0), Some(true));
            assert_eq!(results.get(1), Some(false));
            assert_eq!(results.get(2), None);
        }

        #[test]
        fn test_count_total_and_density() {
            let df = df!(
                "count_m" => [5i64, 10],
                "count_um" => [5i64, 10]
            ).unwrap();

            // Test count_total_col_expr
            let res = df.clone().lazy()
                .with_column(count_total_col_expr().alias("count_total"))
                .collect()
                .unwrap();
            let totals = res.column("count_total").unwrap().i64().unwrap();
            assert_eq!(totals.get(0), Some(10));
            assert_eq!(totals.get(1), Some(20));

            // Test density_col_expr
            let res = df.lazy()
                .with_column(count_total_col_expr().alias("count_total"))
                .with_column(density_col_expr().alias("density"))
                .collect()
                .unwrap();
            let densities = res.column("density").unwrap().f64().unwrap();
            assert_eq!(densities.get(0), Some(0.5));
            assert_eq!(densities.get(1), Some(0.5));
        }

        fn create_bismark_df() -> DataFrame {
            df!(
                "chr" => ["chr1", "chr1"],
                "position" => [100u32, 150],
                "strand" => ["+", "-"], // Decoded uses "+", "-"
                "count_m" => [5u32, 10],
                "count_um" => [5u32, 10],
                "context" => ["CG", "CHG"],
                "trinuc" => ["CGA", "CAG"] // Bismark requires trinuc column
            )
            .unwrap()
        }

        fn create_cgmap_df() -> DataFrame {
            df!(
                "chr" => ["chr2", "chr2"],
                "nuc" => ["C", "G"], // Determines strand ('C' -> '+', 'G' -> '-')
                "position" => [200u32, 250],
                "context" => ["CG", "CHH"],
                "dinuc" => ["CG", "CA"], // CgMap requires dinuc, not pattern
                "density" => [0.8f64, 0.4], // CgMap has density column
                "count_m" => [8u32, 2],
                "count_total" => [10u32, 5]
            )
            .unwrap()
        }

        fn create_coverage_df() -> DataFrame {
            df!(
                "chr" => ["chr3", "chr3"],
                "start" => [300u32, 350], // Becomes position
                "end" => [301u32, 351],
                "density" => [0.75f64, 0.4], // Coverage has density column
                "count_m" => [15u32, 20],
                "count_um" => [5u32, 30] // Used to calculate count_total
            )
            .unwrap()
        }

        fn create_bedgraph_df() -> DataFrame {
            df!(
                "chr" => ["chr4", "chr4"],
                "start" => [400u32, 450], // Becomes position
                "end" => [401u32, 451],
                "density" => [0.75f64, 0.8] // Only density provided
            )
            .unwrap()
        }

        #[rstest]
        #[case::bismark(ReportType::Bismark, create_bismark_df())]
        #[case::cgmap(ReportType::CgMap, create_cgmap_df())]
        #[case::coverage(ReportType::Coverage, create_coverage_df())]
        #[case::bedgraph(ReportType::BedGraph, create_bedgraph_df())]
        fn test_build_decoded_from_report_type(
            #[case] report_type: ReportType,
            #[case] input_df: DataFrame,
        ) -> anyhow::Result<()> {
            let batch = BsxBatchBuilder::all_checks().build_from_report(input_df, report_type)?;
            let _data = batch.data();
            Ok(())
        }
    }
}

pub struct BsxBatch {
    data: DataFrame,
}

macro_rules! get_col_fn {
    ($name: ident, $col: expr, $col_fn: ident, $rettype: ty) => {
        fn $name(&self) -> &$rettype {
            self.data().column($col).unwrap().$col_fn().unwrap()
        }
    };
}

impl BsxBatch {
    // COLUMN GETTERS
    get_col_fn!(chr, BsxCol::Chr.as_str(), categorical, CategoricalChunked);

    get_col_fn!(position, BsxCol::Position.as_str(), u32, UInt32Chunked);

    get_col_fn!(strand, BsxCol::Strand.as_str(), bool, BooleanChunked);

    get_col_fn!(context, BsxCol::Context.as_str(), bool, BooleanChunked);

    get_col_fn!(count_m, BsxCol::CountM.as_str(), u16, UInt16Chunked);

    get_col_fn!(count_total, BsxCol::CountTotal.as_str(), u16, UInt16Chunked);

    get_col_fn!(density, BsxCol::Density.as_str(), f32, Float32Chunked);

    // CONSTRUCTORS
    #[inline(always)]
    pub(crate) unsafe fn new_unchecked(df: DataFrame) -> Self {
        BsxBatch { data: df }
    }

    pub fn try_from_columns(
        chr: &str,
        chr_dtype: Option<DataType>,
        positions: Vec<u32>,
        strand: Vec<bool>,
        context: Vec<Option<bool>>,
        count_m: Vec<u16>,
        count_total: Vec<u16>,
        density: Vec<f32>,
    ) -> PolarsResult<Self> {
        assert!(
            [
                positions.len(),
                strand.len(),
                context.len(),
                count_m.len(),
                count_total.len(),
                density.len()
            ]
            .iter()
            .all_equal(),
            "All input vectors must have the same length"
        );

        let height = positions.len();
        let df = DataFrame::from_iter([
            {
                let name = plsmallstr!(BsxCol::Chr.as_str());
                let data = Scalar::new(
                    chr_dtype.unwrap_or(create_empty_categorical_dtype()),
                    AnyValue::StringOwned(plsmallstr!(chr)),
                );
                Column::Scalar(ScalarColumn::new(name, data, height))
            },
            BsxCol::Position.create_series(positions)?.into(),
            BsxCol::Strand.create_series(strand)?.into(),
            BsxCol::Context.create_series(context)?.into(),
            BsxCol::CountM.create_series(count_m)?.into(),
            BsxCol::CountTotal.create_series(count_total)?.into(),
            BsxCol::Density.create_series(density)?.into(),
        ]);

        Ok(unsafe { BsxBatch::new_unchecked(df) })
    }

    pub fn empty(chr_dtype: Option<&DataType>) -> Self {
        let df = DataFrame::from_iter([
            Series::new_empty(
                plsmallstr!(BsxCol::Chr.as_str()),
                chr_dtype.unwrap_or(&BsxCol::Chr.dtype()),
            ),
            create_empty_series!(Position),
            create_empty_series!(Strand),
            create_empty_series!(Context),
            create_empty_series!(CountM),
            create_empty_series!(CountTotal),
            create_empty_series!(Density),
        ]);
        unsafe { BsxBatch::new_unchecked(df) }
    }

    // CONVERSION
    pub fn data(&self) -> &DataFrame {
        &self.data
    }

    pub fn into_inner(self) -> DataFrame {
        self.data
    }

    pub fn column(
        &self,
        name: &str,
    ) -> Option<&Series> {
        if BsxCol::has_name(name) {
            Some(self.data().column(name).unwrap().as_materialized_series())
        }
        else {
            None
        }
    }

    // CHECKS
    pub fn is_chr_enum(&self) -> bool {
        self.data()
            .column(BsxCol::Chr.as_str())
            .unwrap()
            .dtype()
            .is_enum()
    }

    pub fn is_empty(&self) -> bool {
        self.data().is_empty()
    }

    // OPERATIONS
    pub fn split_at(
        &self,
        index: usize,
    ) -> (Self, Self) {
        let (left, right) = self.data().split_at(index as i64);
        #[allow(unsafe_code)]
        unsafe {
            (Self::new_unchecked(left), Self::new_unchecked(right))
        }
    }

    pub fn add_context_data(self, context_data: ContextData) -> Self {
        todo!()
    }

    #[allow(unsafe_code)]
    pub fn slice(
        &self,
        start: i64,
        length: usize,
    ) -> Self {
        let slice = self.data().slice(start, length);
        unsafe { Self::new_unchecked(slice) }
    }

    pub fn into_report(
        self,
        report_type: ReportType,
    ) -> PolarsResult<DataFrame> {
        let lf = self.data.lazy();
        let mut res = match report_type {
            ReportType::BedGraph => report_type_conversion::bedgraph(lf),
            ReportType::Bismark => report_type_conversion::bismark(lf),
            ReportType::CgMap => report_type_conversion::cgmap(lf),
            ReportType::Coverage => report_type_conversion::coverage(lf),
        };

        let res_schema = res.collect_schema()?;
        let target_schema = SchemaRef::new(report_type.schema());
        let schemas_equal = res_schema.deref().eq(&target_schema.deref());

        if !schemas_equal {
            Err(PolarsError::SchemaMismatch(
                format!("{:?} != {:?}", res_schema, target_schema).into(),
            ))
        }
        else {
            res.collect()
        }
    }

    pub fn get_methylation_stats(&self) -> MethylationStats {
        let nonull = self.density().drop_nulls();
        let mean = nonull.mean().unwrap_or(f64::NAN);
        let var = nonull.into_no_null_iter().map(|x| x as f64).variance();

        MethylationStats::from_data(
            mean,
            var,
            self.get_coverage_dist(),
            self.get_context_stats(),
            self.get_strand_stats(),
        )
    }

    pub fn get_coverage_dist(&self) -> HashMap<u16, u32> {
        izip!(self.count_m(), self.count_total())
            .filter_map(|(k, v)| Option::zip(k, v))
            .into_group_map()
            .into_iter()
            .map(|(k, v)| (k, v.iter().map(|x| *x as u32).sum()))
            .collect()
    }

    /// Returns context -> (sum methylation ratios, total counts)
    pub fn get_context_stats(&self) -> HashMap<Context, (f64, u32)> {
        izip!(self.context(), self.density())
            .filter_map(|(k, v)| {
                v.map(|density| (Context::from_bool(k), density as f64))
            })
            .into_group_map()
            .into_iter()
            .map(|(k, v)| (k, (v.iter().sum(), v.len() as u32)))
            .collect()
    }

    /// Returns strand -> (sum methylation ratios, total counts)
    pub fn get_strand_stats(&self) -> HashMap<Strand, (f64, u32)> {
        izip!(self.strand(), self.density())
            .filter_map(|(k, v)| {
                v.map(|density| (Strand::from_bool(k), density as f64))
            })
            .into_group_map()
            .into_iter()
            .map(|(k, v)| (k, (v.iter().sum(), v.len() as u32)))
            .collect()
    }

    pub fn as_binom(
        self,
        mean: f64,
        pvalue: f64,
    ) -> PolarsResult<Self> {
        let pvalue_vec = izip!(self.count_m(), self.count_total())
            .map(|(m, n)| (m.unwrap_or(0), n.unwrap_or(0)))
            .map(|(m, n)| {
                if n.is_zero() {
                    f64::NAN
                }
                else if m.is_zero() {
                    1.0
                }
                else if n == m {
                    0.0
                }
                else {
                    let binom =
                        statrs::distribution::Binomial::new(mean, n as u64).unwrap();
                    1.0 - binom.cdf(m as u64)
                }
            })
            .collect_vec();

        let binom_count_m: Column = Series::from_any_values(
            plsmallstr!(BsxCol::CountM.as_str()),
            &pvalue_vec
                .iter()
                .cloned()
                .map(|p| {
                    if 0.0 <= p && p <= pvalue {
                        1
                    }
                    else {
                        0
                    }
                })
                .map(|v| AnyValue::UInt16(v))
                .collect_vec(),
            true,
        )?
        .into();

        let binom_count_total: Column = Series::from_any_values(
            plsmallstr!(BsxCol::CountTotal.as_str()),
            &pvalue_vec
                .iter()
                .cloned()
                .map(|p| {
                    if !p.is_nan() {
                        1
                    }
                    else {
                        0
                    }
                })
                .map(|v| AnyValue::UInt16(v))
                .collect_vec(),
            true,
        )?
        .into();

        let density = (&binom_count_m.cast(&DataType::Float32)?
            / &binom_count_total.cast(&DataType::Float32)?)?
            .with_name(plsmallstr!(BsxCol::Density.as_str()));

        let mut new_data = self.data.clone();
        new_data.with_column(binom_count_m)?;
        new_data.with_column(binom_count_total)?;
        new_data.with_column(density)?;

        Ok(unsafe { Self::new_unchecked(new_data) })
    }

    // POSITION
    #[inline]
    pub fn first_pos(&self) -> Option<u32> {
        self.position().first()
    }

    #[inline]
    pub fn last_pos(&self) -> Option<u32> {
        self.position().last()
    }

    pub fn seqname(&self) -> Option<&str> {
        self.chr().iter_str().next().flatten()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.data().height()
    }

    pub fn first_genomic_pos(&self) -> Option<GenomicPosition<BsxSmallStr, u32>> {
        let seqname = self.seqname().map(BsxSmallStr::from);
        let pos = self.first_pos();
        seqname.and_then(|seqname| pos.map(|pos| GenomicPosition::new(seqname, pos)))
    }

    pub fn last_genomic_pos(&self) -> Option<GenomicPosition<BsxSmallStr, u32>> {
        let seqname = self.seqname().map(BsxSmallStr::from);
        let pos = self.last_pos();
        seqname.and_then(|seqname| pos.map(|pos| GenomicPosition::new(seqname, pos)))
    }

    pub fn as_contig(&self) -> Option<Contig<BsxSmallStr, u32>> {
        let seqname = self.seqname().map(BsxSmallStr::from);
        let first = self.first_pos();
        let last = self.last_pos();
        if let (Some(s), Some(f), Some(l)) = (seqname, first, last) {
            Some(Contig::new(s, f, l, Strand::None))
        }
        else {
            None
        }
    }
}

mod report_type_conversion {
    use super::*;

    pub fn bedgraph(lf: LazyFrame) -> LazyFrame {
        lf.select([
            // chr
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            BsxCol::Position.col().alias("start"),
            BsxCol::Position.col().alias("end"),
            BsxCol::Density.col().alias("density"),
        ])
        .drop_nans(Some(vec![BsxCol::Density.col()]))
    }

    pub fn coverage(lf: LazyFrame) -> LazyFrame {
        lf.select([
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            BsxCol::Position.col().alias("start"),
            BsxCol::Position.col().alias("end"),
            BsxCol::CountM.col().alias("count_m"),
            (BsxCol::CountTotal.col() - BsxCol::CountM.col()).alias("count_um"),
            BsxCol::Density.col().alias("density"),
        ])
        .drop_nans(Some(vec![BsxCol::Density.col()]))
    }

    pub fn bismark(lf: LazyFrame) -> LazyFrame {
        lf.with_column(
            when(BsxCol::Context.col() == lit(NULL))
                .then(lit("CHH"))
                .when(BsxCol::Context.col() == lit(false))
                .then(lit("CHG"))
                .otherwise("CG")
                .cast(DataType::String)
                .alias(BsxCol::Context.as_str()),
        )
        .select([
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            BsxCol::Position.col().alias("position"),
            // strand
            when(BsxCol::Strand.col() == lit(true))
                .then(lit("+"))
                .when(BsxCol::Strand.col() == lit(false))
                .then(lit("-"))
                .otherwise(".")
                .cast(DataType::String)
                .alias("strand"),
            BsxCol::CountM.col().alias("count_m"),
            (BsxCol::CountTotal.col() - BsxCol::CountM.col()).alias("count_um"),
            BsxCol::Density.col().alias("density"),
            BsxCol::Context.col().alias("context"),
            BsxCol::Context.col().alias("trinuc"),
        ])
        .drop_nans(Some(vec![BsxCol::Density.col()]))
    }

    pub fn cgmap(lf: LazyFrame) -> LazyFrame {
        lf.with_column(
            when(BsxCol::Context.col() == lit(NULL))
                .then(lit("CHH"))
                .when(BsxCol::Context.col() == lit(false))
                .then(lit("CHG"))
                .otherwise("CG")
                .cast(DataType::String)
                .alias(BsxCol::Context.as_str()),
        )
        .select([
            // chr
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            // nuc
            when(BsxCol::Strand.col() == lit(true))
                .then(lit("C"))
                .when(BsxCol::Strand.col() == lit(false))
                .then(lit("G"))
                .otherwise(lit("."))
                .alias("nuc"),
            // position
            BsxCol::Position.col().alias("position"),
            // context
            BsxCol::Context.col().alias("context"),
            // dinuc
            BsxCol::Context.col().str().head(lit(2)).alias("dinuc"),
            // density
            BsxCol::Density.col().alias("density"),
            // count_m
            BsxCol::CountM.col().alias("count_m"),
            // count_total
            (BsxCol::CountTotal.col() - BsxCol::CountM.col()).alias("count_um"),
        ])
        .drop_nans(Some(vec![BsxCol::Density.col()]))
    }
}

struct LazyBsxBatch {
    data: LazyFrame,
}

impl LazyBsxBatch {
    fn filter(
        self,
        predicate: Expr,
    ) -> Self {
        Self {
            data:     self.data.filter(predicate),
        }
    }

    /// Creates a LazyBsxBatch from an existing LazyFrame.
    fn from_lazy(lazy: LazyFrame) -> Self {
        Self {
            data:     lazy,
        }
    }


    /// Filters positions less than the specified value.
    pub fn filter_pos_lt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(BsxCol::Position.col().lt(lit(value)))
    }

    /// Filters positions greater than the specified value.
    pub fn filter_pos_gt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(BsxCol::Position.col().gt(lit(value)))
    }

    /// Filters entries with coverage less than the specified value.
    pub fn filter_coverage_lt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(BsxCol::CountTotal.col().lt(lit(value)))
    }

    /// Filters entries by strand value.
    pub fn filter_strand(
        self,
        value: Strand,
    ) -> Self {
        self.filter(BsxCol::Strand.col().eq(value
            .to_bool()
            .map(lit)
            .unwrap_or(lit(NULL))))
    }

    /// Filters entries by context value.
    pub fn filter_context(
        self,
        value: Context,
    ) -> Self {
        self.filter(BsxCol::Context.col().eq(value
            .to_bool()
            .map(lit)
            .unwrap_or(lit(NULL))))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::utils::get_categorical_dtype;

    #[test]
    fn test_empty_batch() {
        let batch = BsxBatch::empty(None);
        assert!(batch.is_empty());
    }

    #[rstest]
    #[case::no_chr(None, None)]
    #[case::both_chr(Some(get_categorical_dtype(vec!["chr1".into()])), Some(get_categorical_dtype(vec!["chr1".into()])))]
    #[should_panic]
    #[case::different_types(None, Some(get_categorical_dtype(vec!["chr1".into()])))]
    fn test_can_extend(
        #[case] first_dtype: Option<DataType>,
        #[case] second_dtype: Option<DataType>,
    ) {
        let batch1 = BsxBatch::empty(first_dtype.as_ref());
        let batch2 = BsxBatch::try_from_columns(
            "chr1",
            second_dtype,
            vec![1, 2, 3],
            vec![true, false, true],
            vec![Some(true), Some(false), None],
            vec![1, 2, 3],
            vec![3, 6, 9],
            vec![0.3, 0.3, 0.3],
        )
        .unwrap();

        assert!(matches!(
            batch1.column(BsxCol::Chr.as_str()).unwrap().dtype(),
            DataType::Categorical(_, _) | DataType::Enum(_, _)
        ));
        assert!(matches!(
            batch2.column(BsxCol::Chr.as_str()).unwrap().dtype(),
            DataType::Categorical(_, _) | DataType::Enum(_, _)
        ));

        let vstack = batch1.data().vstack(&batch2.data()).unwrap();
        assert!(matches!(
            vstack.column(BsxCol::Chr.as_str()).unwrap().dtype(),
            DataType::Categorical(_, _) | DataType::Enum(_, _)
        ));
        assert!(vstack
            .column(BsxCol::Chr.as_str())
            .unwrap()
            .categorical()
            .unwrap()
            .get_rev_map()
            .find("chr1")
            .is_some());
    }
}
