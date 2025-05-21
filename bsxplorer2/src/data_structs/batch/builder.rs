use hashbrown::HashMap;
use itertools::Itertools;
use polars::prelude::*;

use super::{create_caregorical_dtype,
            create_empty_categorical_dtype,
            name_dtype_tuple,
            BsxBatch,
            BsxColumns as BsxCol};
use crate::data_structs::coords::Contig;
use crate::io::report::ReportType;
use crate::with_field_fn;

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
            chr_values.into_iter().map(Some).collect_vec(),
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
                        BsxCol::Chr.as_str(),
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
            .select(BsxCol::colnames().into_iter().map(col).collect_vec())
            .collect()
    }

    pub fn build_from_report(
        &self,
        df: DataFrame,
        report_type: ReportType,
    ) -> PolarsResult<BsxBatch> {
        let lf = df.lazy();
        let mut converted = match report_type {
            ReportType::Bismark => build::bismark_expr(lf),
            ReportType::CgMap => build::cgmap_expr(lf),
            ReportType::BedGraph => build::bedgraph_expr(lf),
            ReportType::Coverage => build::coverage_expr(lf),
        }
        .collect()?;

        if self.rechunk {
            converted.rechunk_mut()
        }
        self.checks_only(&converted)?;
        let casted = self.cast_only(converted)?;

        Ok(unsafe { BsxBatch::new_unchecked(casted) })
    }

    /// # Safety
    /// Batch data and schema are assumed to be correct
    pub unsafe fn build_unchecked(df: DataFrame) -> BsxBatch {
        BsxBatch::new_unchecked(df)
    }

    pub fn concat(batches: Vec<BsxBatch>) -> PolarsResult<BsxBatch> {
        let contigs = batches
            .iter()
            .filter(|b| !b.is_empty())
            .map(|b| b.as_contig())
            .collect::<Option<Vec<_>>>();

        if contigs.is_none() {
            return Ok(BsxBatch::empty(None));
        }
        let contigs = contigs.unwrap();

        if !contigs.iter().map(Contig::seqname).all_equal() {
            return Err(PolarsError::ComputeError("Chromosomes do not match".into()));
        }
        if !contigs.windows(2).all(|w| w[1].start() >= w[0].end()) {
            return Err(PolarsError::ComputeError(
                "Batch positions are not sorted".into(),
            ));
        }
        let data = concat(
            batches
                .into_iter()
                .map(|b| b.into_inner().lazy())
                .collect_vec(),
            Default::default(),
        )?
        .collect()?;

        let builder = BsxBatchBuilder::no_checks()
            .with_rechunk(true)
            .with_check_duplicates(true)
            .with_chr_dtype(Some(
                data.schema().get(BsxCol::Chr.as_str()).unwrap().clone(),
            ));

        builder.checks_only(&data)?;
        Ok(unsafe { BsxBatch::new_unchecked(builder.cast_only(data)?) })
    }
}

pub(super) mod build {
    use polars::series::IsSorted;

    use super::*;

    /// Check duplicates in the position column
    pub fn check_no_duplicates(df: &DataFrame) -> PolarsResult<bool> {
        Ok(df.column(BsxCol::Position.as_str())?.n_unique()? == df.height())
    }

    /// Check sortedness of the position column
    pub fn check_sorted(df: &DataFrame) -> PolarsResult<bool> {
        let pos_col = df.column(BsxCol::Position.as_str())?;
        if pos_col.n_chunks() > 1 {
            return Err(PolarsError::ComputeError(
                "Column is chunked. Rechunk before operation".into(),
            ));
        }

        Ok(df
            .column(BsxCol::Position.as_str())?
            .as_materialized_series()
            .iter()
            .map(|v| v.try_extract())
            .collect::<PolarsResult<Vec<u64>>>()?
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
        lf.with_columns([
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
}

#[cfg(test)]
mod tests {
    use polars::series::IsSorted;
    use rstest::{fixture, rstest};

    use super::build::*;
    use super::*;

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
        #[case] check_fn: fn(&DataFrame) -> PolarsResult<bool>,
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
        let res = df
            .lazy()
            .with_column(nuc_to_bool_expr().alias("result"))
            .collect()
            .unwrap();
        let results = res.column("result").unwrap().bool().unwrap();
        assert_eq!(results.get(0), Some(true));
        assert_eq!(results.get(1), Some(false));
        assert_eq!(results.get(2), None);

        // Test strand_to_bool_expr
        let df = df!("strand" => ["+", "-", "."]).unwrap();
        let res = df
            .lazy()
            .with_column(strand_to_bool_expr().alias("result"))
            .collect()
            .unwrap();
        let results = res.column("result").unwrap().bool().unwrap();
        assert_eq!(results.get(0), Some(true));
        assert_eq!(results.get(1), Some(false));
        assert_eq!(results.get(2), None);

        // Test context_to_bool_expr
        let df = df!("context" => ["CG", "CHG", "CHH"]).unwrap();
        let res = df
            .lazy()
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
        )
        .unwrap();

        // Test count_total_col_expr
        let res = df
            .clone()
            .lazy()
            .with_column(count_total_col_expr().alias("count_total"))
            .collect()
            .unwrap();
        let totals = res.column("count_total").unwrap().i64().unwrap();
        assert_eq!(totals.get(0), Some(10));
        assert_eq!(totals.get(1), Some(20));

        // Test density_col_expr
        let res = df
            .lazy()
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
        let batch =
            BsxBatchBuilder::all_checks().build_from_report(input_df, report_type)?;
        let _data = batch.data();
        Ok(())
    }
}
