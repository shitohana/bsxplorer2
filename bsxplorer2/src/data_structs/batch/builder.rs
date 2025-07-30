use hashbrown::HashMap;
use itertools::Itertools;
use polars::prelude::*;

use super::{
    create_caregorical_dtype,
    create_empty_categorical_dtype,
    name_dtype_tuple,
    BsxBatch,
    BsxColumns as BsxCol,
};
use crate::data_structs::coords::Contig;
use crate::io::report::ReportType;
use crate::with_field_fn;

/// Builder for creating `BsxBatch` instances with various configuration
/// options.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BsxBatchBuilder {
    /// Optional categorical dtype for the chromosome column.
    pub chr_dtype:        Option<DataType>,
    /// Flag to enable null checks.
    pub check_nulls:      bool,
    /// Flag to enable sorted checks on the position column.
    pub check_sorted:     bool,
    /// Flag to enable duplicate checks on the position column.
    pub check_duplicates: bool,
    /// Flag to rechunk the resulting DataFrame.
    pub rechunk:          bool,
    /// Flag to check for a single chromosome value per batch.
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

    /// Sets the chromosome values and configures the categorical dtype.
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

    /// Performs validation checks on the DataFrame based on builder
    /// configuration.
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

    /// Casts the DataFrame columns to the expected `BsxBatch` dtypes.
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

    /// Builds a `BsxBatch` from a DataFrame based on a specified report type.
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

    /// Concatenates multiple `BsxBatch` instances into a single batch.
    pub fn concat(mut batches: Vec<BsxBatch>) -> PolarsResult<BsxBatch> {
        if batches.is_empty() {
            return Err(PolarsError::NoData("Vector is empty".into()));
        }
        else if batches.len() == 1 {
            return Ok(batches.pop().unwrap());
        }

        let contigs = batches
            .iter()
            .filter(|b| !b.is_empty())
            .map(|b| b.as_contig())
            .collect::<Option<Vec<_>>>()
            .ok_or(PolarsError::NoData("All batches are empty".into()))?;

        if !contigs.iter().map(Contig::seqname).all_equal() {
            return Err(PolarsError::ComputeError("Chromosomes do not match".into()));
        }
        if !contigs.windows(2).all(|w| w[1].start() >= w[0].end()) {
            return Err(PolarsError::ComputeError(
                "Batch positions are not sorted".into(),
            ));
        }

        let mut batches_data = unsafe {
            batches
                .into_iter()
                .map(BsxBatch::into_inner)
                .reduce(|mut acc, df| {
                    acc.get_columns_mut()
                        .iter_mut()
                        .zip(df.get_columns().iter())
                        .for_each(|(left, right)| {
                            left.append(right).expect("should not fail");
                        });
                    acc.set_height(acc.height() + df.height());
                    acc
                })
                .unwrap_unchecked()
        };
        batches_data.align_chunks_par();
        if batches_data.should_rechunk() {
            batches_data.rechunk_mut();
        };

        Ok(unsafe { BsxBatch::new_unchecked(batches_data) })
    }
}

/// Internal helper functions for building and validating `BsxBatch` DataFrames.
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

    /// Sets the sorted flag on the position column.
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

    /// Creates an expression to convert strand symbols (+/-) to boolean.
    pub fn strand_to_bool_expr() -> Expr {
        when(col("strand").eq(lit("+")))
            .then(lit(true))
            .when(col("strand").eq(lit("-")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
    }

    /// Creates an expression to convert context values (CG/CHG) to boolean.
    pub fn context_to_bool_expr() -> Expr {
        when(col("context").eq(lit("CG")))
            .then(lit(true))
            .when(col("context").eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
    }

    /// Creates an expression to calculate total counts (m + um).
    pub fn count_total_col_expr() -> Expr {
        col("count_m") + col("count_um")
    }

    /// Creates an expression to calculate methylation density (m / total).
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
            context_to_bool_expr().alias(BsxCol::Context.as_str()),
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
