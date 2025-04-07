use crate::data_structs::bsx_batch::{BsxBatch, BsxBatchMethods};
use crate::io::report::schema::ReportTypeSchema;
use anyhow::anyhow;
use itertools::Itertools;
use log::warn;
use polars::prelude::*;

struct BsxBatchBuilder {
    report_type: Option<ReportTypeSchema>,
    batch_size: usize,
    check_nulls: bool,
    check_sorted: bool,
    rechunk: bool,
    partition: bool,
}


impl Default for BsxBatchBuilder {
    fn default() -> Self {
        BsxBatchBuilder {
            report_type: None,
            batch_size: 10_000,
            check_sorted: true,
            check_nulls: true,
            rechunk: true,
            partition: true,
        }
    }
}

/// Builds a standardized BSX batch from various input formats
///
/// Splits input data by chromosomes. Processes input data according to 
/// specified report type, performs validation, and ensures proper sorting of 
/// genomic positions
pub fn build_bsx_batch(
    data: DataFrame,
    report_type: Option<&ReportTypeSchema>,
    check_nulls: bool,
    check_sorted: bool,
    rechunk: bool,
) -> anyhow::Result<Vec<BsxBatch>> {
    let mut casted = if let Some(report_type) = report_type {
        match report_type {
            ReportTypeSchema::Bismark => from_bismark(data)?,
            ReportTypeSchema::CgMap => from_cgmap(data)?,
            ReportTypeSchema::BedGraph => from_bedgraph(data)?,
            ReportTypeSchema::Coverage => from_coverage(data)?,
        }
    } else {
        select_cast_bsx(data.lazy()).collect()?
    };

    if rechunk {
        casted.rechunk_mut()
    }

    let mut res = Vec::new();
    for batch in casted.partition_by_stable(["chr"], true)? {
        if check_nulls {
            check_has_nulls(&casted)?;
        }

        if check_sorted && !check_pos_ascending(&casted) {
            warn!("Position is not sorted");
            casted = casted.sort(
                ["position"],
                SortMultipleOptions::default().with_order_descending(false)
            )?;
        }

        res.push(unsafe { BsxBatch::new_unchecked(batch) });
    }
    Ok(res)
}

/// Selects and casts columns to the proper types for BSX batch processing
///
/// Transforms input data to a standardized schema with properly typed columns
fn select_cast_bsx(lf: LazyFrame) -> LazyFrame {
    lf.select([
        col("chr").cast(BsxBatch::CHR_DTYPE).alias(BsxBatch::CHR_NAME),
        col("position").cast(BsxBatch::POS_DTYPE).alias(BsxBatch::POS_NAME),
        col("strand").cast(BsxBatch::STRAND_DTYPE).alias(BsxBatch::STRAND_NAME),
        col("context").cast(BsxBatch::CONTEXT_DTYPE).alias(BsxBatch::CONTEXT_NAME),
        col("count_m").cast(BsxBatch::COUNT_DTYPE).alias(BsxBatch::COUNT_M_NAME),
        col("count_total").cast(BsxBatch::COUNT_DTYPE).alias(BsxBatch::COUNT_TOTAL_NAME),
        col("density").cast(BsxBatch::DENSITY_DTYPE).alias(BsxBatch::DENSITY_NAME),
    ])
}

/// Creates an expression to compute the total count from methylated and unmethylated counts
fn count_total_col_expr() -> Expr {
    (col("count_m") + col("count_um")).alias("count_total")
}

/// Creates an expression to compute methylation density from count data
fn density_col_expr() -> Expr {
    (col("count_m") / col("count_total")).alias("density")
}

/// Creates an expression to convert nucleotide information to strand designation
fn nuc_to_strand_expr() -> Expr {
    when(col("nuc").eq(lit("C")))
        .then(lit("+"))
        .when(col("nuc").eq(lit("G")))
        .then(lit("-"))
        .otherwise(lit("."))
        .alias("strand")
}

/// Processes Bismark format data into standardized BSX format
fn from_bismark(df: DataFrame) -> anyhow::Result<DataFrame> {
    let col_added = df.lazy()
        .with_column(count_total_col_expr())
        .with_column(density_col_expr());
    let casted = select_cast_bsx(col_added);
    casted.collect().map_err(|e| anyhow!(e))
}

/// Converts CG map format data into standardized BSX format
fn from_cgmap(df: DataFrame) -> anyhow::Result<DataFrame> {
    let col_added = df.lazy()
        .with_column(nuc_to_strand_expr());
    let casted = select_cast_bsx(col_added);
    casted.collect().map_err(|e| anyhow!(e))
}

/// Processes coverage format data into standardized BSX format
fn from_coverage(df: DataFrame) -> anyhow::Result<DataFrame> {
    let col_added = df.lazy()
        .with_column(count_total_col_expr())
        .with_columns([
            col("start").alias("position"),
            lit(".").alias("strand"),
            lit(NULL).alias("context"),
            density_col_expr()
        ]);
    let casted = select_cast_bsx(col_added);
    casted.collect().map_err(|e| anyhow!(e))
}

/// Processes BedGraph format data into standardized BSX format
fn from_bedgraph(df: DataFrame) -> anyhow::Result<DataFrame> {
    let col_added = df.lazy()
        .with_columns([
            col("start").alias("position"),
            lit(".").alias("strand"),
            lit(NULL).alias("context"),
            lit(NULL).alias("count_m"),
            lit(NULL).alias("count_total")
        ]);
    let casted = select_cast_bsx(col_added);
    casted.collect().map_err(|e| anyhow!(e))
}

/// Validates that critical columns do not contain null values
fn check_has_nulls(df: &DataFrame) -> anyhow::Result<()> {
    if df.column("chr").expect("chr column not found").null_count() > 0 {
        return anyhow::bail!("Nulls not allowed in 'chr' column");
    }
    if df.column("position").expect("position column not found").null_count() > 0 {
        return anyhow::bail!("Nulls not allowed in 'position' column");
    }
    Ok(())
}

/// Checks if position column values are in ascending order
fn check_pos_ascending(df: &DataFrame) -> bool {
    let pos = df.column("position").expect("position column not found");
    pos.as_series().unwrap().iter().is_sorted()
}
