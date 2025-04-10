use super::encoded::EncodedBsxBatch;
use crate::data_structs::batch::decoded::BsxBatch;
use crate::data_structs::batch::traits::{BsxBatchMethods, BsxColNames};
use crate::io::report::schema::ReportTypeSchema;
use crate::utils::{decode_context, decode_strand, encode_context, encode_strand};
use anyhow::anyhow;
use log::warn;
use polars::prelude::*;
use polars::series::IsSorted;

/// Builder for constructing and validating BSX batch data
#[derive(Debug, Copy, Clone)]
pub struct BsxBatchBuilder {
    report_type: Option<ReportTypeSchema>,
    check_nulls: bool,
    check_sorted: bool,
    check_duplicates: bool,
    rechunk: bool,
    check_single_chr: bool,
}

impl Default for BsxBatchBuilder {
    fn default() -> Self {
        BsxBatchBuilder {
            report_type: None,
            check_nulls: true,
            check_sorted: true,
            check_duplicates: true,
            rechunk: true,
            check_single_chr: true
        }
    }
}

// Public methods
impl BsxBatchBuilder {
    /// Creates a builder with all data validation checks enabled
    pub fn all_checks() -> Self {
        Self {
            report_type: None,
            check_duplicates: true,
            check_sorted: true,
            rechunk: true,
            check_nulls: true,
            check_single_chr: true
        }
    }

    /// Creates a builder with all validation checks disabled
    pub fn no_checks() -> Self {
        Self {
            report_type: None,
            check_duplicates: false,
            check_sorted: false,
            rechunk: false,
            check_nulls: false,
            check_single_chr: false
        }
    }

    /// Sets the report type schema for data conversion
    pub fn with_report_type(mut self, report_type: ReportTypeSchema) -> Self {
        self.report_type = Some(report_type);
        self
    }

    /// Sets whether to check for null values in critical columns
    pub fn with_check_nulls(mut self, check_nulls: bool) -> Self {
        self.check_nulls = check_nulls;
        self
    }

    /// Sets whether to check and ensure positions are sorted
    pub fn with_check_sorted(mut self, check_sorted: bool) -> Self {
        self.check_sorted = check_sorted;
        self
    }

    /// Sets whether to rechunk the data for memory efficiency
    pub fn with_rechunk(mut self, rechunk: bool) -> Self {
        self.rechunk = rechunk;
        self
    }

    /// Sets whether to check if all data is from a single chromosome
    pub fn with_check_single_chr(mut self, check: bool) -> Self {
        self.check_single_chr = check;
        self
    }

    /// Sets whether to check for duplicate positions
    pub fn with_check_duplicates(mut self, check_duplicates: bool) -> Self {
        self.check_duplicates = check_duplicates;
        self
    }

    /// Builds a decoded BSX batch from a DataFrame
    pub fn build_decoded(self, data: DataFrame) -> anyhow::Result<BsxBatch> {
        let casted = self.cast_decoded(data)?;
        self.run_checks(&casted)?;
        let casted = self.sort(casted)?;
        let casted = self.rechunk(casted);

        Ok(unsafe { BsxBatch::new_unchecked(casted) })
    }

    /// Builds an encoded BSX batch from a DataFrame
    pub fn build_encoded(self, data: DataFrame, chr_dtype: DataType) -> anyhow::Result<EncodedBsxBatch> {
        let casted = self.cast_encoded(data, chr_dtype)?;
        self.run_checks(&casted)?;
        let casted = self.sort(casted)?;
        let casted = self.rechunk(casted);

        Ok(unsafe { EncodedBsxBatch::new_unchecked(casted) })
    }

    /// Validates and optimizes a DataFrame according to builder settings
    pub fn check(&self, data: DataFrame) -> anyhow::Result<DataFrame> {
        self.run_checks(&data)?;
        let sorted = self.sort(data)?;
        let res = self.rechunk(sorted);
        Ok(res)
    }

    /// Converts a decoded batch to an encoded batch format
    pub fn encode_batch(batch: BsxBatch, chr_dtype: DataType) -> anyhow::Result<EncodedBsxBatch> {
        let chr = batch.chr_val()?.to_string();
        let start = batch.start_pos().ok_or("no data").unwrap();
        let end = batch.end_pos().ok_or("no data").unwrap();
        let batch_data = DataFrame::from(batch);
        let mut batch_lazy = batch_data.lazy();

        batch_lazy = encode_context(batch_lazy, "context");
        batch_lazy = encode_strand(batch_lazy, "strand");

        let casted = encoded::select_cast(batch_lazy, chr_dtype).collect()?;
        let encoded = EncodedBsxBatch::new_with_fields(
            casted, chr
        );
        Ok(encoded)
    }

    /// Converts an encoded batch to a decoded batch format
    pub fn decode_batch(batch: EncodedBsxBatch) -> anyhow::Result<BsxBatch> {
        let batch_data = DataFrame::from(batch);
        let mut batch_lazy = batch_data.lazy();
        batch_lazy = decode_context(batch_lazy, "context", "context");
        batch_lazy = decode_strand(batch_lazy, "strand", "strand");
        let casted = decoded::select_cast(batch_lazy);

        let result = casted.collect()?;
        Ok( unsafe { BsxBatch::new_unchecked(result) } )
    }
}

// Private methods
impl BsxBatchBuilder {
    /// Casts data to decoded format using appropriate schema conversion
    fn cast_decoded(&self, data: DataFrame) -> anyhow::Result<DataFrame> {
        let casted = if let Some(report_type) = self.report_type {
            match report_type {
                ReportTypeSchema::Bismark => decoded::from_bismark(data)?,
                ReportTypeSchema::CgMap => decoded::from_cgmap(data)?,
                ReportTypeSchema::BedGraph => decoded::from_bedgraph(data)?,
                ReportTypeSchema::Coverage => decoded::from_coverage(data)?,
            }
        } else {
            decoded::select_cast(data.lazy()).collect()?
        };
        Ok(casted)
    }

    /// Casts data to encoded format using appropriate schema conversion
    fn cast_encoded(&self, data: DataFrame, chr_dtype: DataType) -> anyhow::Result<DataFrame> {
        let casted = if let Some(report_type) = self.report_type {
            match report_type {
                ReportTypeSchema::Bismark => encoded::from_bismark(data, chr_dtype)?,
                ReportTypeSchema::CgMap => encoded::from_cgmap(data, chr_dtype)?,
                ReportTypeSchema::BedGraph => encoded::from_bedgraph(data, chr_dtype)?,
                ReportTypeSchema::Coverage => encoded::from_coverage(data, chr_dtype)?,
            }
        } else {
            encoded::select_cast(data.lazy(), chr_dtype).collect()?
        };
        Ok(casted)
    }

    /// Rechunks the data if enabled in the builder
    fn rechunk(&self, data: DataFrame) -> DataFrame {
        let mut data = data;
        if self.rechunk {
            data.rechunk_mut()
        }
        data
    }

    /// Performs data validation based on builder settings
    fn run_checks(&self, data: &DataFrame) -> anyhow::Result<()> {
        if self.check_nulls {
            check_has_nulls(data)?;
        }
        if self.check_duplicates && !(data.column(BsxBatch::POS_NAME)?.n_unique()? != data.height()) {
            anyhow::bail!("Duplicated positions")
        }
        Ok(())
    }

    /// Sorts data by position if needed and updates sorting flags
    fn sort(&self, data: DataFrame) -> anyhow::Result<DataFrame> {
        let mut data = data;
        if self.check_sorted && !check_pos_ascending(&data) {
            warn!("Position is not sorted");
            data = data.sort(
                ["position"],
                SortMultipleOptions::default().with_order_descending(false)
            )?;
        } else {
            let pos_col_idx = data.get_column_index(BsxBatch::POS_NAME).unwrap();
            unsafe {
                let columns = data.get_columns_mut();
                let pos_col = columns.get_mut(pos_col_idx).unwrap();
                pos_col.set_sorted_flag(IsSorted::Ascending);
            }
        }
        Ok(data)
    }
}


mod decoded {
    use super::*;

    /// Selects and casts columns to the proper types for BSX batch processing
    ///
    /// Transforms input data to a standardized schema with properly typed columns
    pub(crate) fn select_cast(lf: LazyFrame) -> LazyFrame {
        lf.select([
            col("chr").cast(BsxBatch::chr_type()).alias(BsxBatch::CHR_NAME),
            col("position").cast(BsxBatch::pos_type()).alias(BsxBatch::POS_NAME),
            col("strand").cast(BsxBatch::strand_type()).alias(BsxBatch::STRAND_NAME),
            col("context").cast(BsxBatch::context_type()).alias(BsxBatch::CONTEXT_NAME),
            col("count_m").cast(BsxBatch::count_type()).alias(BsxBatch::COUNT_M_NAME),
            col("count_total").cast(BsxBatch::count_type()).alias(BsxBatch::COUNT_TOTAL_NAME),
            col("density").cast(BsxBatch::density_type()).alias(BsxBatch::DENSITY_NAME),
        ])
    }

    /// Creates an expression to convert nucleotide values to strand symbols
    fn nuc_to_strand_expr() -> Expr {
        when(col("nuc").eq(lit("C")))
            .then(lit("+"))
            .when(col("nuc").eq(lit("G")))
            .then(lit("-"))
            .otherwise(lit("."))
            .alias("strand")
    }

    /// Processes Bismark format data into standardized BSX format
    pub(crate) fn from_bismark(df: DataFrame) -> anyhow::Result<DataFrame> {
        let col_added = df.lazy()
            .with_column(count_total_col_expr())
            .with_column(density_col_expr());
        let casted = select_cast(col_added);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Converts CG map format data into standardized BSX format
    pub(crate) fn from_cgmap(df: DataFrame) -> anyhow::Result<DataFrame> {
        let col_added = df.lazy()
            .with_column(nuc_to_strand_expr());
        let casted = select_cast(col_added);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Processes coverage format data into standardized BSX format
    pub(crate) fn from_coverage(df: DataFrame) -> anyhow::Result<DataFrame> {
        let col_added = df.lazy()
            .with_column(count_total_col_expr())
            .with_columns([
                col("start").alias("position"),
                lit(".").alias("strand"),
                lit(NULL).alias("context"),
                density_col_expr()
            ]);
        let casted = select_cast(col_added);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Processes BedGraph format data into standardized BSX format
    pub(crate) fn from_bedgraph(df: DataFrame) -> anyhow::Result<DataFrame> {
        let col_added = df.lazy()
            .with_columns([
                col("start").alias("position"),
                lit(".").alias("strand"),
                lit(NULL).alias("context"),
                lit(NULL).alias("count_m"),
                lit(NULL).alias("count_total")
            ]);
        let casted = select_cast(col_added);
        casted.collect().map_err(|e| anyhow!(e))
    }

}

mod encoded {
    use super::*;
    use crate::data_structs::batch::encoded::EncodedBsxBatch;

    /// Selects and casts columns to the proper types for BSX batch processing
    ///
    /// Transforms input data to a standardized schema with properly typed columns
    pub(crate) fn select_cast(lf: LazyFrame, chr_dtype: DataType) -> LazyFrame {
        lf.select([
            col("chr").cast(chr_dtype).alias(EncodedBsxBatch::CHR_NAME),
            col("position").cast(EncodedBsxBatch::pos_type()).alias(EncodedBsxBatch::POS_NAME),
            col("strand").cast(EncodedBsxBatch::strand_type()).alias(EncodedBsxBatch::STRAND_NAME),
            col("context").cast(EncodedBsxBatch::context_type()).alias(EncodedBsxBatch::CONTEXT_NAME),
            col("count_m").cast(EncodedBsxBatch::count_type()).alias(EncodedBsxBatch::COUNT_M_NAME),
            col("count_total").cast(EncodedBsxBatch::count_type()).alias(EncodedBsxBatch::COUNT_TOTAL_NAME),
            col("density").cast(EncodedBsxBatch::density_type()).alias(EncodedBsxBatch::DENSITY_NAME),
        ])
    }

    fn nuc_to_strand_expr() -> Expr {
        when(col("nuc").eq(lit("C")))
            .then(lit(true))
            .when(col("nuc").eq(lit("G")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .alias("strand")
    }

    /// Processes Bismark format data into standardized BSX format
    pub(crate) fn from_bismark(df: DataFrame, chr_dtype: DataType) -> anyhow::Result<DataFrame> {
        let col_added = df.lazy()
            .with_column(count_total_col_expr())
            .with_column(density_col_expr());
        let casted = select_cast(col_added, chr_dtype);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Converts CG map format data into standardized BSX format
    pub(crate) fn from_cgmap(df: DataFrame, chr_dtype: DataType) -> anyhow::Result<DataFrame> {
        let col_added = df.lazy()
            .with_column(nuc_to_strand_expr());
        let casted = select_cast(col_added, chr_dtype);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Processes coverage format data into standardized BSX format
    pub(crate) fn from_coverage(df: DataFrame, chr_dtype: DataType) -> anyhow::Result<DataFrame> {
        let col_added = df.lazy()
            .with_column(count_total_col_expr())
            .with_columns([
                col("start").alias("position"),
                lit(NULL).alias("strand"),
                lit(NULL).alias("context"),
                density_col_expr()
            ]);
        let casted = select_cast(col_added, chr_dtype);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Processes BedGraph format data into standardized BSX format
    pub(crate) fn from_bedgraph(df: DataFrame, chr_dtype: DataType) -> anyhow::Result<DataFrame> {
        let col_added = df.lazy()
            .with_columns([
                col("start").alias("position"),
                lit(NULL).alias("strand"),
                lit(NULL).alias("context"),
                lit(NULL).alias("count_m"),
                lit(NULL).alias("count_total")
            ]);
        let casted = select_cast(col_added, chr_dtype);
        casted.collect().map_err(|e| anyhow!(e))
    }
}

/// Validates that critical columns do not contain null values
fn check_has_nulls(df: &DataFrame) -> anyhow::Result<()> {
    if df.column("chr").expect("chr column not found").null_count() > 0 {
        anyhow::bail!("Nulls not allowed in 'chr' column");
    }
    if df.column("position").expect("position column not found").null_count() > 0 {
        anyhow::bail!("Nulls not allowed in 'position' column");
    }
    Ok(())
}

/// Checks if position column values are in ascending order
fn check_pos_ascending(df: &DataFrame) -> bool {
    let pos = df.column("position").expect("position column not found");
    pos.as_series().unwrap().iter().is_sorted()
}


/// Creates an expression to compute the total count from methylated and unmethylated counts
fn count_total_col_expr() -> Expr {
    (col("count_m") + col("count_um")).alias("count_total")
}

/// Creates an expression to compute methylation density from count data
fn density_col_expr() -> Expr {
    (col("count_m") / col("count_total")).alias("density")
}
