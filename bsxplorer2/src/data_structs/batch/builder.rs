use super::colnames::*;
use super::encoded::EncodedBsxBatch;
use crate::data_structs::batch::decoded::BsxBatch;
use crate::data_structs::batch::traits::{
    BatchType, BsxBatchMethods, BsxTypeTag,
};
use crate::data_structs::context_data::ContextData;
use crate::data_structs::coords::Contig;
use crate::io::report::ReportTypeSchema;
use anyhow::{anyhow, bail};
use itertools::Itertools;
use log::warn;
use polars::prelude::*;
use polars::series::IsSorted;

/// Encodes strand information to boolean ( "+" to true, "-" to false).
pub fn encode_strand(
    lazy_frame: LazyFrame,
    strand_col: &str,
) -> LazyFrame {
    lazy_frame.with_column(
        when(col(strand_col).eq(lit("+")))
            .then(lit(true))
            .when(col(strand_col).eq(lit("-")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
            .alias("strand"),
    )
}

/// Encodes methylation context information as boolean values ("CG" to true, "CHG" to false).
pub fn encode_context(
    lazy_frame: LazyFrame,
    context_col: &str,
) -> LazyFrame {
    lazy_frame.with_column(
        when(col(context_col).eq(lit("CG")))
            .then(lit(true))
            .when(col(context_col).eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
            .alias(context_col),
    )
}

/// Decodes boolean strand information back to string representation (true to "+", false to "-", null to ".").
pub fn decode_strand(
    lazy_frame: LazyFrame,
    strand_col: &str,
    result_name: &str,
) -> LazyFrame {
    lazy_frame.with_column(
        when(col(strand_col).eq(lit(true)))
            .then(lit("+"))
            .when(col(strand_col).eq(lit(false)))
            .then(lit("-"))
            .otherwise(lit("."))
            .cast(DataType::String)
            .alias(result_name),
    )
}

/// Decodes boolean context information back to string representation (true to "CG", false to "CHG", null to "CHH").
pub fn decode_context(
    lazy_frame: LazyFrame,
    context_col: &str,
    result_name: &str,
) -> LazyFrame {
    lazy_frame.with_column(
        when(col(context_col).eq(lit(true)))
            .then(lit("CG"))
            .when(col(context_col).eq(lit(false)))
            .then(lit("CHG"))
            .otherwise(lit("CHH"))
            .cast(DataType::String)
            .alias(result_name),
    )
}

/// Builder for constructing and validating BSX batch data.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BsxBatchBuilder {
    report_type: Option<ReportTypeSchema>,
    check_nulls: bool,
    check_sorted: bool,
    check_duplicates: bool,
    rechunk: bool,
    check_single_chr: bool,
    context_data: Option<ContextData>,
    chr_dtype: Option<DataType>,
}

impl Default for BsxBatchBuilder {
    fn default() -> Self {
        Self::all_checks()
    }
}

// Public methods
impl BsxBatchBuilder {
    /// Creates a builder with all data validation checks enabled.
    pub fn all_checks() -> Self {
        Self {
            report_type: None,
            check_duplicates: true,
            check_sorted: true,
            rechunk: true,
            check_nulls: true,
            check_single_chr: true,
            context_data: None,
            chr_dtype: None,
        }
    }

    /// Creates a builder with all validation checks disabled.
    pub fn no_checks() -> Self {
        Self {
            report_type: None,
            check_duplicates: false,
            check_sorted: false,
            rechunk: false,
            check_nulls: false,
            check_single_chr: false,
            context_data: None,
            chr_dtype: None,
        }
    }

    /// Sets the report type schema for data conversion.
    pub fn with_report_type(
        mut self,
        report_type: ReportTypeSchema,
    ) -> Self {
        self.report_type = Some(report_type);
        self
    }

    /// Sets whether to check for null values in critical columns.
    pub fn with_check_nulls(
        mut self,
        check_nulls: bool,
    ) -> Self {
        self.check_nulls = check_nulls;
        self
    }

    /// Sets whether to check and ensure positions are sorted.
    pub fn with_check_sorted(
        mut self,
        check_sorted: bool,
    ) -> Self {
        self.check_sorted = check_sorted;
        self
    }

    /// Sets whether to rechunk the data for memory efficiency.
    pub fn with_rechunk(
        mut self,
        rechunk: bool,
    ) -> Self {
        self.rechunk = rechunk;
        self
    }

    /// Sets whether to check if all data is from a single chromosome.
    pub fn with_check_single_chr(
        mut self,
        check: bool,
    ) -> Self {
        self.check_single_chr = check;
        self
    }

    /// Sets whether to check for duplicate positions.
    pub fn with_check_duplicates(
        mut self,
        check_duplicates: bool,
    ) -> Self {
        self.check_duplicates = check_duplicates;
        self
    }

    /// Sets the data type for the chromosome column.
    pub fn with_chr_dtype(
        mut self,
        chr_dtype: Option<DataType>,
    ) -> Self {
        self.chr_dtype = chr_dtype;
        self
    }

    /// Sets the context data.
    pub fn with_context_data(
        mut self,
        context_data: Option<ContextData>,
    ) -> Self {
        self.context_data = context_data;
        self
    }

    /// Builds a BSX batch from the provided DataFrame.
    pub fn build<B: BsxBatchMethods + BsxTypeTag>(
        &self,
        data: DataFrame,
    ) -> anyhow::Result<B> {
        let casted: DataFrame = self.cast::<B>(data)?;
        Ok(unsafe { B::new_unchecked(self.check_modify(casted)?) })
    }

    /// Validates and optimizes a DataFrame according to builder settings.
    pub fn check_modify(
        &self,
        data: DataFrame,
    ) -> anyhow::Result<DataFrame> {
        if self.check_single_chr && data.column(CHR_NAME)?.n_unique()? != 1 {
            bail!("Multiple chromosomes")
        }
        if self.check_nulls {
            check_has_nulls(&data)?;
        }
        if self.check_duplicates
            && (data.column(POS_NAME)?.n_unique()? != data.height())
        {
            bail!("Duplicated positions")
        };
        let sorted = self.sort(data)?;
        let res = self.rechunk(sorted);
        Ok(res)
    }

    /// Checks the data for validity based on the builder settings.
    pub fn check(
        &self,
        data: &DataFrame,
    ) -> anyhow::Result<()> {
        if self.check_single_chr && data.column(CHR_NAME)?.n_unique()? != 1 {
            bail!("Multiple chromosomes")
        }
        if self.check_nulls {
            check_has_nulls(data)?;
        }
        if self.check_duplicates
            && (data.column(POS_NAME)?.n_unique()? != data.height())
        {
            bail!("Duplicated positions")
        };
        if !check_pos_ascending(data) {
            bail!("Data not sorted")
        };
        Ok(())
    }

    /// Converts a decoded batch to an encoded batch format.
    pub fn encode_batch(
        batch: BsxBatch,
        chr_dtype: DataType,
    ) -> anyhow::Result<EncodedBsxBatch> {
        let chr = batch.chr_val()?.to_string();
        let batch_data = DataFrame::from(batch);
        let mut batch_lazy = batch_data.lazy();

        let casted =
            encoded::select_cast(batch_lazy, chr_dtype, true).collect()?;
        let encoded = EncodedBsxBatch::new_with_fields(casted, chr);
        Ok(encoded)
    }

    /// Converts an encoded batch to a decoded batch format.
    pub fn decode_batch(batch: EncodedBsxBatch) -> anyhow::Result<BsxBatch> {
        let batch_data = DataFrame::from(batch);
        let mut batch_lazy = batch_data.lazy();
        batch_lazy = decode_context(batch_lazy, CONTEXT_NAME, CONTEXT_NAME);
        batch_lazy = decode_strand(batch_lazy, STRAND_NAME, STRAND_NAME);
        let casted = decoded::select_cast(batch_lazy);

        let result = casted.collect()?;
        Ok(unsafe { BsxBatch::new_unchecked(result) })
    }

    pub fn concat<B: BsxBatchMethods>(batches: Vec<B>) -> anyhow::Result<B> {
        let contigs = batches
            .iter()
            .map(|b| b.as_contig())
            .collect::<anyhow::Result<Vec<_>>>()?;
        if !contigs
            .iter()
            .map(Contig::seqname)
            .all_equal()
        {
            bail!("Chromosomes do not match")
        }
        if !contigs
            .windows(2)
            .map(|w| w[1].start() >= w[0].end())
            .all(|a| a)
        {
            bail!("Batch positions are not sorted")
        }
        let data = concat(
            batches
                .into_iter()
                .map(|b| b.take().lazy())
                .collect_vec(),
            Default::default(),
        )?
        .collect()?;
        let res: B = BsxBatchBuilder::no_checks()
            .with_rechunk(true)
            .with_check_duplicates(true)
            .with_chr_dtype(Some(
                data.schema()
                    .get(CHR_NAME)
                    .unwrap()
                    .clone(),
            ))
            .build(data)?;
        Ok(res)
    }
}

// Private methods
impl BsxBatchBuilder {
    /// Casts the DataFrame to the correct types based on the report type.
    fn cast<B: BsxBatchMethods + BsxTypeTag>(
        &self,
        data: DataFrame,
    ) -> anyhow::Result<DataFrame> {
        use BatchType as BT;
        use ReportTypeSchema as RS;

        let res = match (
            B::type_enum(),
            self.report_type,
            self.chr_dtype.clone(),
        ) {
            // BsxBatch
            (BT::Decoded, Some(RS::Bismark), _) => decoded::from_bismark(data)?,
            (BT::Decoded, Some(RS::CgMap), _) => decoded::from_cgmap(data)?,
            (BT::Decoded, Some(RS::BedGraph), _) => {
                decoded::from_bedgraph(data)?
            },
            (BT::Decoded, Some(RS::Coverage), _) => {
                decoded::from_coverage(data)?
            },
            (BT::Decoded, None, _) => decoded::from_df(data)?,
            // EncodedBsxBatch
            (BT::Encoded, Some(RS::Bismark), Some(dtype)) => {
                encoded::from_bismark(data, dtype)?
            },
            (BT::Encoded, Some(RS::CgMap), Some(dtype)) => {
                encoded::from_cgmap(data, dtype)?
            },
            (BT::Encoded, Some(RS::BedGraph), Some(dtype)) => {
                encoded::from_bedgraph(data, dtype)?
            },
            (BT::Encoded, Some(RS::Coverage), Some(dtype)) => {
                encoded::from_coverage(data, dtype)?
            },
            (BT::Encoded, None, Some(dtype)) => encoded::from_df(data, dtype)?,
            (BT::Encoded, _, None) => {
                bail!("Chr dtype must be specified for encoded conversion")
            },
            (batch_type, report_type, _) => {
                unimplemented!("Conversion from {batch_type:?} to {report_type:?} is not supported")
            },
        };
        Ok(res)
    }

    /// Rechunks the data if enabled in the builder.
    fn rechunk(
        &self,
        data: DataFrame,
    ) -> DataFrame {
        let mut data = data;
        if self.rechunk {
            data.rechunk_mut()
        }
        data
    }

    /// Sorts data by position if needed and updates sorting flags.
    fn sort(
        &self,
        data: DataFrame,
    ) -> anyhow::Result<DataFrame> {
        let mut data = data;
        if self.check_sorted && !check_pos_ascending(&data) {
            warn!("Position is not sorted");
            data = data.sort(
                ["position"],
                SortMultipleOptions::default().with_order_descending(false),
            )?;
        } else {
            let pos_col_idx = data.get_column_index(POS_NAME).unwrap();
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

    /// Selects and casts columns to the proper types for BSX batch processing.
    ///
    /// Transforms input data to a standardized schema with properly typed columns.
    pub(crate) fn select_cast(lf: LazyFrame) -> LazyFrame {
        lf.select([
            col(CHR_NAME)
                .cast(BsxBatch::chr_type())
                .alias(CHR_NAME),
            col(POS_NAME)
                .cast(BsxBatch::pos_type())
                .alias(POS_NAME),
            col(STRAND_NAME)
                .cast(BsxBatch::strand_type())
                .alias(STRAND_NAME),
            col(CONTEXT_NAME)
                .cast(BsxBatch::context_type())
                .alias(CONTEXT_NAME),
            col(COUNT_M_NAME)
                .cast(BsxBatch::count_type())
                .alias(COUNT_M_NAME),
            col(COUNT_TOTAL_NAME)
                .cast(BsxBatch::count_type())
                .alias(COUNT_TOTAL_NAME),
            col(DENSITY_NAME)
                .cast(BsxBatch::density_type())
                .alias(DENSITY_NAME),
        ])
    }

    /// Creates an expression to convert nucleotide values to strand symbols.
    fn nuc_to_strand_expr() -> Expr {
        when(col("nuc").eq(lit("C")))
            .then(lit("+"))
            .when(col("nuc").eq(lit("G")))
            .then(lit("-"))
            .otherwise(lit("."))
    }

    /// Creates a DataFrame from a decoded source
    pub(crate) fn from_df(df: DataFrame) -> anyhow::Result<DataFrame> {
        let casted = select_cast(df.lazy());
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Processes Bismark format data into standardized BSX format.
    pub(crate) fn from_bismark(df: DataFrame) -> anyhow::Result<DataFrame> {
        let col_added = df
            .lazy()
            .with_column(count_total_col_expr())
            .with_column(density_col_expr());
        let casted = select_cast(col_added);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Converts CG map format data into standardized BSX format.
    pub(crate) fn from_cgmap(df: DataFrame) -> anyhow::Result<DataFrame> {
        let mut col_added = df
            .lazy()
            .with_column(nuc_to_strand_expr().alias(STRAND_NAME));
        let casted = select_cast(col_added);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Processes coverage format data into standardized BSX format.
    pub(crate) fn from_coverage(df: DataFrame) -> anyhow::Result<DataFrame> {
        let col_added = df
            .lazy()
            .with_column(count_total_col_expr())
            .with_columns([
                col("start").alias("position"),
                lit(".").alias("strand"),
                lit(NULL).alias("context"),
            ]);
        let casted = select_cast(col_added);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Processes BedGraph format data into standardized BSX format.
    pub(crate) fn from_bedgraph(df: DataFrame) -> anyhow::Result<DataFrame> {
        let col_added = df.lazy().with_columns([
            col("start").alias("position"),
            lit(".").alias("strand"),
            lit(NULL).alias("context"),
            lit(NULL).alias("count_m"),
            lit(NULL).alias("count_total"),
        ]);
        let casted = select_cast(col_added);
        casted.collect().map_err(|e| anyhow!(e))
    }
}

mod encoded {
    use super::*;
    use crate::data_structs::batch::encoded::EncodedBsxBatch;

    /// Encodes context information as boolean values ("CG" to true, "CHG" to false).
    fn encode_context() -> Expr {
        when(col(CONTEXT_NAME).eq(lit("CG")))
            .then(lit(true))
            .when(col(CONTEXT_NAME).eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
            .alias(CONTEXT_NAME)
    }

    /// Encodes strand information to boolean ( "+" to true, "-" to false).
    fn encode_strand() -> Expr {
        when(col(STRAND_NAME).eq(lit("+")))
            .then(lit(true))
            .when(col(STRAND_NAME).eq(lit("-")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
            .alias(STRAND_NAME)
    }

    /// Creates a DataFrame from an encoded source
    pub(crate) fn from_df(
        df: DataFrame,
        chr_dtype: DataType,
    ) -> anyhow::Result<DataFrame> {
        let casted = select_cast(df.lazy(), chr_dtype, false);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Selects and casts columns to the proper types for BSX batch processing.
    ///
    /// Transforms input data to a standardized schema with properly typed columns.
    pub(crate) fn select_cast(
        lf: LazyFrame,
        chr_dtype: DataType,
        encode: bool,
    ) -> LazyFrame {
        lf.select([
            col("chr")
                .cast(chr_dtype)
                .alias(CHR_NAME),
            col("position")
                .cast(EncodedBsxBatch::pos_type())
                .alias(POS_NAME),
            if encode {
                encode_strand()
            } else {
                col("strand")
                    .cast(EncodedBsxBatch::strand_type())
                    .alias(STRAND_NAME)
            },
            if encode {
                encode_context()
            } else {
                col("context")
                    .cast(EncodedBsxBatch::context_type())
                    .alias(CONTEXT_NAME)
            },
            col("count_m")
                .cast(EncodedBsxBatch::count_type())
                .alias(COUNT_M_NAME),
            col("count_total")
                .cast(EncodedBsxBatch::count_type())
                .alias(COUNT_TOTAL_NAME),
            col("density")
                .cast(EncodedBsxBatch::density_type())
                .alias(DENSITY_NAME),
        ])
    }

    /// Creates an expression to convert nucleotide values to strand symbols.
    fn nuc_to_strand_expr() -> Expr {
        when(col("nuc").eq(lit("C")))
            .then(lit("+"))
            .when(col("nuc").eq(lit("G")))
            .then(lit("-"))
            .otherwise(lit(NULL))
            .alias(STRAND_NAME)
    }

    /// Processes Bismark format data into standardized BSX format.
    pub(crate) fn from_bismark(
        df: DataFrame,
        chr_dtype: DataType,
    ) -> anyhow::Result<DataFrame> {
        let col_added = df
            .lazy()
            .with_column(count_total_col_expr())
            .with_column(density_col_expr());
        let casted = select_cast(col_added, chr_dtype, true);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Converts CG map format data into standardized BSX format.
    pub(crate) fn from_cgmap(
        df: DataFrame,
        chr_dtype: DataType,
    ) -> anyhow::Result<DataFrame> {
        let col_added = df
            .lazy()
            .with_column(nuc_to_strand_expr());
        let casted = select_cast(col_added, chr_dtype, true);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Processes coverage format data into standardized BSX format.
    pub(crate) fn from_coverage(
        df: DataFrame,
        chr_dtype: DataType,
    ) -> anyhow::Result<DataFrame> {
        let col_added = df
            .lazy()
            .with_column(count_total_col_expr())
            .with_columns([
                col("start").alias("position"),
                lit(NULL).alias("strand"),
                lit(NULL).alias("context"),
            ]);
        let casted = select_cast(col_added, chr_dtype, true);
        casted.collect().map_err(|e| anyhow!(e))
    }

    /// Processes BedGraph format data into standardized BSX format.
    pub(crate) fn from_bedgraph(
        df: DataFrame,
        chr_dtype: DataType,
    ) -> anyhow::Result<DataFrame> {
        let col_added = df.lazy().with_columns([
            col("start").alias("position"),
            lit(NULL).alias("strand"),
            lit(NULL).alias("context"),
            lit(NULL).alias("count_m"),
            lit(NULL).alias("count_total"),
        ]);
        let casted = select_cast(col_added, chr_dtype, true);
        casted.collect().map_err(|e| anyhow!(e))
    }
}

/// Validates that critical columns do not contain null values.
fn check_has_nulls(df: &DataFrame) -> anyhow::Result<()> {
    if df
        .column(CHR_NAME)
        .expect("chr column not found")
        .null_count()
        > 0
    {
        bail!("Nulls not allowed in 'chr' column");
    }
    if df
        .column(POS_NAME)
        .expect("position column not found")
        .null_count()
        > 0
    {
        bail!("Nulls not allowed in 'position' column");
    }
    Ok(())
}

/// Checks if position column values are in ascending order.
fn check_pos_ascending(df: &DataFrame) -> bool {
    let pos = df
        .column(POS_NAME)
        .expect("position column not found");
    pos.as_series()
        .unwrap()
        .iter()
        .is_sorted()
}

/// Creates an expression to compute the total count from methylated and unmethylated counts.
fn count_total_col_expr() -> Expr {
    (col("count_m") + col("count_um")).alias("count_total")
}

/// Creates an expression to compute methylation density from count data.
fn density_col_expr() -> Expr {
    (col("count_m").cast(DataType::Float64)
        / col("count_total").cast(DataType::Float64))
    .alias("density")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_structs::batch::{
        decoded::BsxBatch, traits::BsxBatchMethods,
    };
    use polars::df;
    use polars::prelude::NamedFrom;
    use rstest::rstest;

    // Helper function to create a basic DataFrame for testing
    fn create_test_df() -> DataFrame {
        df!(
            "chr" => ["chr1", "chr1", "chr1"],
            "position" => [100u32, 200, 150], // Intentionally unsorted
            "strand" => ["+", "-", "+"],
            "context" => ["CG", "CHG", "CHH"],
            "count_m" => [10u32, 5, 8],
            "count_total" => [20u32, 15, 10],
            "density" => [0.5f32, 0.333, 0.8]
        )
        .unwrap()
    }

    // Helper function to create a basic DataFrame with nulls
    fn create_test_df_with_nulls() -> DataFrame {
        df!(
            "chr" => [Some("chr1"), None, Some("chr1")], // Null in chr
            "position" => [Some(100u32), Some(200), None], // Null in position
            "strand" => ["+", "-", "+"],
            "context" => ["CG", "CHG", "CHH"],
            "count_m" => [10u32, 5, 8],
            "count_total" => [20u32, 15, 10],
            "density" => [0.5f32, 0.333, 0.8]
        )
        .unwrap()
    }

    // Helper function to create a basic DataFrame with duplicates
    fn create_test_df_with_duplicates() -> DataFrame {
        df!(
            "chr" => ["chr1", "chr1", "chr1"],
            "position" => [100u32, 200, 100], // Duplicate position
            "strand" => ["+", "-", "+"],
            "context" => ["CG", "CHG", "CHH"],
            "count_m" => [10u32, 5, 8],
            "count_total" => [20u32, 15, 10],
            "density" => [0.5f32, 0.333, 0.8]
        )
        .unwrap()
    }

    #[test]
    fn test_check_modify_sorts_when_needed() -> anyhow::Result<()> {
        let builder = BsxBatchBuilder::all_checks();
        let df = create_test_df(); // Unsorted
        let modified_df = builder.check_modify(df)?;

        let pos_series = modified_df.column(POS_NAME)?;
        assert!(check_pos_ascending(&modified_df));
        assert_eq!(
            pos_series
                .u32()?
                .into_no_null_iter()
                .collect::<Vec<_>>(),
            vec![100, 150, 200]
        );
        Ok(())
    }

    #[test]
    fn test_check_modify_detects_nulls() {
        let builder = BsxBatchBuilder::all_checks();
        let df = create_test_df_with_nulls();
        let result = builder.check_modify(df);
        assert!(result.is_err());
    }

    #[test]
    fn test_check_modify_ignores_nulls_when_disabled() -> anyhow::Result<()> {
        let builder = BsxBatchBuilder::no_checks(); // Disable null check
        let df = create_test_df_with_nulls();
        // Should not error, even with nulls, because check is off
        let result = builder.check_modify(df);
        assert!(result.is_ok());
        Ok(())
    }

    #[test]
    fn test_check_modify_detects_duplicates() {
        let builder = BsxBatchBuilder::all_checks();
        let df = create_test_df_with_duplicates();
        let result = builder.check_modify(df);
        assert!(result.is_err());
    }

    #[test]
    fn test_check_modify_ignores_duplicates_when_disabled() -> anyhow::Result<()>
    {
        let builder = BsxBatchBuilder::no_checks(); // Disable duplicate check
        let df = create_test_df_with_duplicates();
        // Should not error, even with duplicates, because check is off
        let result = builder.check_modify(df);
        assert!(result.is_ok());
        Ok(())
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

    // --- Test Decoded Conversions ---

    #[rstest]
    #[case(ReportTypeSchema::Bismark, create_bismark_df())]
    #[case(ReportTypeSchema::CgMap, create_cgmap_df())]
    #[case(ReportTypeSchema::Coverage, create_coverage_df())]
    #[case(ReportTypeSchema::BedGraph, create_bedgraph_df())]
    fn test_build_decoded_from_report_type(
        #[case] report_type: ReportTypeSchema,
        #[case] input_df: DataFrame,
    ) -> anyhow::Result<()> {
        let builder =
            BsxBatchBuilder::no_checks().with_report_type(report_type);
        let batch: BsxBatch = builder.build(input_df)?;
        let df = DataFrame::from(batch);

        // Check common columns and types
        assert_eq!(df.column(CHR_NAME)?.dtype(), &BsxBatch::chr_type());
        assert_eq!(df.column(POS_NAME)?.dtype(), &BsxBatch::pos_type());
        assert_eq!(df.column(STRAND_NAME)?.dtype(), &BsxBatch::strand_type());
        assert_eq!(df.column(CONTEXT_NAME)?.dtype(), &BsxBatch::context_type());
        assert_eq!(df.column(COUNT_M_NAME)?.dtype(), &BsxBatch::count_type());
        assert_eq!(
            df.column(COUNT_TOTAL_NAME)?.dtype(),
            &BsxBatch::count_type()
        );
        assert_eq!(df.column(DENSITY_NAME)?.dtype(), &BsxBatch::density_type());

        // Check specific transformations (example for CgMap strand)
        if report_type == ReportTypeSchema::CgMap {
            let expected_strand = Series::new(STRAND_NAME.into(), ["+", "-"]);
            assert!(df
                .column(STRAND_NAME)?
                .equals(&expected_strand.into_column()));
        }
        // Check specific transformations (example for Bismark/Coverage count_total)
        if report_type == ReportTypeSchema::Bismark
            || report_type == ReportTypeSchema::Coverage
        {
            let count_m = df.column(COUNT_M_NAME)?.u32()?;
            let count_total = df.column(COUNT_TOTAL_NAME)?.u32()?;
            assert!(count_m
                .into_iter()
                .zip(count_total.into_iter())
                .all(|(m, t)| m <= t));
        }
        // Check specific transformations (example for BedGraph nulls)
        if report_type == ReportTypeSchema::BedGraph {
            assert!(df.column(COUNT_M_NAME)?.null_count() > 0);
            assert!(
                df.column(COUNT_TOTAL_NAME)?
                    .null_count()
                    > 0
            );
            assert!(df.column(CONTEXT_NAME)?.null_count() > 0);
        }

        Ok(())
    }

    // --- Test Encoded Conversions ---
    #[rstest]
    #[case(ReportTypeSchema::Bismark, create_bismark_df())]
    #[case(ReportTypeSchema::CgMap, create_cgmap_df())]
    #[case(ReportTypeSchema::Coverage, create_coverage_df())]
    #[case(ReportTypeSchema::BedGraph, create_bedgraph_df())]
    fn test_build_encoded_from_report_type(
        #[case] report_type: ReportTypeSchema,
        #[case] input_df: DataFrame,
    ) -> anyhow::Result<()> {
        let chr_dtype = DataType::Categorical(None, Default::default()); // Use Categorical for encoded chr
        let builder = BsxBatchBuilder::no_checks()
            .with_report_type(report_type)
            .with_chr_dtype(Some(chr_dtype.clone())); // Must provide chr_dtype for encoded

        let batch: EncodedBsxBatch = builder.build(input_df)?;
        let df = DataFrame::from(batch);

        // Check common columns and types
        assert_eq!(df.column(CHR_NAME)?.dtype(), &chr_dtype); // Check specified dtype
        assert_eq!(df.column(POS_NAME)?.dtype(), &EncodedBsxBatch::pos_type());
        assert_eq!(
            df.column(STRAND_NAME)?.dtype(),
            &EncodedBsxBatch::strand_type()
        );
        assert_eq!(
            df.column(CONTEXT_NAME)?.dtype(),
            &EncodedBsxBatch::context_type()
        );
        assert_eq!(
            df.column(COUNT_M_NAME)?.dtype(),
            &EncodedBsxBatch::count_type()
        );
        assert_eq!(
            df.column(COUNT_TOTAL_NAME)?.dtype(),
            &EncodedBsxBatch::count_type()
        );
        assert_eq!(
            df.column(DENSITY_NAME)?.dtype(),
            &EncodedBsxBatch::density_type()
        );

        // Check specific transformations (example for CgMap strand - encoded as bool)
        if report_type == ReportTypeSchema::CgMap {
            let expected_strand =
                Series::new(STRAND_NAME.into(), [Some(true), Some(false)]); // C -> true, G -> false
            assert!(df
                .column(STRAND_NAME)?
                .equals(&expected_strand.into_column()));
        }
        // Check specific transformations (example for BedGraph nulls - strand/context should be null)
        if report_type == ReportTypeSchema::BedGraph {
            assert_eq!(df.column(STRAND_NAME)?.null_count(), df.height());
            assert_eq!(df.column(CONTEXT_NAME)?.null_count(), df.height());
        }

        Ok(())
    }

    #[test]
    fn test_build_encoded_fails_without_chr_dtype() {
        let builder = BsxBatchBuilder::no_checks()
            .with_report_type(ReportTypeSchema::Bismark); // No chr_dtype
        let df = create_bismark_df();
        let result: anyhow::Result<EncodedBsxBatch> = builder.build(df);
        assert!(result.is_err());
    }

    // --- Test Encode/Decode Batch ---
    fn create_decoded_test_df() -> DataFrame {
        df!(
            CHR_NAME => ["chrTest", "chrTest", "chrTest"],
            POS_NAME => [100u32, 150, 200], // Sorted
            STRAND_NAME => ["+", "-", "+"],
            CONTEXT_NAME => ["CG", "CHG", "CHH"], // Mix of contexts
            COUNT_M_NAME => [10u32, 5, 8],
            COUNT_TOTAL_NAME => [20u32, 15, 10],
            DENSITY_NAME => [0.5f32, 0.333, 0.8]
        )
        .unwrap()
    }

    #[test]
    fn test_encode_decode_batch_cycle() -> anyhow::Result<()> {
        // 1. Create a valid decoded batch
        let decoded_df = create_decoded_test_df();
        let original_batch = unsafe { BsxBatch::new_unchecked(decoded_df) };

        // 2. Encode the batch
        let chr_dtype = DataType::Categorical(None, Default::default());
        let encoded_batch = BsxBatchBuilder::encode_batch(
            original_batch.clone(),
            chr_dtype.clone(),
        )?;
        let encoded_df = DataFrame::from(encoded_batch.clone());

        // 3. Check encoded batch properties
        assert_eq!(encoded_batch.chr_val()?, "chrTest");
        assert_eq!(encoded_df.column(CHR_NAME)?.dtype(), &chr_dtype);
        assert_eq!(
            encoded_df.column(STRAND_NAME)?.dtype(),
            &EncodedBsxBatch::strand_type()
        ); // Boolean
        assert_eq!(
            encoded_df.column(CONTEXT_NAME)?.dtype(),
            &EncodedBsxBatch::context_type()
        ); // Boolean

        // Check encoded values: strand ["+", "-", "+"] -> [true, false, true]
        let expected_strand = Series::new(
            STRAND_NAME.into(),
            [Some(true), Some(false), Some(true)],
        );
        assert_eq!(
            encoded_df
                .column(STRAND_NAME)?
                .as_materialized_series(),
            &expected_strand
        );

        // Check encoded values: context ["CG", "CHG", "CHH"] -> [true, false, null] (CHH -> null)
        let expected_context = Series::new(
            CONTEXT_NAME.into(),
            [Some(true), Some(false), None::<bool>],
        );
        assert_eq!(
            encoded_df
                .column(CONTEXT_NAME)?
                .as_materialized_series(),
            &expected_context
        );

        // 4. Decode the batch
        let decoded_batch = BsxBatchBuilder::decode_batch(encoded_batch)?;

        // 5. Verify the decoded batch matches the original
        // Compare the DataFrames directly
        let original_df = DataFrame::from(original_batch);
        let final_df = DataFrame::from(decoded_batch);

        let cols_to_compare = [
            CHR_NAME,
            POS_NAME,
            STRAND_NAME,
            CONTEXT_NAME,
            COUNT_M_NAME,
            COUNT_TOTAL_NAME,
        ];
        assert_eq!(
            original_df.select(cols_to_compare)?,
            final_df.select(cols_to_compare)?
        );

        Ok(())
    }
}
