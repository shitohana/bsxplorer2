use anyhow::{Context as _, Result, anyhow};
use polars::{prelude::*, series::IsSorted};
use std::collections::HashSet;

// This import is marked as unused in the original code and isn't used here.
// use super::bsx_dtype::BsxPolarsDtype;

/// Expected column names in the DataFrame.
const COL_NAMES: [&str; 7] = [
    "chr",
    "position",
    "strand",
    "context",
    "count_m",
    "count_total",
    "density",
];

/// Trait representing a batch of BSX data backed by a Polars DataFrame.
/// Provides validated access to specific columns and common operations.
trait BsxDataBatch {
    // Associated types defining the expected Polars DataType for each column.
    type Chr: PolarsDataType;
    type Position: PolarsDataType;
    type Strand: PolarsDataType;
    type Context: PolarsDataType;
    type Count: PolarsDataType;
    type Density: PolarsDataType;

    /// Creates a new instance from a DataFrame without performing checks.
    ///
    /// # Safety
    /// The caller must ensure that the DataFrame adheres to the required schema
    /// (column names, types, non-null constraints) before calling this function.
    unsafe fn new_unchecked(data: DataFrame) -> Self;

    /// Returns an immutable reference to the underlying DataFrame.
    fn data(&self) -> &DataFrame;

    /// Returns a mutable reference to the underlying DataFrame.
    fn data_mut(&mut self) -> &mut DataFrame;

    /// Consumes the `BsxDataBatch` and returns a LazyFrame.
    fn to_lazy(self) -> LazyFrame where Self: Sized;

    /// Attempts to create a new `BsxDataBatch` from a DataFrame, performing validations.
    ///
    /// Checks for:
    /// - Non-empty DataFrame.
    /// - Presence of all required columns (`COL_NAMES`).
    /// - Correct data types for required columns (deferred to accessor methods).
    /// - Non-null values in 'chr' and 'position' columns.
    ///
    /// Aligns chunks in the DataFrame for potentially better performance in subsequent operations.
    fn try_new(df: DataFrame) -> Result<Self> where Self: Sized {
        if df.height() == 0 {
            return Err(anyhow!("Cannot create BsxDataBatch from an empty DataFrame"));
        }

        // Validate column presence first.
        if !Self::validate_columns(&df) {
            return Err(anyhow!("DataFrame is missing required columns"));
        }

        // Check essential columns for nulls.
        if Self::check_col_has_null(&df, "chr")? {
            return Err(anyhow!("Column 'chr' contains unexpected null values"));
        }
        if Self::check_col_has_null(&df, "position")? {
            return Err(anyhow!("Column 'position' contains unexpected null values"));
        }

        // Safety: We've performed the necessary checks (column presence, non-null core columns).
        // Type checks are implicitly handled by the accessor methods later.
        let mut new_self = unsafe {
            Self::new_unchecked(df)
        };

        // Align chunks for potentially better performance downstream.
        // This can be computationally intensive for large DataFrames.
        new_self.align_chunks();

        Ok(new_self)
    }

    /// Gets the 'chr' column as a Series, checking its data type.
    fn chr(&self) -> Result<&Series> {
        let series = self.data()
            .column("chr")
            .context("Failed to access 'chr' column")?;
        // `as_materialized_series` just provides a reference, no computation needed here usually.
        let series = series.as_materialized_series();
        let expected_dtype = Self::Chr::get_dtype();
        if series.dtype() != &expected_dtype {
            return Err(anyhow!(
                "Column 'chr' has incorrect dtype (Expected: {}, Found: {})",
                expected_dtype, series.dtype()
            ));
        }
        Ok(series)
    }

    /// Gets the 'position' column as a Series, checking its data type.
    fn position(&self) -> Result<&Series> {
        let series = self.data()
            .column("position")
            .context("Failed to access 'position' column")?;
        let series = series.as_materialized_series();
        let expected_dtype = Self::Position::get_dtype();
        if series.dtype() != &expected_dtype {
            return Err(anyhow!(
                "Column 'position' has incorrect dtype (Expected: {}, Found: {})",
                expected_dtype, series.dtype()
            ));
        }
        Ok(series)
    }

    /// Gets the 'strand' column as a Series, checking its data type.
    fn strand(&self) -> Result<&Series> {
        let series = self.data()
            .column("strand")
            .context("Failed to access 'strand' column")?;
        let series = series.as_materialized_series();
        let expected_dtype = Self::Strand::get_dtype();
        if series.dtype() != &expected_dtype {
            return Err(anyhow!(
                "Column 'strand' has incorrect dtype (Expected: {}, Found: {})",
                expected_dtype, series.dtype()
            ));
        }
        Ok(series)
    }

    /// Gets the 'context' column as a Series, checking its data type.
    fn context(&self) -> Result<&Series> {
        let series = self.data()
            .column("context")
            .context("Failed to access 'context' column")?;
        let series = series.as_materialized_series();
        let expected_dtype = Self::Context::get_dtype();
        if series.dtype() != &expected_dtype {
            return Err(anyhow!(
                "Column 'context' has incorrect dtype (Expected: {}, Found: {})",
                expected_dtype, series.dtype()
            ));
        }
        Ok(series)
    }

    /// Gets the 'count_m' column as a Series, checking its data type.
    fn count_m(&self) -> Result<&Series> {
        let series = self.data()
            .column("count_m")
            .context("Failed to access 'count_m' column")?;
        let series = series.as_materialized_series();
        let expected_dtype = Self::Count::get_dtype();
        if series.dtype() != &expected_dtype {
            return Err(anyhow!(
                "Column 'count_m' has incorrect dtype (Expected: {}, Found: {})",
                expected_dtype, series.dtype()
            ));
        }
        Ok(series)
    }

    /// Gets the 'count_total' column as a Series, checking its data type.
    fn count_total(&self) -> Result<&Series> {
        let series = self.data()
            .column("count_total")
            .context("Failed to access 'count_total' column")?;
        let series = series.as_materialized_series();
        let expected_dtype = Self::Count::get_dtype();
        if series.dtype() != &expected_dtype {
            return Err(anyhow!(
                "Column 'count_total' has incorrect dtype (Expected: {}, Found: {})",
                expected_dtype, series.dtype()
            ));
        }
        Ok(series)
    }

    /// Gets the 'density' column as a Series, checking its data type.
    fn density(&self) -> Result<&Series> {
        let series = self.data()
            .column("density")
            .context("Failed to access 'density' column")?;
        let series = series.as_materialized_series();
        let expected_dtype = Self::Density::get_dtype();
        if series.dtype() != &expected_dtype {
            return Err(anyhow!(
                "Column 'density' has incorrect dtype (Expected: {}, Found: {})",
                expected_dtype, series.dtype()
            ));
        }
        Ok(series)
    }

    /// Validates that the DataFrame contains all required columns.
    fn validate_columns(df: &DataFrame) -> bool
    where
        Self: Sized, {
        df.schema()
            .iter_names()
            .all(|name| COL_NAMES.contains(&name.as_str()))
    }

    /// Checks if a specific column in the DataFrame contains any null values.
    /// Returns an error if the column doesn't exist.
    fn check_col_has_null(df: &DataFrame, col: &str) -> Result<bool> {
        df.column(col)
            .with_context(|| format!("Failed to access column '{}' for null check", col))
            .map(|s| s.null_count() > 0)
    }

    /// Checks if the 'chr' column contains only a single unique value.
    /// Returns an error if the column doesn't exist or unique check fails.
    /// Note: This method might be unused depending on the specific implementation using the trait.
    #[allow(dead_code)] // Mark as potentially unused
    fn check_chr_unique(df: &DataFrame) -> Result<bool> {
        let chr_series = df.column("chr")
            .context("Failed to access 'chr' column for unique check")?;
        let n_unique = chr_series.n_unique()
            .context("Failed to count unique values in 'chr' column")?;
        Ok(n_unique == 1)
    }

    fn check_sorted(&mut self) -> Result<bool> {
        let was_sorted = self.data().column("position").unwrap().is_sorted_flag();

        match was_sorted {
            IsSorted::Ascending => Ok(true),
            _ => {
                let pos_series = self.position()
                    .context("Failed to get 'position' series for sorted check")?;
                let pos_arr = pos_series.to_physical_repr();

                let is_sorted_new = match pos_arr.dtype() {
                    DataType::UInt32 => {
                        pos_arr.u32().expect("Cast failed").iter().is_sorted()
                    },
                    DataType::UInt64 => {
                        pos_arr.u64().expect("Cast failed").iter().is_sorted()
                    },
                    _ => unreachable!("Unexpected DataType of position column")
                };

                if is_sorted_new {
                    let pos_index = self.data().get_column_index("position").unwrap();
                    let cols = unsafe { self.data_mut().get_columns_mut() };
                    let pos_col = &mut cols[pos_index];
                    pos_col.set_sorted_flag(IsSorted::Ascending);
                }
                Ok(is_sorted_new)
            }
        }
    }


    /// Rechunk the DataFrame to have contiguous memory layout for each column.
    /// Uses parallel processing if the "parallel" feature is enabled in Polars.
    fn align_chunks(&mut self) {
        // Consider making this optional or configurable, as it can be costly.
        self.data_mut().align_chunks_par();
    }

    /// Filters the `BsxDataBatch` based on a Polars expression predicate.
    ///
    /// Returns a new `BsxDataBatch` containing only the rows that satisfy the predicate.
    ///
    /// # Safety
    /// This method uses `new_unchecked` internally, assuming that the filter operation
    /// preserves the necessary schema invariants (column existence, types, non-null constraints).
    /// Polars filters generally maintain these properties.
    fn filter(self, predicate: Expr) -> Result<Self> where Self: Sized {
        let ldf = self.to_lazy().filter(predicate);
        let df = ldf.collect()
            .context("Failed to collect DataFrame after applying filter")?;

        // Safety: Assuming the filter operation doesn't invalidate the DataFrame structure
        // (e.g., remove required columns, change types, introduce nulls in chr/position).
        Ok(unsafe { Self::new_unchecked(df) })
    }

    /// Returns the number of rows in the underlying DataFrame.
    fn len(&self) -> usize {
        self.data().height()
    }

    /// Checks if the DataFrame is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
