use itertools::Itertools;
use polars::frame::DataFrame;
use polars::prelude::*;

use super::{
    BsxBatch,
    BsxColumns,
};
use crate::plsmallstr;

/// Merges multiple `BsxBatch` replicates into a single batch.
///
/// Assumes all input batches have the same length and identical 'chr', 'pos',
/// 'strand', and 'context' columns.
///
/// # Arguments
///
/// * `batches` - A vector of `BsxBatch` instances to merge.
/// * `count_agg` - A function that takes a vector of `Column` references (for
///   CountM or CountTotal) and returns a single aggregated `Column`.
/// * `density_agg` - A function that takes a vector of `Column` references (for
///   Density) and returns a single aggregated `Column`.
///
/// # Returns
///
/// A `Result` containing the merged `BsxBatch` or an `anyhow::Error`.
#[allow(clippy::type_complexity)]
pub fn merge_replicates(
    mut batches: Vec<BsxBatch>,
    count_agg: Box<dyn Fn(Vec<&Column>) -> Column>,
    density_agg: Box<dyn Fn(Vec<&Column>) -> Column>,
) -> PolarsResult<BsxBatch> {
    if batches.is_empty() {
        polars_bail!(InvalidOperation: "batches cannot be empty");
    }
    else if batches.len() == 1 {
        return Ok(batches.pop().unwrap());
    }
    else {
        // Assume all batches have the same length and identical chr, pos,
        // strand, context This should be guaranteed by the caller or
        // previous steps (e.g., alignment)
        let first_batch = unsafe { batches.get_unchecked(0) };
        let chr_col = first_batch.data().column(BsxColumns::Chr.as_str())?;
        let pos_col = first_batch.data().column(BsxColumns::Position.as_str())?;
        let strand_col = first_batch.data().column(BsxColumns::Strand.as_str())?;
        let context_col = first_batch.data().column(BsxColumns::Context.as_str())?;

        let count_m_col = count_agg(
            batches
                .iter()
                .map(|b| b.data().column(BsxColumns::CountM.as_str()).unwrap())
                .collect_vec(),
        )
        .with_name(plsmallstr!(BsxColumns::CountM.as_str())); // Ensure correct name

        let count_total_col = count_agg(
            batches
                .iter()
                .map(|b| b.data().column(BsxColumns::CountTotal.as_str()).unwrap())
                .collect_vec(),
        )
        .with_name(plsmallstr!(BsxColumns::CountTotal.as_str())); // Ensure correct name

        let density_col = density_agg(
            batches
                .iter()
                .map(|b| b.data().column(BsxColumns::Density.as_str()).unwrap())
                .collect_vec(),
        )
        .with_name(plsmallstr!(BsxColumns::Density.as_str())); // Ensure correct name

        // Ensure columns have the correct data type for the target batch type B
        let count_m_col = count_m_col.cast(&BsxColumns::CountM.dtype())?;
        let count_total_col = count_total_col.cast(&BsxColumns::CountTotal.dtype())?;
        let density_col = density_col.cast(&BsxColumns::Density.dtype())?;

        let df = DataFrame::from_iter([
            chr_col.to_owned(),
            pos_col.to_owned(),
            strand_col.to_owned(),
            context_col.to_owned(),
            count_m_col,
            count_total_col,
            density_col,
        ]);

        // Create the batch without checks, assuming input batches were valid
        // and aligned

        Ok(unsafe { BsxBatch::new_unchecked(df) })
    }
}

/// Creates a Polars Categorical `DataType` from a list of chromosome values.
///
/// # Arguments
///
/// * `chr_values` - A reference to a slice of optional string-like values
///   representing chromosomes.
pub fn create_caregorical_dtype<S, P>(chr_values: P) -> DataType
where
    S: AsRef<str>,
    P: AsRef<[Option<S>]>, {
    use polars::export::arrow::array::Utf8ViewArray;
    let categories = Utf8ViewArray::from_slice(chr_values);
    let rev_mapping = Arc::new(RevMapping::build_local(categories));
    DataType::Categorical(Some(rev_mapping), CategoricalOrdering::Physical)
}

/// Creates an empty Polars Categorical `DataType`.
pub const fn create_empty_categorical_dtype() -> DataType {
    DataType::Categorical(None, CategoricalOrdering::Physical)
}

#[macro_export]
macro_rules! name_dtype_tuple {
    ($enum_var: expr) => {
        ($enum_var.as_str().into(), $enum_var.dtype())
    };
}
#[macro_export]
macro_rules! create_empty_series {
    ($col: ident) => {
        Series::new_empty(plsmallstr!(BsxCol::$col.as_str()), &BsxCol::$col.dtype())
    };
}
#[macro_export]
macro_rules! get_col_fn {
    ($name: ident, $col: expr, $col_fn: ident, $rettype: ty) => {
        pub fn $name(&self) -> &$rettype {
            self.data().column($col).unwrap().$col_fn().unwrap()
        }
    };
}

pub(crate) use {
    create_empty_series,
    get_col_fn,
    name_dtype_tuple,
};
