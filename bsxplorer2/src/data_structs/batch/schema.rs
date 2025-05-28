use std::any::Any;

use polars::prelude::*;

use super::{
    create_empty_categorical_dtype,
    name_dtype_tuple,
};
use crate::plsmallstr;

/// Represents the columns expected in BSX data.
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum BsxColumns {
    Chr,
    Position,
    Strand,
    Context,
    CountM,
    CountTotal,
    Density,
}


impl BsxColumns {
    /// Returns the Polars Schema for the BSX columns.
    pub fn schema() -> Schema {
        Schema::from_iter([
            name_dtype_tuple!(BsxColumns::Chr),
            name_dtype_tuple!(BsxColumns::Position),
            name_dtype_tuple!(BsxColumns::Strand),
            name_dtype_tuple!(BsxColumns::Context),
            name_dtype_tuple!(BsxColumns::CountM),
            name_dtype_tuple!(BsxColumns::CountTotal),
            name_dtype_tuple!(BsxColumns::Density),
        ])
    }

    /// Returns the string representation of the column name.
    pub const fn as_str(&self) -> &'static str {
        match self {
            BsxColumns::Chr => "chr",
            BsxColumns::Position => "position",
            BsxColumns::Strand => "strand",
            BsxColumns::Context => "context",
            BsxColumns::CountM => "count_m",
            BsxColumns::CountTotal => "count_total",
            BsxColumns::Density => "density",
        }
    }

    /// Returns the Polars DataType for the column.
    pub const fn dtype(&self) -> DataType {
        match self {
            BsxColumns::Chr => create_empty_categorical_dtype(),
            BsxColumns::Position => DataType::UInt32,
            BsxColumns::Strand => DataType::Boolean,
            BsxColumns::Context => DataType::Boolean,
            BsxColumns::CountM => DataType::UInt16,
            BsxColumns::CountTotal => DataType::UInt16,
            BsxColumns::Density => DataType::Float32,
        }
    }

    /// Returns an array containing all BSX column names as strings.
    pub const fn colnames() -> [&'static str; 7] {
        [
            BsxColumns::Chr.as_str(),
            BsxColumns::Position.as_str(),
            BsxColumns::Strand.as_str(),
            BsxColumns::Context.as_str(),
            BsxColumns::CountM.as_str(),
            BsxColumns::CountTotal.as_str(),
            BsxColumns::Density.as_str(),
        ]
    }

    /// Checks if the given string matches any of the BSX column names.
    pub fn has_name(name: &str) -> bool {
        Self::colnames().contains(&name)
    }

    /// Creates a Polars expression (Expr) referencing this column.
    #[inline(always)]
    pub fn col(&self) -> Expr {
        col(self.as_str())
    }

    /// Attempts to create a Polars AnyValue from a boxed value based on the
    /// column's expected type.
    ///
    /// # Arguments
    ///
    /// * `value`: A boxed value that is expected to match the column's
    ///   DataType.
    ///
    /// # Returns
    ///
    /// Returns `Some(AnyValue)` if the value can be downcast to the expected
    /// type, otherwise returns `None`.
    pub fn create_anyvalue(
        &self,
        value: Box<dyn Any>,
    ) -> Option<AnyValue> {
        match self {
            BsxColumns::Chr => {
                value
                    .downcast_ref::<String>()
                    .map(|v| AnyValue::StringOwned(v.into()))
            },
            BsxColumns::Position => {
                value.downcast_ref::<u32>().map(|v| AnyValue::UInt32(*v))
            },
            BsxColumns::Strand => {
                value.downcast_ref::<bool>().map(|v| AnyValue::Boolean(*v))
            },
            BsxColumns::Context => {
                value
                    .downcast_ref::<Option<bool>>()
                    .map(|v| v.map(AnyValue::Boolean).unwrap_or(AnyValue::Null))
            },
            BsxColumns::CountM => {
                value.downcast_ref::<u16>().map(|v| AnyValue::UInt16(*v))
            },
            BsxColumns::CountTotal => {
                value.downcast_ref::<u16>().map(|v| AnyValue::UInt16(*v))
            },
            BsxColumns::Density => {
                value.downcast_ref::<f32>().map(|v| AnyValue::Float32(*v))
            },
        }
    }

    /// Creates a Polars Series from a vector of values, attempting to cast them
    /// to the column's expected type.
    ///
    /// # Type Parameters
    ///
    /// * `T`: The type of the elements in the input vector. Must be `Sized` and
    ///   have `'static` lifetime.
    ///
    /// # Arguments
    ///
    /// * `data`: A vector of values to be converted into a Series.
    ///
    /// # Returns
    ///
    /// Returns a `PolarsResult<Series>` containing the created Series or an
    /// error if the conversion fails (e.g., due to incorrect type).
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
