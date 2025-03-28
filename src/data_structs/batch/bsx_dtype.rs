use itertools::Itertools;
use polars::prelude::*;
use anyhow::{Result, anyhow};

pub trait BsxPolarsDtype: Clone + Sized + ToString + PartialEq {
    type RustType;
    type PolarsType: PolarsDataType;

    fn from_polars(polars_val: Self::PolarsType) -> Result<Option<Self::RustType>>;
    fn to_polars(val: Self::RustType) -> Result<Self::PolarsType>;
}

// impl BsxPolarsDtype for String {
//     type Output = Self;

//     fn convert_value(val: AnyValue) -> Result<Option<Self::Output>> {
//         match val {
//             AnyValue::Null => Ok(None),
//             AnyValue::String(v) => Ok(Some(v.to_string())),
//             _ => Err(anyhow!("Unsupported value type for String")),
//         }
//     }

//     fn convert_column(column: Column) -> Result<Vec<Option<Self::Output>>> {
//         column
//             .str()
//             .map(|arr| arr.into_iter().map(|v| v.map(|s| s.to_string())).collect_vec())
//             .map_err(|e| anyhow!(e)) // Simplified error handling
//     }
// }

// impl BsxPolarsDtype for u64 {
//     type Output = Self;

//     fn convert_value(val: AnyValue) -> Result<Option<Self::Output>> {
//         match val {
//             AnyValue::Null => Ok(None),
//             AnyValue::UInt64(v) => Ok(Some(v)),
//             _ => Err(anyhow!("Unsupported value type for u64")),
//         }
//     }

//     fn convert_column(column: Column) -> Result<Vec<Option<Self::Output>>> {
//         column
//             .u64()
//             .map(|arr| arr.into_iter().collect_vec())
//             .map_err(|e| anyhow!(e)) // Simplified error handling
//     }
// }

// impl BsxPolarsDtype for f64 {
//     type Output = Self;

//     fn convert_value(val: AnyValue) -> Result<Option<Self::Output>> {
//         match val {
//             AnyValue::Null => Ok(None),
//             AnyValue::Float64(v) => Ok(Some(v)),
//             _ => Err(anyhow!("Unsupported value type for f64")),
//         }
//     }

//     fn convert_column(column: Column) -> Result<Vec<Option<Self::Output>>> {
//         column
//             .f64()
//             .map(|arr| arr.into_iter().collect_vec())
//             .map_err(|e| anyhow!(e))
//     }
// }

// impl BsxPolarsDtype for f32 {
//     type Output = Self;

//     fn convert_value(val: AnyValue) -> Result<Option<Self::Output>> {
//         match val {
//             AnyValue::Null => Ok(None),
//             AnyValue::Float32(v) => Ok(Some(v)),
//             _ => Err(anyhow!("Unsupported value type for f32")),
//         }
//     }

//     fn convert_column(column: Column) -> Result<Vec<Option<Self::Output>>> {
//         column
//             .f32()
//             .map(|arr| arr.into_iter().collect_vec())
//             .map_err(|e| anyhow!(e))
//     }
// }

// impl BsxPolarsDtype for i16 {
//     type Output = Self;

//     fn convert_value(val: AnyValue) -> Result<Option<Self::Output>> {
//         match val {
//             AnyValue::Null => Ok(None),
//             AnyValue::Int16(v) => Ok(Some(v)),
//             _ => Err(anyhow!("Unsupported value type for i16")),
//         }
//     }

//     fn convert_column(column: Column) -> Result<Vec<Option<Self::Output>>> {
//         column
//             .i16()
//             .map(|arr| arr.into_iter().collect_vec())
//             .map_err(|e| anyhow!(e))
//     }
// }

// impl BsxPolarsDtype for u16 {
//     type Output = Self;

//     fn convert_value(val: AnyValue) -> Result<Option<Self::Output>> {
//         match val {
//             AnyValue::Null => Ok(None),
//             AnyValue::UInt16(v) => Ok(Some(v)),
//             _ => Err(anyhow!("Unsupported value type for u16")),
//         }
//     }

//     fn convert_column(column: Column) -> Result<Vec<Option<Self::Output>>> {
//         column
//             .u16()
//             .map(|arr| arr.into_iter().collect_vec())
//             .map_err(|e| anyhow!(e))
//     }
// }

// impl BsxPolarsDtype for u32 {
//     type Output = Self;

//     fn convert_value(val: AnyValue) -> Result<Option<Self::Output>> {
//         match val {
//             AnyValue::Null => Ok(None),
//             AnyValue::UInt32(v) => Ok(Some(v)),
//             _ => Err(anyhow!("Unsupported value type for u32")),
//         }
//     }

//     fn convert_column(column: Column) -> Result<Vec<Option<Self::Output>>> {
//         column
//             .u32()
//             .map(|arr| arr.into_iter().collect_vec())
//             .map_err(|e| anyhow!(e))
//     }
// }

// impl BsxPolarsDtype for bool {
//     type Output = Self;

//     fn convert_value(val: AnyValue) -> Result<Option<Self::Output>> {
//         match val {
//             AnyValue::Null => Ok(None),
//             AnyValue::Boolean(v) => Ok(Some(v)),
//             _ => Err(anyhow!("Unsupported value type for bool")),
//         }
//     }

//     fn convert_column(column: Column) -> Result<Vec<Option<Self::Output>>> {
//         column
//             .bool()
//             .map(|arr| arr.into_iter().collect_vec())
//             .map_err(|e| anyhow!(e))
//     }
// }
