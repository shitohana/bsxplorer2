use std::{io::{BufReader, Read}, sync::Arc};

use itertools::Itertools;
use log::warn;
use noodles::fasta::io::Indexer;
use num::{Float, PrimInt, Unsigned};
use polars::prelude::*;
use paste::paste;

mod stats;
pub use stats::*;

/// Creates a categorical data_structs type from a list of categories.
pub fn get_categorical_dtype(categories: Vec<String>) -> DataType {
    let categories = polars::export::arrow::array::Utf8ViewArray::from_vec(
        categories
            .iter()
            .map(String::as_str)
            .collect_vec(),
        ArrowDataType::Utf8View,
    );
    let rev_mapping = Arc::new(RevMapping::build_local(categories));
    DataType::Enum(Some(rev_mapping), CategoricalOrdering::Physical)
}

/// Creates a schema from separate arrays of names and data_structs types.
pub(crate) fn schema_from_arrays(
    names: &[&str],
    dtypes: &[DataType],
) -> Schema {
    if names.len() != dtypes.len() {
        warn!(
            "Mismatch between names and dtypes array lengths: {} vs {}",
            names.len(),
            dtypes.len()
        );
    }
    Schema::from_iter(
        names
            .iter()
            .cloned()
            .map_into()
            .zip(dtypes.iter().cloned()),
    )
}

/// Creates a hashmap from separate arrays of names and data_structs types.
pub(crate) fn hashmap_from_arrays<'a>(
    names: &[&'a str],
    dtypes: &[DataType],
) -> PlHashMap<&'a str, DataType> {
    if names.len() != dtypes.len() {
        warn!(
            "Mismatch between names and dtypes array lengths: {} vs {}",
            names.len(),
            dtypes.len()
        );
    }
    PlHashMap::from_iter(
        names
            .iter()
            .cloned()
            .map_into()
            .zip(dtypes.iter().cloned()),
    )
}

#[macro_export]
macro_rules! plsmallstr {
    ($string: expr) => {
        PlSmallStr::from($string)
    };
    () => {
        PlSmallStr::from("")
    };
    () => {
        PlSmallStr::from("")
    };
}

#[macro_export]
macro_rules! with_field_fn {
    ($field_name: ident, $field_type: ty) => {
        paste::paste! {
            pub fn [<with_$field_name>](mut self, value: $field_type) -> Self {
            self.$field_name = value;
            self
            }
        }
    };
}

pub fn read_chrs_from_fai<R: Read>(reader: R) -> anyhow::Result<Vec<String>> {
    let records: Vec<noodles::fasta::fai::Record> =
        noodles::fasta::fai::io::Reader::new(BufReader::new(reader)).read_index()?.into();
    Ok(records
        .into_iter()
        .map(|r| String::from_utf8_lossy(r.name()).to_string())
        .collect())
}

pub fn read_chrs_from_fa<R: Read>(reader: R) -> anyhow::Result<Vec<String>> {
    let mut indexer = Indexer::new(BufReader::new(reader));
    let mut records = Vec::new();

    while let Some(record) = indexer.index_record()? {
        records.push(record);
    }

    Ok(records
        .into_iter()
        .map(|r| String::from_utf8_lossy(r.name()).to_string())
        .collect())
}

pub(crate) fn float2int<F: Float, U: Unsigned + PrimInt>(
    value: F
) -> anyhow::Result<U> {
    if value < F::from(0.0).unwrap() || value > F::from(1.0).unwrap() {
        anyhow::bail!("Value should be in [0,1] interval")
    }
    else if value == F::from(0.0).unwrap() {
        Ok(U::min_value())
    }
    else if value == F::from(1.0).unwrap() {
        Ok(U::max_value())
    }
    else {
        Ok(
            U::from((value * F::from(U::max_value()).unwrap()).floor())
                .unwrap(),
        )
    }
}

pub(crate) fn int2float<F: Float, U: Unsigned + PrimInt>(value: U) -> F {
    if value == U::max_value() {
        return F::from(1.0).unwrap();
    }
    else if value == U::min_value() {
        return F::from(0.0).unwrap();
    }
    else {
        F::from(value).unwrap() / F::from(U::max_value()).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(0.0, 0)]
    #[case(1.0, u32::MAX)]
    fn test_float2int_edge_cases(
        #[case] input: f64,
        #[case] expected: u32,
    ) {
        // Test boundary values
        assert_eq!(float2int::<f64, u32>(input).unwrap(), expected);
    }

    #[rstest]
    #[case(0, 0.0)]
    #[case(u32::MAX, 1.0)]
    fn test_int2float_edge_cases(
        #[case] input: u32,
        #[case] expected: f64,
    ) {
        // Test boundary values
        assert_eq!(int2float::<f64, u32>(input), expected);
    }

    #[rstest]
    #[case(-0.1)]
    #[case(1.1)]
    fn test_float2int_errors(#[case] input: f64) {
        // Test out of range values
        assert!(float2int::<f64, u32>(input).is_err());
    }

    #[test]
    fn test_roundtrip_conversion() {
        // Test that converting from float to int and back preserves the
        // relation
        let values = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0];

        for &val in &values {
            let int_val = float2int::<f64, u32>(val).unwrap();
            let float_val = int2float::<f64, u32>(int_val);

            // For exact boundary cases
            if val == 0.0 || val == 1.0 {
                assert_eq!(float_val, val);
            }
            else {
                // For other cases, values should be approximately equal
                // Due to precision limitations, we use approximate equality
                assert_approx_eq!(float_val, val, 0.001);
            }
        }
    }

    #[test]
    fn test_relative_ordering() {
        // Test that the relative ordering is preserved
        let values = [0.1, 0.2, 0.3, 0.4, 0.5];

        for i in 0..values.len() - 1 {
            let int_val1 = float2int::<f64, u32>(values[i]).unwrap();
            let int_val2 = float2int::<f64, u32>(values[i + 1]).unwrap();

            // The ordering should be preserved
            assert!(int_val1 < int_val2);

            // Test the reverse conversion maintains ordering
            let float_val1 = int2float::<f64, u32>(int_val1);
            let float_val2 = int2float::<f64, u32>(int_val2);

            assert!(float_val1 < float_val2);
        }
    }
}
