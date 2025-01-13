pub(crate) use report_static::*;

mod report_static {
    use crate::io::new::data_schema::DataBatchSchema;
    use itertools::{izip, Itertools};
    use polars::prelude::*;
    use std::ops::*;

    pub(crate) fn slices_to_schema<'a, K, V>(col_names: K, col_types: V) -> Schema
    where
        K: IntoIterator<Item = &'a &'a str>,
        V: IntoIterator<Item = &'a DataType>,
    {
        Schema::from_iter(izip!(
            col_names.into_iter().map(|s| PlSmallStr::from_str(s)),
            col_types.into_iter().cloned()
        ))
    }

    pub(crate) fn get_bismark() -> DataBatchSchema {
        use bismark::*;

        DataBatchSchema {
            schema: slices_to_schema(COL_NAMES, COL_TYPES),
            need_align: false,
            strand_col: Some("strand".into()),
            mutate_method: bsx_mutate,
            ..Default::default()
        }
    }

    pub(crate) fn get_cgmap() -> DataBatchSchema {
        use cgmap::*;

        DataBatchSchema {
            schema: slices_to_schema(COL_NAMES, COL_TYPES),
            need_align: false,
            strand_col: None,
            mutate_method: bsx_mutate,
            ..Default::default()
        }
    }

    pub(crate) fn get_bedgraph() -> DataBatchSchema {
        use bedgraph::*;

        DataBatchSchema {
            schema: slices_to_schema(COL_NAMES, COL_TYPES),
            position_col: POS_COL.to_string(),
            need_align: true,
            strand_col: None,
            mutate_method: bsx_mutate,
            ..Default::default()
        }
    }

    pub(crate) fn get_coverage() -> DataBatchSchema {
        use coverage::*;

        DataBatchSchema {
            schema: slices_to_schema(COL_NAMES, COL_TYPES),
            position_col: POS_COL.to_string(),
            need_align: true,
            strand_col: None,
            mutate_method: bsx_mutate,
            ..Default::default()
        }
    }

    mod bismark {
        use super::*;
        use crate::io::new::bsx_batch::get_bsx_schema;
        pub(super) const COL_NAMES: &[&str] = &[
            "chr", "position", "strand", "count_m", "count_um", "context", "trinuc",
        ];
        pub(super) const COL_TYPES: &[DataType] = &[
            DataType::String,
            DataType::UInt64,
            DataType::String,
            DataType::UInt32,
            DataType::UInt32,
            DataType::String,
            DataType::String,
        ];
        pub(super) fn bsx_mutate(lf: LazyFrame) -> LazyFrame {
            let bsx_schema = get_bsx_schema();
            let mut _hash_map = bsx_schema.hash_map();
            let hash_map: PlHashMap<&str, DataType> = _hash_map
                .iter_mut()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();

            lf.with_column(col("count_m").add(col("count_um")).alias("count_total"))
                .with_column(
                    col("count_m")
                        .div(col("count_total").cast(DataType::Float64))
                        .alias("density"),
                )
                .select(bsx_schema.col_names().into_iter().map(col).collect_vec())
                .cast(hash_map, true)
        }
    }

    mod cgmap {
        use super::*;
        use crate::io::new::bsx_batch::get_bsx_schema;
        pub(super) const COL_NAMES: &[&str] = &[
            "chr",
            "nuc",
            "position",
            "context",
            "dinuc",
            "count_m",
            "count_total",
        ];
        pub(super) const COL_TYPES: &[DataType] = &[
            DataType::String,
            DataType::UInt64,
            DataType::UInt64,
            DataType::Float64,
            DataType::UInt32,
            DataType::UInt32,
        ];
        pub(super) fn bsx_mutate(lf: LazyFrame) -> LazyFrame {
            let bsx_schema = get_bsx_schema();
            let mut _hash_map = bsx_schema.hash_map();
            let hash_map: PlHashMap<&str, DataType> = _hash_map
                .iter_mut()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();
            lf.with_columns([
                when(col("nuc").eq("C"))
                    .then(lit("+"))
                    .when(col("nuc").eq("G"))
                    .then(lit("-"))
                    .otherwise(lit("."))
                    .alias("strand"),
                col("count_m")
                    .div(col("count_total").cast(DataType::Float64))
                    .alias("density"),
            ])
            .select(bsx_schema.col_names().into_iter().map(col).collect_vec())
            .cast(hash_map, true)
        }
    }

    mod bedgraph {
        use super::*;
        use crate::io::new::bsx_batch::get_bsx_schema;

        pub(super) const COL_NAMES: &[&str] = &["chr", "start", "end", "density"];
        pub(super) const COL_TYPES: &[DataType] = &[
            DataType::String,
            DataType::UInt64,
            DataType::UInt64,
            DataType::Float64,
        ];
        pub(super) const POS_COL: &str = "start";
        pub(super) fn bsx_mutate(lf: LazyFrame) -> LazyFrame {
            let bsx_schema = get_bsx_schema();
            let mut _hash_map = bsx_schema.hash_map();
            let hash_map: PlHashMap<&str, DataType> = _hash_map
                .iter_mut()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();
            lf.with_columns([
                col(POS_COL).alias("position"),
                lit(NULL).alias("context"),
                lit(NULL).alias("strand"),
                lit(NULL).alias("count_total"),
                lit(NULL).alias("count_m"),
            ])
            .select(bsx_schema.col_names().into_iter().map(col).collect_vec())
            .cast(hash_map, true)
        }
    }

    mod coverage {
        use super::*;
        use crate::io::new::bsx_batch::get_bsx_schema;
        pub(super) const COL_NAMES: &[&str] =
            &["chr", "start", "end", "density", "count_m", "count_um"];
        pub(super) const POS_COL: &str = "start";
        pub(super) const COL_TYPES: &[DataType] = &[
            DataType::String,
            DataType::UInt64,
            DataType::UInt64,
            DataType::Float64,
            DataType::UInt32,
            DataType::UInt32,
        ];
        pub(super) fn bsx_mutate(lf: LazyFrame) -> LazyFrame {
            let bsx_schema = get_bsx_schema();
            let mut _hash_map = bsx_schema.hash_map();
            let hash_map: PlHashMap<&str, DataType> = _hash_map
                .iter_mut()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();
            lf.with_columns([
                col(POS_COL).alias("position"),
                lit(NULL).alias("context"),
                lit(NULL).alias("strand"),
                col("count_m").add(col("count_um")).alias("count_total"),
            ])
            .select(bsx_schema.col_names().into_iter().map(col).collect_vec())
            .cast(hash_map, true)
        }
    }
}
