use std::borrow::Cow;
use itertools::Itertools;
use polars::prelude::*;
use crate::io::new::report_schema::{DataSchemaConst, DataSchemaInterface, ReportSchema};

pub struct BsxSchema;

impl DataSchemaConst for BsxSchema {
    const NEED_ALIGN: bool = false;
    const COL_NAMES: &'static [&'static str] = &[
        "chr",
        "strand",
        "position",
        "context",
        "count_m",
        "count_um",
        "density",
    ];
    const COL_TYPES: &'static [DataType] = &[
        DataType::String,
        DataType::Boolean,
        DataType::UInt64,
        DataType::Boolean,
        DataType::UInt32,
        DataType::UInt32,
        DataType::Float64,
    ];
}

impl DataSchemaInterface for BsxSchema {
    fn bsx_schema_mutate(&self, lf: LazyFrame) -> LazyFrame
    where
        Self: Sized
    {
        lf
    }

    fn chr_col(&self) -> &'static str {
        &Self::CHR_COL
    }
    fn position_col(&self) -> &'static str {
        &Self::POSITION_COL
    }
    fn strand_col(&self) -> &'static Option<&'static str> {
        &Self::STRAND_COL
    }
    fn need_align(&self) -> &'static bool {
        &Self::NEED_ALIGN
    }
    fn strand_encoded(&self) -> &'static bool {
        &Self::STRAND_ENCODED
    }
    fn context_encoded(&self) -> &'static bool {
        &Self::CONTEXT_ENCODED
    }
    fn col_names(&self) -> &'static [&'static str] {
        &Self::COL_NAMES
    }
    fn col_types(&self) -> &'static [DataType] {
        &Self::COL_TYPES
    }
}

pub struct BsxBatch(DataFrame);

impl BsxBatch {
    const SCHEMA: BsxSchema = BsxSchema;
    
    const CHR_SORT_OPTIONS: SortOptions = SortOptions {
        descending: false,
        nulls_last: false,
        multithreaded: true,
        maintain_order: true
    };

    const POS_SORT_OPTIONS: SortOptions = SortOptions {
        descending: false,
        nulls_last: false,
        multithreaded: true,
        maintain_order: true
    };
    
    pub fn get_data(&self) -> &DataFrame { &self.0 }
    
    pub fn get_data_mut(&mut self) -> &mut DataFrame { &mut self.0 }
    
    /// Check if "position" column sorted in ascending order
    pub fn is_pos_sorted(&self) -> bool {
        self.0.column("position").unwrap()
            .as_materialized_series()
            .is_sorted(Self::POS_SORT_OPTIONS).unwrap()
    }
    /// Check if "chr" column sorted in ascending order
    pub fn is_chr_sorted(&self) -> bool {
        self.0.column("chr").unwrap()
            .as_materialized_series()
            .is_sorted(Self::CHR_SORT_OPTIONS).unwrap()
    }
    
    pub fn partition_by<I>(&self, cols: I, include_key: bool) -> PolarsResult<Vec<BsxBatch>>
    where I: IntoIterator<Item=String> {
        match self.0.partition_by(cols, include_key) {
            Ok(dfs) => Ok(
                dfs.into_iter().map(BsxBatch).collect_vec()
            ),
            Err(e) => Err(e),
        }
    }
    
    pub fn vstack(&self, other: &BsxBatch) ->PolarsResult<Self> {
        Ok(BsxBatch(self.0.vstack(&other.0)?))
    }
    
    pub fn extend(&mut self, other: &BsxBatch) -> PolarsResult<()> {
        self.0.extend(&other.0)
    }
    pub fn filter(self, predicate: Expr) -> PolarsResult<BsxBatch> {
        let new = self.0.lazy().filter(predicate).collect()?;
        Ok(BsxBatch(new))
    }

    // fn cast_data(df: DataFrame) -> PolarsResult<DataFrame> {
    //     let data_schema = df.schema();
    //     let target_schema = self.schema();
    // 
    //     for field in data_schema.iter_fields() {
    //         if let Some(target_field) = target_schema.get_field(field.name().as_str()) {
    // 
    //             if let Err(e) = field.dtype.matches_schema_type(target_field.dtype()) {
    //                 return Err(e);
    //             };
    // 
    //         } else {
    //             return Err(
    //                 PolarsError::ColumnNotFound(Cow::from(field.name().to_string()).into())
    //             );
    //         }
    //     }
    // 
    //     let df = df.select_with_schema(
    //         target_schema.clone().iter_names_cloned(),
    //         &target_schema.into(),
    //     )?;
    // 
    //     Ok(df)
    // }
}