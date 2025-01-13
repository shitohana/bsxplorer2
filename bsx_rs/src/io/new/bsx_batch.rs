use crate::io::new::data_schema::DataBatchSchema;
use crate::io::new::report_schema::slices_to_schema;
use polars::prelude::*;

const COL_NAMES: &[&str] = &[
    "chr", "strand", "position", "context", "count_m", "count_um", "density",
];
const COL_TYPES: &[DataType] = &[
    DataType::String,
    DataType::Boolean,
    DataType::UInt64,
    DataType::Boolean,
    DataType::UInt32,
    DataType::UInt32,
    DataType::Float64,
];

pub(crate) fn get_bsx_schema() -> DataBatchSchema {
    DataBatchSchema {
        schema: slices_to_schema(COL_NAMES, COL_TYPES),
        need_align: false,
        strand_col: Some("strand".into()),
        mutate_method: |lf| lf,
        ..Default::default()
    }
}

#[derive(Debug, Clone)]
pub struct BsxBatch {
    schema: DataBatchSchema,
    data: DataFrame,
}

impl BsxBatch {
    const CHR_SORT_OPTIONS: SortOptions = SortOptions {
        descending: false,
        nulls_last: false,
        multithreaded: true,
        maintain_order: true,
        limit: None,
    };

    const POS_SORT_OPTIONS: SortOptions = SortOptions {
        descending: false,
        nulls_last: false,
        multithreaded: true,
        maintain_order: true,
        limit: None,
    };
}

impl BsxBatch {
    // pub fn try_new(data: DataFrame) -> PolarsResult<Self> {
    //     let casted = data.select_with_schema(
    //         Self::SCHEMA.col_names().iter().copied(),
    //         &SchemaRef::from(Self::SCHEMA.schema())
    //     )?;
    //     Ok(BsxBatch(casted))
    // }
    //
    // pub fn partition_by<I>(&self, cols: I, include_key: bool) -> PolarsResult<Vec<BsxBatch>>
    // where I: IntoIterator<Item=String> {
    //     match self.0.partition_by(cols, include_key) {
    //         Ok(dfs) => Ok(
    //             dfs.into_iter().map(BsxBatch).collect_vec()
    //         ),
    //         Err(e) => Err(e),
    //     }
    // }
    //
    // pub fn vstack(&self, other: &BsxBatch) ->PolarsResult<Self> {
    //     Ok(BsxBatch(self.0.vstack(&other.0)?))
    // }
    //
    // pub fn extend(&mut self, other: &BsxBatch) -> PolarsResult<()> {
    //     self.0.extend(&other.0)
    // }
    //
    // pub fn filter(self, predicate: Expr) -> PolarsResult<BsxBatch> {
    //     let new = self.0.lazy().filter(predicate).collect()?;
    //     Ok(BsxBatch(new))
    // }

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
//
// impl DataBatch for BsxBatch {
//     fn data(&self) -> &DataFrame {
//         &self.0
//     }
//
//     fn data_mut(&mut self) -> &mut DataFrame {
//         &mut self.0
//     }
//
//     fn schema(&self) -> Box<dyn DataSchemaInterface> {
//         Box::from(Self::SCHEMA)
//     }
// }

impl From<BsxBatch> for DataFrame {
    fn from(value: BsxBatch) -> Self {
        value.data
    }
}
