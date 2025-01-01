use crate::utils::types::IPCEncodedEnum;
use itertools::Itertools;
use polars::prelude::*;
use std::fmt::Debug;
use std::ops::{Add, Div};
use crate::io::new::bsx_batch::BsxSchema;

/// Interface trait for [ReportSchema] values. Consists of
/// getters for constant values (see [DataSchemaConst]).
/// Implementation of this trait is just a boilerplate, if 
/// you've implemented [DataSchemaConst].
pub trait DataSchemaInterface {
    /// Get [Schema] instance
    fn schema(&self) -> Arc<Schema> {
        Arc::from(Schema::from_iter(
            self.col_names().iter().cloned()
                .map(|s| PlSmallStr::from_str(s))
                .zip(self.col_types().iter().cloned())
                .collect_vec()
        ))
    }
    /// Get [PlHashMap] instance
    fn hash_map(&self) -> PlHashMap<&str, DataType> {
        PlHashMap::from_iter(
            self.col_names().into_iter().cloned().zip(self.col_types().iter().cloned())
        )
    }
    
    /// Get [DataType] by column name
    fn get(&self, column: &str) -> Option<&DataType> {
        if let Some(index) = self.col_names().iter().position(|c| *c == column) {
            self.col_types().get(index)
        } else { None }
    }

    /// Returns columns from self if they are not present
    /// in other.
    fn compare_cols(&self, other: &Self) -> Vec<&str> where Self: Sized {
        let other_cols = other.col_names();
        self.col_names().iter().cloned()
            .filter(|col| !other_cols.contains(col))
            .collect_vec()
    }
    
    fn check_compatible(&self, other: &Self) -> bool where Self: Sized {
        if !other.col_names().iter().all(|col| self.col_names().contains(col)) {
            return false;
        }
        for colname in other.col_names() {
            if self.get(colname).unwrap().matches_schema_type(other.get(colname).unwrap()).is_err() {
                return false;
            }
        }
        true
    }
    
    fn bsx_schema_mutate(&self, lf: LazyFrame) -> LazyFrame;
    
    fn chr_col(&self) -> &'static str;
    fn position_col(&self) -> &'static str;
    fn strand_col(&self) -> &'static Option<&'static str>;
    fn need_align(&self) -> &'static bool;
    fn strand_encoded(&self) -> &'static bool;
    fn context_encoded(&self) -> &'static bool;
    fn col_names(&self) -> &'static [&'static str];
    fn col_types(&self) -> &'static [DataType];
}

pub(crate) trait DataSchemaConst {
    /// Name of the chromosome column. Obligatory
    const CHR_COL: &'static str = "chr";
    /// Name of the position column. Obligatory
    const POSITION_COL: &'static str = "position";
    /// Name of the strand column
    const STRAND_COL: Option<&'static str> = Some("strand");
    /// Does report contain all possible cytosines
    const NEED_ALIGN: bool = true;
    /// Is strand encoded with boolean keys
    /// - +: Some(true)
    /// - -: Some(false)
    /// - _: None
    const STRAND_ENCODED: bool = false;
    /// Is context encoded with boolean keys:
    /// - CG: Some(True)
    /// - CHG: Some(False)
    /// - CHH: None
    const CONTEXT_ENCODED: bool = false;
    /// List of column names in the schema
    const COL_NAMES: &'static [&'static str];
    /// List of column [DataType] in the schema
    const COL_TYPES: &'static [DataType];
}

struct BismarkSchema;

impl DataSchemaInterface for BismarkSchema {
    fn bsx_schema_mutate(&self, lf: LazyFrame) -> LazyFrame
    where
        Self: Sized
    {
        let bsx_schema = BsxSchema;
        lf
            .with_column(col("count_m").add(col("count_um")).alias("count_total"))
            .with_column(
                col("count_m")
                    .div(col("count_total").cast(DataType::Float64))
                    .alias("density")
            )
            .select(bsx_schema.col_names().into_iter().map(|s| col(*s)).collect_vec())
            .cast(bsx_schema.hash_map(), true)
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

impl DataSchemaConst for BismarkSchema {
    const NEED_ALIGN: bool = false;
    const COL_NAMES: &'static [&'static str] = &[
        "chr", 
        "position", 
        "strand", 
        "count_m", 
        "count_um", 
        "context", 
        "trinuc"
    ];
    const COL_TYPES: &'static [DataType] = &[
        DataType::String,
        DataType::UInt64,
        DataType::String,
        DataType::UInt32,
        DataType::UInt32,
        DataType::String,
        DataType::String,
    ];
}

struct CgMapSchema;

impl DataSchemaInterface for CgMapSchema {
    fn bsx_schema_mutate(&self, lf: LazyFrame) -> LazyFrame
    where
        Self: Sized
    {
        let bsx_schema = BsxSchema;
        lf
            .with_columns([
                when(col("nuc").eq("C")).then(lit("+"))
                    .when(col("nuc").eq("G")).then(lit("-"))
                    .otherwise(lit("."))
                    .alias("strand"),
                
                col("count_m")
                    .div(col("count_total")
                    .cast(DataType::Float64))
                    .alias("density")
                
            ])
            .select(bsx_schema.col_names().into_iter().map(|s| col(*s)).collect_vec())
            .cast(bsx_schema.hash_map(), true)
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

impl DataSchemaConst for CgMapSchema {
    const STRAND_COL: Option<&'static str> = Some("nuc");
    const NEED_ALIGN: bool = false;
    const COL_NAMES: &'static [&'static str] = &[
        "chr",
        "nuc",
        "position",
        "context",
        "dinuc",
        "count_m",
        "count_total",
    ];
    const COL_TYPES: &'static [DataType] = &[
        DataType::String,
        DataType::UInt64,
        DataType::UInt64,
        DataType::Float64,
        DataType::UInt32,
        DataType::UInt32,
    ];
}

struct CoverageSchema;

impl DataSchemaInterface for CoverageSchema {
    fn bsx_schema_mutate(&self, lf: LazyFrame) -> LazyFrame
    where
        Self: Sized
    {
        let bsx_schema = BsxSchema;
        lf
            .with_columns([
                col(Self::POSITION_COL).alias("position"),
                lit(NULL).alias("context"),
                lit(NULL).alias("strand"),
                col("count_m").add(col("count_um")).alias("count_total"),
            ])
            .select(bsx_schema.col_names().into_iter().map(|s| col(*s)).collect_vec())
            .cast(bsx_schema.hash_map(), true)
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

impl DataSchemaConst for CoverageSchema {
    const POSITION_COL: &'static str = "start";
    const COL_NAMES: &'static [&'static str] = &[
        "chr", "start", "end", "density", "count_m", "count_um"
    ];
    const COL_TYPES: &'static [DataType] = &[
        DataType::String,
        DataType::UInt64,
        DataType::UInt64,
        DataType::Float64,
        DataType::UInt32,
        DataType::UInt32,
    ];
}

struct BedGraphSchema;

impl DataSchemaInterface for BedGraphSchema {
    fn bsx_schema_mutate(&self, lf: LazyFrame) -> LazyFrame
    where
        Self: Sized
    {
        let bsx_schema = BsxSchema;
        lf
            .with_columns([
                col(Self::POSITION_COL).alias("position"),
                lit(NULL).alias("context"),
                lit(NULL).alias("strand"),
                lit(NULL).alias("count_total"),
                lit(NULL).alias("count_m"),
            ])
            .select(bsx_schema.col_names().into_iter().map(|s| col(*s)).collect_vec())
            .cast(bsx_schema.hash_map(), true)
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

impl DataSchemaConst for BedGraphSchema {
    const CHR_COL: &'static str = "chr";
    const COL_NAMES: &'static [&'static str] = &[
        "chr", "start", "end", "density"
    ];
    const COL_TYPES: &'static [DataType] = &[
        DataType::String,
        DataType::UInt64,
        DataType::UInt64,
        DataType::Float64
    ];
}

/// Enum with implemented in BSXplorer report schemas
#[derive(Copy, Clone, Debug)]
pub enum ReportSchema {
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
}

impl ReportSchema {
    const BISMARK: BismarkSchema = BismarkSchema;
    const CGMAP: CgMapSchema = CgMapSchema;
    const BEDGRAPH: BedGraphSchema = BedGraphSchema;
    const COVERAGE: CoverageSchema = CoverageSchema;
    
    /// Get the underlying [DataSchemaInterface] for 
    /// the variant.
    pub fn get(&self) -> Box<dyn DataSchemaInterface> {
        match self { 
            ReportSchema::Bismark => Box::from(Self::BISMARK),
            ReportSchema::CgMap => Box::from(Self::CGMAP),
            ReportSchema::BedGraph => Box::from(Self::BEDGRAPH),
            ReportSchema::Coverage => Box::from(Self::COVERAGE),
        }
    }
    
    pub(crate) fn parse_options(&self) -> CsvParseOptions {
        // Default
        let mut options = CsvParseOptions {
            separator: b'\t',
            quote_char: Some(b'"'),
            eol_char: b'\n',
            encoding: CsvEncoding::Utf8,
            null_values: None,
            missing_is_null: true,
            truncate_ragged_lines: false,
            comment_prefix: None,
            try_parse_dates: false,
            decimal_comma: false,
        };

        // This exists because if more report types
        // will be added, they may need different
        // parse options
        match self {
            _ => {}
        };

        options
    }
    
    pub(crate) fn read_options(&self) -> CsvReadOptions {
        let parse_options = self.parse_options();

        let mut read_options = CsvReadOptions {
            path: None,
            rechunk: false,
            n_threads: None,
            low_memory: false,
            n_rows: None,
            row_index: None,
            columns: None,
            projection: None,
            schema: None,
            schema_overwrite: None,
            dtype_overwrite: None,
            parse_options: Arc::from(parse_options),
            has_header: false,
            chunk_size: 0,
            skip_rows: 0,
            skip_rows_after_header: 0,
            infer_schema_length: None,
            raise_if_empty: false,
            ignore_errors: false,
            fields_to_cast: vec![],
        };

        match self {
            Self::BedGraph => {
                read_options.skip_rows = 1;
            }
            _ => {}
        }

        read_options
    }
}


// struct ReportBatch {
//     data: DataFrame,
//     report_type: ReportType,
// }
// 
// impl ReportBatch {
//     pub(crate) fn try_new(data: DataFrame, report_type: ReportType) -> PolarsResult<Self> {
//         let data_schema = data.schema();
//         let target_schema = report_type.schema();
//         
//         let mut need_cast = false;
//         for target_field in target_schema.iter_fields() {
//             
//             if let Some(data_field) = data_schema.get_field(target_field.name().into()) {
//                 let does_match = data_field.dtype.matches_schema_type(&target_field.dtype)?;
//                 if !does_match {need_cast = true;}
//             } else {
//                 return Err(PolarsError::ColumnNotFound(target_field.name().into()));
//             }
//         }
//         
//         let mut new_data = data;
//         if need_cast {
//             new_data = new_data
//                 .lazy().cast(report_type.hash_map(), true)
//                 .collect()?
//         }
//         
//         Ok(Self { data: new_data, report_type })
//     }
//     
//     pub fn get_data(&self) -> &DataFrame { &self.data }
//     pub fn get_data_mut(&mut self) -> &mut DataFrame { &mut self.data }
//     pub fn get_report_type(&self) -> ReportType { self.report_type }
//     
//     pub(crate) fn into_bsx(self) -> PolarsResult<Self> {
//         let new_data = match self.report_type { 
//             ReportType::Bsx => self.data,
// 
//             ReportType::Bismark => {
//                 self.data.lazy()
//                     .with_column(
//                         col("count_m")
//                             .add(col("count_um"))
//                             .alias("count_total")
//                     )
//                     .with_column(
//                         col("count_m")
//                             .div(col("count_total"))
//                             .alias("density")
//                     )
//                     .collect()?
//             }
//             
//             ReportType::CgMap => {
//                 self.data.lazy()
//                     .with_column(
//                         col("count_m")
//                             .div(col("count_total"))
//                             .alias("density")
//                 ).collect()?
//             },
//             
//             ReportType::BedGraph => {
//                 self.data.lazy()
//                     .with_columns([
//                         col("start").alias("position"),
//                         lit(NULL).alias("strand"),
//                         lit(NULL).alias("count_m"),
//                         lit(NULL).alias("count_total"),
//                     ])
//                     .collect()?
//             }
//             
//             ReportType::Coverage => {
//                 self.data.lazy()
//                     .with_columns([
//                         col("start").alias("position"),
//                         lit(NULL).alias("strand"),
//                         col("count_m").add(col("count_um")).alias("count_total"),
//                     ])
//                     .collect()?
//             }
//         };
//         let data = new_data.select_with_schema(
//             ReportType::Bsx.col_names(), 
//             &SchemaRef::from(ReportType::Bsx.schema())
//         )?;
//         
//         Ok(Self {data, report_type: ReportType::Bsx })
//     }
//     // 
//     // fn into_aligned(self, context_data: ContextData) -> PolarsResult<Self>
//     // where
//     //     Self: Sized,
//     // {
//     //     let mut df: DataFrame = self.into_data();
//     //     let context_df = context_data.0;
//     // 
//     //     let mut join_args = JoinArgs::default().with_coalesce(JoinCoalesce::CoalesceColumns);
//     //     join_args.how = JoinType::Left;
//     // 
//     //     for name in ["context", "strand"] {
//     //         match df.drop_in_place(name) {
//     //             Ok(_) => {
//     //                 debug!("Dropping name: {}", name);
//     //             }
//     //             Err(_) => {}
//     //         };
//     //     }
//     // 
//     //     // TODO fix position
//     //     let aligned = context_df.join(
//     //         &df,
//     //         [PlSmallStr::from_str("position")],
//     //         [PlSmallStr::from("position")],
//     //         join_args,
//     //     )?;
//     // 
//     //     Ok(aligned)
//     // }
// }
