use std::path::Path;

use itertools::Itertools;
use polars::prelude::*;

use crate::utils::array_to_schema;

pub const ANNOTATION_SCHEMA: [(&str, DataType); 6] = [
    ("chr", DataType::String),
    ("strand", DataType::String),
    ("start", DataType::UInt64),
    ("end", DataType::UInt64),
    ("type", DataType::String),
    ("id", DataType::String),
];

/// Struct of settings for annotation file reading.
#[derive(Default, Clone)]
pub struct AnnotationFileConf {
    chr_col:        usize,
    start_col:      usize,
    end_col:        usize,
    id_col:         Option<usize>,
    strand_col:     Option<usize>,
    type_col:       Option<usize>,
    comment_prefix: Option<String>,
    separator:      Option<u8>,
    has_header:     Option<bool>,
    read_filters:   Option<Expr>,
    regex:          Option<Expr>,
}

impl AnnotationFileConf {
    fn all_colnames() -> Vec<String> {
        AnnotationBuilder::schema()
            .iter_names()
            .map(|name| name.to_owned().into_string())
            .collect()
    }

    fn indexes(&self) -> [Option<usize>; 6] {
        [
            Some(self.chr_col),
            self.strand_col,
            Some(self.start_col),
            Some(self.end_col),
            self.type_col,
            self.id_col,
        ]
    }

    /// Filter non specified column indices and create a vector
    pub fn col_indices(&self) -> Vec<usize> {
        self.indexes()
            .iter()
            .flatten()
            .cloned()
            .collect::<Vec<_>>()
    }

    /// Filter unspecified columns and get their names vector
    pub fn col_names(&self) -> Vec<String> {
        let all = Self::all_colnames();
        self.indexes()
            .iter()
            .positions(|val| val.is_some())
            .map(|i| all[i].clone())
            .collect()
    }

    pub fn into_reader(
        self,
        file: String,
    ) -> LazyCsvReader {
        LazyCsvReader::new(file)
            .with_separator(self.separator.unwrap_or(b'\t'))
            .with_comment_prefix(Some(PlSmallStr::from(
                self.comment_prefix
                    .clone()
                    .unwrap_or("#".to_string()),
            )))
            .with_has_header(self.has_header.unwrap_or(false))
    }
}

pub struct AnnotationBuilder {
    raw: LazyFrame,
}

impl AnnotationBuilder {
    /// Schema of internal representation of annotation file
    pub fn schema() -> Schema { array_to_schema(&ANNOTATION_SCHEMA) }

    pub fn new(mut raw: LazyFrame) -> AnnotationBuilder {
        let schema = raw
            .collect_schema()
            .expect("schema generation failed");
        if schema.get("chr").is_none()
            || schema.get("start").is_none()
            || schema.get("end").is_none()
        {
            panic!("schema doesn't have chr, start or end field");
        };
        AnnotationBuilder { raw }
    }

    /// Finish mutating LazyFrame. Consumes self and collects DataFrame.
    pub fn finish(self) -> PolarsResult<DataFrame> { self.raw.collect() }

    pub fn filter(
        mut self,
        predicate: Expr,
    ) -> AnnotationBuilder {
        self.raw = self.raw.filter(predicate);
        self
    }

    pub fn with_columns<E: AsRef<[Expr]>>(
        mut self,
        exprs: E,
    ) -> AnnotationBuilder {
        self.raw = self.raw.with_columns(exprs);
        self
    }

    /// Read annotation from custom annotation format.
    pub fn from_custom(
        file: &str,
        configuration: AnnotationFileConf,
    ) -> Self {
        assert!(Path::new(file).exists());
        let self_schema = Self::schema();

        let mut raw = configuration
            .clone()
            .into_reader(file.to_string())
            .finish()
            .expect("could not open file");

        let schema_cols = raw
            .collect_schema()
            .expect("schema generation failed")
            .iter_names()
            .cloned()
            .map(|x| x.into_string())
            .collect::<Vec<_>>();

        let raw = {
            // Rearrange columns
            let mut df = raw.select({
                itertools::izip!(
                    &configuration.col_indices(),
                    &configuration.col_names()
                )
                .map(|(idx, name)| -> Expr {
                    let old_name = schema_cols.get(*idx).unwrap().clone();
                    let target_type = {
                        self_schema
                            .get_field(name)
                            .unwrap()
                            .dtype
                    };
                    let expr = col(old_name)
                        .cast(target_type)
                        .alias(name);
                    expr
                })
                .collect::<Vec<_>>()
            });

            if configuration.read_filters.is_some() {
                // Apply filters if present
                df = df.filter(
                    configuration
                        .read_filters
                        .unwrap_or(lit(true)),
                )
            }

            if configuration.regex.is_some() {
                df = df.with_column(
                    col("id")
                        .str()
                        .extract(configuration.regex.unwrap(), 1),
                )
            }

            df
        };

        Self { raw }
    }

    pub fn from_gff(file: &str) -> Self {
        Self::from_custom(file, AnnotationFileConf {
            chr_col:        0,
            start_col:      3,
            end_col:        4,
            id_col:         Some(8),
            strand_col:     Some(6),
            type_col:       Some(2),
            comment_prefix: Some("#".to_string()),
            separator:      Some(b'\t'),
            has_header:     Some(false),
            read_filters:   None,
            regex:          Some(lit("^ID=([^;]+)")),
        })
    }
}
