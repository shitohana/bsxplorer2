use polars::prelude::{
    col, lit, Expr, LazyCsvReader, LazyFileListReader, LazyFrame,
};
use polars_core::prelude::*;
use std::path::Path;

pub struct Annotation {
    raw: LazyFrame,
}

/// Struct of settings for annotation file reading.
#[derive(Default)]
pub struct AnnotationFileConf {
    chr_col: usize,
    start_col: usize,
    end_col: usize,
    id_col: Option<usize>,
    strand_col: Option<usize>,
    type_col: Option<usize>,
    comment_prefix: Option<String>,
    separator: Option<u8>,
    has_header: Option<bool>,
    read_filters: Option<Expr>,
}

impl AnnotationFileConf {
    /// Filter non specified column indices and create a vector
    pub fn col_indices(&self) -> Vec<usize> {
        [
            Some(self.chr_col),
            self.strand_col,
            Some(self.start_col),
            Some(self.end_col),
            self.type_col,
            self.id_col,
        ]
            .iter()
            .flatten()
            .cloned()
            .collect::<Vec<_>>()
    }
    /// Filter unspecified columns and get their names vector
    pub fn col_names(&self) -> Vec<String> {
        [
            Some("chr"),
            if self.strand_col.is_some() {
                Some("strand")
            }
            else {
                None
            },
            Some("start"),
            Some("end"),
            if self.type_col.is_some() {
                Some("type")
            }
            else {
                None
            },
            if self.id_col.is_some() {
                Some("id")
            }
            else {
                None
            },
        ]
            .iter()
            .flatten()
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
    }
}

impl Annotation {
    /// Schema of internal representation of annotation file
    pub fn schema() -> Schema {
        let mut schema = Schema::default();
        schema.insert("chr".into(), DataType::String);
        schema.insert("strand".into(), DataType::String);
        schema.insert("start".into(), DataType::UInt64);
        schema.insert("end".into(), DataType::UInt64);
        schema.insert("type".into(), DataType::String);
        schema.insert("id".into(), DataType::String);
        schema
    }

    pub fn new(mut raw: LazyFrame) -> Annotation {
        let schema = raw
            .collect_schema()
            .expect("schema generation failed");
        if schema.get("chr").is_none()
            || schema.get("start").is_none()
            || schema.get("end").is_none()
        {
            panic!("schema doesn't have chr, start or end field");
        };
        Annotation { raw }
    }

    /// Finish mutating LazyFrame. Consumes self and collects DataFrame.
    pub fn finish(self) -> PolarsResult<DataFrame> {
        self.raw.collect()
    }

    pub fn filter(mut self, predicate: Expr) -> Annotation {
        self.raw = self.raw.filter(predicate);
        self
    }

    pub fn with_columns<E: AsRef<[Expr]>>(mut self, exprs: E) -> Annotation {
        self.raw = self.raw.with_columns(exprs);
        self
    }

    /// Read annotation from custom annotation format.
    //  TODO: add regex for ID col.
    pub fn from_custom(file: &str, configuration: AnnotationFileConf) -> Self {
        assert!(Path::new(file).exists());
        let self_schema = Annotation::schema();
        let mut data = LazyCsvReader::new(file)
            .with_separator(configuration.separator.unwrap_or(b'\t'))
            .with_comment_prefix(Some(PlSmallStr::from(
                configuration
                    .comment_prefix
                    .clone()
                    .unwrap_or("#".to_string()),
            )))
            .with_try_parse_dates(false)
            .with_has_header(
                configuration
                    .has_header
                    .unwrap_or(false),
            )
            .finish()
            .expect("could not open file");

        let schema_cols = data
            .collect_schema()
            .expect("schema generation failed");
        let selector = configuration
            .col_indices()
            .iter()
            .zip(configuration.col_names())
            .map(|(idx, name)| -> Expr {
                let old_name = schema_cols
                    .get_at_index(*idx)
                    .unwrap()
                    .0
                    .to_owned();
                col(old_name)
                    .cast(
                        self_schema
                            .get_field(&name)
                            .unwrap()
                            .dtype,
                    )
                    .alias(name)
            })
            .collect::<Vec<_>>();
        let raw = data.select(&selector).filter(
            configuration
                .read_filters
                .unwrap_or(lit(true)),
        );
        Self { raw }
    }

    pub fn from_gff(file: &str) -> Self {
        Self::from_custom(
            file,
            AnnotationFileConf {
                chr_col: 0,
                start_col: 3,
                end_col: 4,
                id_col: Some(8),
                strand_col: Some(6),
                type_col: Some(2),
                comment_prefix: Some("#".to_string()),
                separator: Some(b'\t'),
                has_header: Some(false),
                read_filters: None,
            },
        )
    }
}
