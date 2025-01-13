use itertools::Itertools;
use polars::datatypes::{CategoricalOrdering, DataType, PlHashMap, PlSmallStr};
use polars::error::PolarsError;
use polars::error::PolarsError::ColumnNotFound;
use polars::prelude::{LazyFrame, PolarsResult, Schema};
use std::error::Error;

pub enum ReportTypeSchema {
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
    Bsx,
    BsxEncoded,
}

impl ReportTypeSchema {
    const fn col_names(&self) -> &[&'static str] {
        match self {
            Self::Bismark => &[
                "chr", "position", "strand", "count_m", "count_um", "context", "trinuc",
            ],
            Self::Coverage => &[
                "chr",
                "nuc",
                "position",
                "context",
                "dinuc",
                "count_m",
                "count_total",
            ],
            Self::CgMap => &[
                "chr",
                "nuc",
                "position",
                "context",
                "dinuc",
                "count_m",
                "count_total",
            ],
            Self::BedGraph => &["chr", "start", "end", "density"],
            Self::Bsx | Self::BsxEncoded => &[
                "chr", "strand", "position", "context", "count_m", "count_um", "density",
            ],
        }
    }

    const fn col_types(&self) -> &[DataType] {
        match self {
            Self::Bismark => &[
                DataType::String,
                DataType::UInt64,
                DataType::String,
                DataType::UInt32,
                DataType::UInt32,
                DataType::String,
                DataType::String,
            ],
            Self::CgMap => &[
                DataType::String,
                DataType::UInt64,
                DataType::UInt64,
                DataType::Float64,
                DataType::UInt32,
                DataType::UInt32,
            ],
            Self::BedGraph => &[
                DataType::String,
                DataType::UInt64,
                DataType::UInt64,
                DataType::Float64,
            ],
            Self::Coverage => &[
                DataType::String,
                DataType::UInt64,
                DataType::UInt64,
                DataType::Float64,
                DataType::UInt32,
                DataType::UInt32,
            ],
            Self::Bsx => &[
                DataType::String,
                DataType::String,
                DataType::UInt64,
                DataType::String,
                DataType::UInt32,
                DataType::UInt32,
                DataType::Float64,
            ],
            Self::BsxEncoded => &[
                DataType::Enum(None, CategoricalOrdering::Physical), // Chr
                DataType::Boolean,                                   // Strand
                DataType::UInt64,
                DataType::Boolean, // Context
                DataType::UInt32,
                DataType::UInt32,
                DataType::Float64,
            ],
        }
    }

    const fn chr_col(&self) -> &'static str {
        match self {
            Self::Bismark
            | Self::BedGraph
            | Self::Coverage
            | Self::Bsx
            | Self::BsxEncoded
            | Self::CgMap => "chr",
        }
    }

    const fn position_col(&self) -> &'static str {
        match self {
            Self::Coverage | Self::BedGraph => "start",

            Self::Bsx | Self::BsxEncoded | Self::Bismark | Self::CgMap => "position",
        }
    }

    const fn context_col(&self) -> Option<&'static str> {
        match self {
            Self::Bsx | Self::BsxEncoded | Self::Bismark | Self::CgMap => Some("context"),

            Self::BedGraph | Self::Coverage => None,
        }
    }

    const fn strand_col(&self) -> Option<&'static str> {
        match self {
            Self::Bsx | Self::BsxEncoded | Self::Bismark | Self::CgMap => Some("strand"),

            Self::BedGraph | Self::Coverage => None,
        }
    }

    pub(crate) const fn need_align(&self) -> bool {
        match self {
            Self::Bsx | Self::BsxEncoded | Self::Bismark | Self::CgMap => false,

            Self::BedGraph | Self::Coverage => true,
        }
    }

    fn schema(&self) -> Schema {
        Schema::from_iter(
            self.col_names()
                .iter()
                .cloned()
                .map_into()
                .zip(self.col_types().iter().cloned()),
        )
    }

    fn hashmap(&self) -> PlHashMap<&str, DataType> {
        PlHashMap::from_iter(
            self.col_names()
                .iter()
                .cloned()
                .map_into()
                .zip(self.col_types().iter().cloned()),
        )
    }
}

#[derive(Debug, Clone)]
pub struct DataBatchSchema {
    /// Target [Schema] for data
    pub(crate) schema: Schema,
    /// Name of the chromosome column
    pub(crate) chr_col: String,
    /// Name of position column
    pub(crate) position_col: String,
    /// Name of strand column (if present)
    pub(crate) strand_col: Option<String>,
    /// Does incoming data contain all possible
    /// cytosines from the reference genome
    pub(crate) need_align: bool,
    /// Is strand encoded as [DataType::Boolean] where
    ///     - "+" => Some(true)
    ///     - "-" => Some(false)
    ///     - "<other>" => None
    pub(crate) strand_encoded: bool,
    /// Method which mutates data described by current
    /// schema into BSX schema ([super::bsx_batch::get_bsx_schema])
    pub(crate) mutate_method: fn(LazyFrame) -> LazyFrame,
}

impl Default for DataBatchSchema {
    fn default() -> Self {
        Self {
            schema: Default::default(),
            chr_col: "chr".to_string(),
            position_col: "position".to_string(),
            strand_col: None,
            need_align: false,
            strand_encoded: false,
            // TODO remove from default batch schema
            mutate_method: |_| LazyFrame::default(),
        }
    }
}

impl DataBatchSchema {
    /// Try to construct [DataBatchSchema] from [Schema] by checking
    /// default names of columns:
    ///
    /// - chromosome - "chr"
    /// - position - "position"
    /// - strand - "strand"
    ///
    /// If either "chr" or "position" are not found, Error is returned
    fn try_new(schema: &Schema) -> PolarsResult<Self> {
        let chr_col = if schema.contains("chr") {
            Some("chr")
        } else {
            None
        };
        let pos_col = if schema.contains("position") {
            Some("position")
        } else {
            None
        };

        if chr_col.is_none() {
            return Err(ColumnNotFound("chr".into()));
        }
        if pos_col.is_none() {
            return Err(ColumnNotFound("position".into()));
        }

        let strand_col = if schema.contains("position") {
            Some("position")
        } else {
            None
        };
        let strand_encoded = if let Some(strand_col) = strand_col {
            matches!(schema.get(strand_col).unwrap(), DataType::Boolean)
        } else {
            false
        };

        let need_align = true;
        let mutate_method = |_| unimplemented!();

        let data_schema = Self {
            schema: schema.clone(),
            chr_col: chr_col.unwrap().to_string(),
            position_col: pos_col.unwrap().to_string(),
            strand_col: strand_col.map(|v| v.to_string()),
            strand_encoded,
            need_align,
            mutate_method,
        };

        data_schema.check_validity()?;
        Ok(data_schema)
    }

    /// Check if other fields of [DataBatchSchema] match schema
    pub(crate) fn check_validity(&self) -> PolarsResult<()> {
        self._check_chr()?;
        self._check_position()?;
        self._check_strand()?;
        Ok(())
    }

    fn _check_chr(&self) -> PolarsResult<()> {
        if self.get(self.chr_col.as_str()).is_none() {
            Err(ColumnNotFound(self.chr_col.clone().into()))
        } else {
            Ok(())
        }
    }

    fn _check_position(&self) -> PolarsResult<()> {
        if self.get(self.position_col.as_str()).is_none() {
            Err(ColumnNotFound(self.position_col.clone().into()))
        } else {
            Ok(())
        }
    }

    fn _check_strand(&self) -> PolarsResult<()> {
        if let Some(strand_col) = &self.strand_col {
            if self.get(strand_col.as_str()).is_none() {
                return Err(ColumnNotFound(strand_col.clone().into()));
            };
            if self.strand_encoded {
                if !self.get(strand_col.as_str()).unwrap().is_bool() {
                    return Err(PolarsError::SchemaMismatch(
                        format!(
                            "Strand column {} is marked as encoded, but DataType is not Boolean!",
                            self.strand_col.clone().unwrap()
                        )
                        .into(),
                    ));
                };
            } else if !self.get(strand_col.as_str()).unwrap().is_string() {
                return Err(PolarsError::SchemaMismatch(
                    format!(
                        "Strand column {} is marked as not encoded, but DataType is not String!",
                        self.strand_col.clone().unwrap()
                    )
                    .into(),
                ));
            }
        }
        Ok(())
    }
    /// Get vector of column names
    pub fn col_names(&self) -> Vec<String> {
        self.schema
            .iter_names_cloned()
            .map(|s| s.to_string())
            .collect_vec()
    }
    /// Get vector of column types
    pub fn col_types(&self) -> Vec<DataType> {
        self.schema.iter_values().cloned().collect_vec()
    }
    pub(crate) fn hash_map(&self) -> PlHashMap<PlSmallStr, DataType> {
        PlHashMap::from_iter(
            self.col_names()
                .into_iter()
                .map(|s| s.into())
                .zip(self.col_types().iter().cloned()),
        )
    }
    fn contains(&self, name: &str) -> bool {
        self.schema.contains(name)
    }
    /// Get [DataType] by column name
    pub fn get(&self, column: &str) -> Option<&DataType> {
        self.schema.get(column)
    }
    /// Returns columns from self if they are not present
    /// in other.
    pub(crate) fn compare_cols(&self, other: &Schema) -> Vec<String>
    where
        Self: Sized,
    {
        let other_cols = other.iter_names().collect_vec();
        self.col_names()
            .iter()
            .filter(|&col| !other_cols.contains(&&PlSmallStr::from(col)))
            .cloned()
            .collect_vec()
    }

    /// Check if other schema has all necessary fields.
    pub(crate) fn check_compatible(&self, other: &Schema) -> Result<(), Box<dyn Error>>
    where
        Self: Sized,
    {
        let missing = self.compare_cols(other);
        if missing.is_empty() {
            Ok(())
        } else {
            Err(Box::from(format!("Missing columns {:?}", missing)))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_schema() {
        let schema = DataBatchSchema {
            schema: Schema::from_iter([
                ("chr".into(), DataType::String),
                ("position".into(), DataType::UInt64),
                ("strand".into(), DataType::String),
            ]),
            strand_col: Some("strand".into()),
            strand_encoded: false,
            ..Default::default()
        };
        assert!(schema.check_validity().is_ok());
    }

    #[test]
    #[should_panic]
    fn invalid_schema_strand() {
        let schema = DataBatchSchema {
            schema: Schema::from_iter([
                ("chr".into(), DataType::String),
                ("position".into(), DataType::UInt64),
                ("strand".into(), DataType::String),
            ]),
            strand_col: Some("strand".into()),
            strand_encoded: true,
            ..Default::default()
        };
        assert!(schema.check_validity().is_ok());
    }
}
