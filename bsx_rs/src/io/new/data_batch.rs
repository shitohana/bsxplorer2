use crate::io::new::data_schema::DataBatchSchema;
use crate::io::new::utils::{encode_context, GenomicPosition};
use crate::region::RegionCoordinates;
use itertools::Itertools;
use polars::frame::DataFrame;
use polars::prelude::*;

/// Main trait for genomic data. Describes traits to retrieve
/// genomic position information.
///
/// # Restrictions
/// 1. Chromosome ('chr') and position ('position') columns must
///    be present (specified in [DataBatchSchema]).
/// 2. Chromosome and position columns must not have nulls.
/// 3. Chromosome column must be compatible for cast to
///    [DataType::String] and position to [DataType::UInt64]
/// 4. Implementation of [DataBatch] must use [DataBatch::check_self]
///    in their constructor
///
/// If those restrictions are not satisfied, unwrap on [None] value
/// can occur. Otherwise, all methods are safe.
pub(crate) trait DataBatch {
    fn data(&self) -> &DataFrame;
    fn data_mut(&mut self) -> &mut DataFrame;

    fn schema(&self) -> &DataBatchSchema;

    fn schema_mut(&mut self) -> &mut DataBatchSchema;

    fn pos_sort_options() -> SortOptions
    where
        Self: Sized,
    {
        SortOptions::default()
    }

    fn check_self(&mut self)
    where
        Self: Sized,
    {
        match self.schema().check_compatible(&self.data().schema()) {
            Ok(_) => {}
            Err(e) => panic!("{}", e),
        }

        // Check columns non-null
        assert_eq!(
            self.data()
                .column(self.schema().position_col.as_str())
                .unwrap()
                .null_count(),
            0,
            "Position column '{}' contains nulls!",
            self.schema().position_col
        );
        assert_eq!(
            self.data()
                .column(self.schema().chr_col.as_str())
                .unwrap()
                .null_count(),
            0,
            "Chromosome column '{}' contains nulls!",
            self.schema().chr_col
        );

        // Cast
        let hash_map = self.schema().hash_map();
        let colnames = self
            .data()
            .get_column_names_owned()
            .into_iter()
            .map(|s| s.to_string())
            .collect_vec();
        let data = self.data_mut();

        for name in colnames {
            // Cast column in-place
            if let Some(dtype) = hash_map.get(name.as_str()) {
                data.apply(name.as_str(), |c| c.cast(dtype).unwrap())
                    .unwrap_or_else(|_| {
                        panic!("Failed to cast column '{}' to datatype {}", name, dtype)
                    });
            } else {
                // Drop column in-place
                let _ = data.drop_in_place(name.as_str());
            }
        }
    }

    fn chr_col(&self) -> &Series {
        self.data()
            .column(&self.schema().chr_col)
            .expect("Missing 'chr' column. Check if DataBatch schema is valid")
            .as_series()
            .expect("Could not convert 'chr' column to Series")
    }

    fn pos_col(&self) -> &Series {
        self.data()
            .column(&self.schema().position_col)
            .expect("Missing 'pos' column. Check if DataBatch schema is valid")
            .as_series()
            .expect("Could not convert 'pos' column to Series")
    }

    fn first_pos(&self) -> GenomicPosition {
        let chr: String = {
            let first = self.chr_col().first();
            let first = first.as_any_value().cast(&DataType::String);
            first.clone().str_value().to_string()
        };
        let pos: u64 = {
            let first = self.pos_col().first();
            first
                .as_any_value()
                .cast(&DataType::UInt64)
                .try_extract()
                .unwrap()
        };

        GenomicPosition::new(chr, pos)
    }

    fn last_pos(&self) -> GenomicPosition {
        let chr: String = {
            let first = self.chr_col().last();
            let first = first.as_any_value().cast(&DataType::String);
            first.clone().str_value().to_string()
        };
        let pos: u64 = {
            let first = self.pos_col().last();
            first
                .as_any_value()
                .cast(&DataType::UInt64)
                .try_extract()
                .expect("Failed to extract numeric value")
        };

        GenomicPosition::new(chr, pos)
    }

    fn region(&self) -> Option<RegionCoordinates> {
        self.first_pos() >> self.last_pos()
    }

    fn check_position_sorted(&self) -> Option<bool>
    where
        Self: Sized,
    {
        if self.check_chr_unique() {
            self.pos_col()
                .is_sorted(Self::pos_sort_options())
                .expect("Failed to check if position is sorted")
                .into()
        } else {
            None
        }
    }

    fn get_chroms(&self) -> Vec<String> {
        self.chr_col()
            .unique_stable()
            .unwrap()
            .cast(&DataType::String)
            .unwrap()
            .iter()
            .map(|v| v.str_value().to_string())
            .collect_vec()
    }

    fn check_chr_unique(&self) -> bool {
        self.chr_col()
            .unique()
            .expect("Failed to get unique chromosomes")
            .len()
            == 1
    }

    fn extend_data(&mut self, other: &DataFrame) {
        self.data_mut()
            .extend(other)
            .expect("Failed to extend data");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::new::bsx_batch::get_bsx_schema;

    fn get_dummy_schema() -> DataBatchSchema {
        DataBatchSchema {
            schema: Schema::from_iter([
                ("chr".into(), DataType::String),
                ("position".into(), DataType::UInt64),
                ("strand".into(), DataType::String),
            ]),
            chr_col: "chr".into(),
            position_col: "position".into(),
            strand_col: Some("strand".into()),
            need_align: true,
            strand_encoded: false,
            mutate_method: |lf| {
                let bsx_schema = get_bsx_schema();
                let mut _hash_map = bsx_schema.hash_map();
                let hash_map: PlHashMap<&str, DataType> = _hash_map
                    .iter_mut()
                    .map(|(k, v)| (k.as_str(), v.clone()))
                    .collect();
                lf.with_columns([
                    col("position").alias("position"),
                    lit(NULL).alias("context"),
                    lit(NULL).alias("strand"),
                    lit(NULL).alias("count_total"),
                    lit(NULL).alias("count_m"),
                ])
                .select(bsx_schema.col_names().into_iter().map(col).collect_vec())
                .cast(hash_map, true)
            },
        }
    }

    struct DummyBatch {
        schema: DataBatchSchema,
        data: DataFrame,
    }

    impl DataBatch for DummyBatch {
        fn data(&self) -> &DataFrame {
            &self.data
        }
        fn data_mut(&mut self) -> &mut DataFrame {
            &mut self.data
        }
        fn schema(&self) -> &DataBatchSchema {
            &self.schema
        }
        fn schema_mut(&mut self) -> &mut DataBatchSchema {
            &mut self.schema
        }
    }

    impl DummyBatch {
        fn new(data: DataFrame, schema: DataBatchSchema) -> Self {
            let mut new = Self { schema, data };
            new.check_self();
            new
        }
    }

    // Test if constructor casts types correctly
    #[test]
    fn check_new() {
        let test_df = df![
            "chr" => [1, 2],
            "position" => ["1", "2"],
            "strand" => ["-", "+"]
        ]
        .unwrap();
        let schema = get_dummy_schema();
        let batch = DummyBatch::new(test_df.clone(), schema.clone());
        assert_eq!(
            batch
                .data()
                .column(schema.chr_col.as_str())
                .unwrap()
                .dtype(),
            schema.get(schema.chr_col.as_str()).unwrap()
        );
        assert_eq!(
            batch
                .data()
                .column(schema.position_col.as_str())
                .unwrap()
                .dtype(),
            schema.get(schema.position_col.as_str()).unwrap()
        );
        assert_eq!(
            batch
                .data()
                .column(schema.strand_col.clone().unwrap().as_str())
                .unwrap()
                .dtype(),
            schema
                .get(schema.strand_col.clone().unwrap().as_str())
                .unwrap()
        );
    }

    #[test]
    #[should_panic]
    fn check_new_invalid_nulls_chr() {
        let test_df = df![
            "chr" => [Some("1"), None],
            "position" => ["1", "2"],
            "strand" => ["-", "+"]
        ]
        .unwrap();
        let schema = get_dummy_schema();
        let batch = DummyBatch::new(test_df.clone(), schema.clone());
    }
    #[test]
    #[should_panic]
    fn check_new_invalid_nulls_pos() {
        let test_df = df![
            "chr" => ["1", "2"],
            "position" => [Some(1), None],
            "strand" => ["-", "+"]
        ]
        .unwrap();
        let schema = get_dummy_schema();
        let batch = DummyBatch::new(test_df.clone(), schema.clone());
    }

    #[test]
    fn check_pos() {
        let test_df = df![
            "chr" => [1, 2],
            "position" => ["1", "2"],
            "strand" => ["-", "+"]
        ]
        .unwrap();
        let schema = get_dummy_schema();
        let batch = DummyBatch::new(test_df.clone(), schema.clone());
        assert_eq!(batch.first_pos(), GenomicPosition::new("1".into(), 1));
        assert_eq!(batch.last_pos(), GenomicPosition::new("2".into(), 2));
    }

    #[test]
    fn check_region() {
        let test_df = df![
            "chr" => [1, 2],
            "position" => ["1", "2"],
            "strand" => ["-", "+"]
        ]
        .unwrap();
        let schema = get_dummy_schema();
        let batch = DummyBatch::new(test_df.clone(), schema.clone());
        assert_eq!(batch.region(), None);

        let test_df = df![
            "chr" => [1, 1],
            "position" => ["1", "2"],
            "strand" => ["-", "+"]
        ]
        .unwrap();
        let batch = DummyBatch::new(test_df.clone(), schema.clone());
        assert_eq!(
            batch.region(),
            Some(RegionCoordinates::new("1".into(), 1, 2))
        );
    }

    // Check both check_chr_unique and check_position_sorted
    #[test]
    fn check_sorted() {
        let test_df = df![
            "chr" => [1, 2],
            "position" => ["1", "2"],
            "strand" => ["-", "+"]
        ]
        .unwrap();
        let schema = get_dummy_schema();
        let batch = DummyBatch::new(test_df.clone(), schema.clone());
        assert_eq!(batch.check_position_sorted(), None);

        let test_df = df![
            "chr" => [1, 1],
            "position" => ["1", "2"],
            "strand" => ["-", "+"]
        ]
        .unwrap();
        let batch = DummyBatch::new(test_df.clone(), schema.clone());
        assert_eq!(batch.check_position_sorted(), Some(true));

        let test_df = df![
            "chr" => [1, 1],
            "position" => ["3", "2"],
            "strand" => ["-", "+"]
        ]
        .unwrap();
        let batch = DummyBatch::new(test_df.clone(), schema.clone());
        assert_eq!(batch.check_position_sorted(), Some(false));
    }

    // Check right unique chromosomes as string vec and
    // check they maintain order
    #[test]
    fn check_chroms_vals() {
        let test_df = df![
            "chr" => ["1", "2"],
            "position" => ["1", "2"],
            "strand" => ["-", "+"]
        ]
        .unwrap();
        let schema = get_dummy_schema();
        let batch = DummyBatch::new(test_df.clone(), schema.clone());
        assert_eq!(batch.get_chroms(), vec!["1", "2"]);
    }
}
