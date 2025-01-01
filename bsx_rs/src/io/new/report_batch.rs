use std::fmt::Debug;
use std::ops::Deref;
use polars::frame::DataFrame;
use polars::prelude::*;
use crate::io::new::report_schema::{DataSchemaInterface, ReportSchema};
use crate::io::new::utils::GenomicPosition;
use crate::region::RegionCoordinates;

pub(crate) trait DataBatch
where Self: Sized + Clone + Debug + Into<DataFrame> {
    fn data(&self) -> &DataFrame;
    fn data_mut(&mut self) -> &mut DataFrame;
    fn schema(&self) -> Box<dyn DataSchemaInterface>;
    
    fn chr_col(&self) -> &Series {
        self.data().column(self.schema().chr_col()).unwrap().as_series().unwrap()
    }

    fn pos_col(&self) -> &Series {
        self.data().column(self.schema().position_col()).unwrap().as_series().unwrap()
    }

    fn strand_col(&self) -> Option<&Series> {
        if let Some(colname) = self.schema().strand_col() {
            Some(self.data().column(colname).unwrap().as_series().unwrap())
        } else { None }
    }
    
    fn first_pos(&self) -> GenomicPosition {
        GenomicPosition::new(
            self.chr_col().first().as_any_value().str_value().into(),
            self.pos_col().first().as_any_value().into_static().try_extract().unwrap(),
        )
    }
    fn last_pos(&self) -> GenomicPosition {
        GenomicPosition::new(
            self.chr_col().last().as_any_value().str_value().into(),
            self.pos_col().last().as_any_value().into_static().try_extract().unwrap(),
        )
    }
    
    fn region(&self) -> Option<RegionCoordinates> {
        self.last_pos() - self.first_pos()
    }

    fn check_position_sorted(&self) -> PolarsResult<bool> {
        self.pos_col().is_sorted(SortOptions::default())
    }
    fn check_chr_unique(&self) -> PolarsResult<usize> {
        Ok(self.chr_col().unique()?.len())
    }
    
    #[cfg(not(all()))]
    // This method casts one ReportBatch to another.
    // But as we can not determine constant transformation methods
    // This method is unneeded
    fn into_batch_with_schema(self, schema: ReportSchema) -> PolarsResult<Self>
    where Self: Sized {
        let from_schema = self.schema().get();
        let target_schema = schema.get();

        let mut lf = DataFrame::from(self).lazy();

        // Rename columns
        lf = lf.rename(from_schema.chr_col(), target_schema.chr_col(), true);
        lf = lf.rename(from_schema.position_col(), target_schema.position_col(), true);
        if target_schema.strand_col() != from_schema.strand_col() {
            lf = lf.rename(from_schema.strand_col().unwrap(), target_schema.strand_col().unwrap(), true);
        };
        
        // Create missing columns
        target_schema.compare_cols(&*from_schema).iter()
            .for_each(|colname| lf = lf.with_column(lit(NULL).alias(colname)));
        // Remove redundant columns
        lf = lf.drop(from_schema.compare_cols(&*target_schema));

        // Encode/decode context
        if from_schema.get("context").is_some() && target_schema.get("context").is_some() {
            match (from_schema.context_encoded(), target_schema.context_encoded()) {
                (true, false) => { lf = decode_context(lf, "context") },
                (false, true) => { lf = encode_context(lf, "context") },
                _ => {}
            }
        }
        // Encode/decode strand
        if from_schema.strand_col().is_some() && target_schema.strand_col().is_some() {
            match (from_schema.strand_encoded(), target_schema.strand_encoded()) {
                (true, false) => { lf = decode_strand(lf, target_schema.strand_col().unwrap()) },
                (false, true) => { lf = encode_strand(lf, target_schema.strand_col().unwrap()) },
                _ => {}
            }
        }
        
        // Final cast
        lf = lf.select(target_schema.col_names().iter().map(col));
        lf = lf.cast(target_schema.hash_map(), true);
        
        lf.collect()
    }
}
