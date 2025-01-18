use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::BitAnd;
use crate::region::{GenomicPosition, RegionCoordinates};
use crate::utils::types::{IPCEncodedEnum, Strand, Context};
use itertools::Itertools;
use num::{PrimInt, Unsigned};
use polars::prelude::*;
use crate::io::report::report_batch_utils::{decode_context, decode_strand, encode_context, encode_strand};

#[derive(Debug)]
/// DataFrame with
/// 1. Non-null chr and position
/// 2. Single unique chr
/// 3. Ascending sorted positions
#[derive(Clone)]
#[derive(PartialEq)]
pub struct BsxBatch(DataFrame);

impl BsxBatch {
    pub(crate) unsafe fn new_unchecked(data_frame: DataFrame) -> Self {
        BsxBatch(data_frame)
    }


    pub(crate) const fn col_names() -> &'static [&'static str] {
        &[
            "chr",
            "position",
            "strand",
            "context",
            "count_m",
            "count_total",
            "density",
        ]
    }

    pub(crate) const fn col_types() -> &'static [DataType] {
        &[
            DataType::String,
            DataType::UInt64,
            DataType::String,
            DataType::String,
            DataType::UInt32,
            DataType::UInt32,
            DataType::Float64,
        ]
    }

    pub(crate) const fn chr_col() -> &'static str {
        "chr"
    }

    pub(crate) const fn pos_col() -> &'static str {
        "position"
    }

    pub fn schema() -> Schema {
        use crate::io::report::report_batch_utils::schema_from_arrays;
        schema_from_arrays(Self::col_names(), Self::col_types())
    }

    pub fn hashmap() -> PlHashMap<&'static str, DataType> {
        use crate::io::report::report_batch_utils::hashmap_from_arrays;
        hashmap_from_arrays(Self::col_names(), Self::col_types())
    }
}

impl BsxBatchMethods for BsxBatch {
    fn filter(self, context: Option<Context>, strand: Option<Strand>) -> Self
    where
        Self: Sized
    {
        if context.is_none() && strand.is_none() {
            self
        } else {
            let mut batch_lazy = DataFrame::from(self).lazy();
            if let Some(strand) = strand {
                batch_lazy = batch_lazy.filter(col("strand").eq(lit(strand.to_string())));
            }
            if let Some(context) = context {
                batch_lazy = batch_lazy.filter(col("context").eq(lit(context.to_string())));
            }
            Self(batch_lazy.collect().unwrap())
        }
    }

    fn data(&self) -> &DataFrame {
        &self.0
    }

    fn data_mut(&mut self) -> &mut DataFrame {
        &mut self.0
    }
}

impl TryFrom<DataFrame> for BsxBatch {
    type Error = PolarsError;

    fn try_from(value: DataFrame) -> Result<Self, Self::Error> {
        // Check empty
        if value.is_empty() {
            return Err(PolarsError::NoData("Data is empty".into()));
        }
        // Try cast
        let data_casted = value.lazy().cast(Self::hashmap(), true).collect()?;
        // Chr null check
        if data_casted.column(Self::chr_col())?.null_count() != 0 {
            return Err(PolarsError::ComputeError("Nulls in chr col".into()));
        }
        // Pos null check
        if data_casted.column(Self::pos_col())?.null_count() != 0 {
            return Err(PolarsError::ComputeError("Nulls in pos col".into()));
        }
        // Check sorted
        if !data_casted
            .column(Self::pos_col())?
            .as_series()
            .unwrap()
            .is_sorted(SortOptions::default().with_order_descending(false))?
        {
            return Err(PolarsError::ComputeError("Data is not sorted".into()));
        }

        Ok(BsxBatch(data_casted))
    }
}


impl From<BsxBatch> for DataFrame {
    fn from(b: BsxBatch) -> Self {
        b.0
    }
}   

#[derive(Debug, Clone, PartialEq)]
pub struct EncodedBsxBatch(DataFrame);

impl EncodedBsxBatch {
    pub(crate) unsafe fn new_unchecked(data_frame: DataFrame) -> Self {
        EncodedBsxBatch(data_frame)
    }
    pub fn encode(batch: BsxBatch, chr_dtype: &DataType) -> PolarsResult<Self> {
        let batch_data = DataFrame::from(batch);
        let target_schema = Self::get_hashmap(chr_dtype);
        
        let mut batch_lazy = batch_data.lazy();
        batch_lazy = encode_context(batch_lazy, "context");
        batch_lazy = encode_strand(batch_lazy, "strand");
        
        batch_lazy = batch_lazy
            .cast(target_schema, true)
            .select(Self::col_names().iter().cloned().map(col).collect_vec());
        let mut result = batch_lazy.collect()?;
        result.rechunk_mut();
        
        Ok(Self(result))
    }
    pub fn decode(self) -> PolarsResult<BsxBatch> {
        let batch_data = DataFrame::from(self);
        let target_schema = BsxBatch::hashmap();
        let mut batch_lazy = batch_data.lazy();
        batch_lazy = decode_context(batch_lazy, "context", "context");
        batch_lazy = decode_strand(batch_lazy, "strand", "strand");
        batch_lazy = batch_lazy
            .cast(target_schema, true)
            .select(BsxBatch::col_names().iter().cloned().map(col).collect_vec());

        let result = batch_lazy.collect()?;
        Ok(BsxBatch(result))
    }
    
    pub(crate) fn get_schema(chr_dtype: &DataType) -> Schema {
        Schema::from_iter([
            ("chr".into(), chr_dtype.clone()),
            ("position".into(), DataType::UInt32),
            ("strand".into(), DataType::Boolean),
            ("context".into(), DataType::Boolean),
            ("count_m".into(), DataType::Int16),
            ("count_total".into(), DataType::Int16),
            ("density".into(), DataType::Float32),
        ])
    }

    pub(crate) fn get_hashmap(chr_dtype: &DataType) -> PlHashMap<&str, DataType> {
        PlHashMap::from_iter([
            ("chr", chr_dtype.clone()),
            ("position", DataType::UInt32),
            ("strand", DataType::Boolean),
            ("context", DataType::Boolean),
            ("count_m", DataType::Int16),
            ("count_total", DataType::Int16),
            ("density", DataType::Float32)
        ])
    }

    pub fn trim_region<N: PrimInt + Unsigned + Display>(&self, region_coordinates: &RegionCoordinates<u64>) -> PolarsResult<Self> where Self: Sized {
        let batch_first = self.first_position()?;
        let pos_col = self.data().column("position")?.u32()?;
        let height = self.height();

        match batch_first.partial_cmp(&region_coordinates.end_gpos()) {
            None => Err(PolarsError::ComputeError("Chromosome does not match".into())),
            Some(Ordering::Greater) => Err(PolarsError::ComputeError("Batch does not contain region information".into())),
            _ => {
                let start_shift = pos_col.iter().position(|v| v.unwrap() > region_coordinates.start() as u32).unwrap().saturating_sub(1);
                let end_shift = pos_col.iter().skip(start_shift).position(|v| v.unwrap() > region_coordinates.end() as u32).map(|val| val + start_shift).unwrap_or(height);
                let mask = BooleanChunked::from_iter({
                    let mut start = vec![false; start_shift];
                    start.extend(vec![true; end_shift - start_shift]);
                    start.extend(vec![false; height - end_shift]);
                    start
                });
                self.data().filter(&mask).map(|df| unsafe {Self::new_unchecked(df)} )
            }
        }
    }
    
    pub(crate) fn schema(&self) -> Schema {
        self.data().schema()
    }
    pub(crate) fn col_names() -> &'static [&'static str] {
        BsxBatch::col_names()
    }
}

impl BsxBatchMethods for EncodedBsxBatch {
    fn filter(self, context: Option<Context>, strand: Option<Strand>) -> Self where Self: Sized {
        if context.is_none() && strand.is_none() {
            self
        } else {
            let mut boolean_chunked = BooleanChunked::new("mask".into(), vec![false; self.0.height()]);
            if let Some(strand) = strand {
                let mask = BooleanChunked::new(
                    "strand".into(),
                    self.0
                        .column("strand")
                        .unwrap()
                        .bool()
                        .unwrap()
                        .into_iter()
                        .map(|v| v == strand.to_bool())
                        .collect_vec(),
                );
                boolean_chunked = boolean_chunked.bitand(mask);
            }
            if let Some(context) = context {
                let mask = BooleanChunked::new(
                    "context".into(),
                    self.0
                        .column("context")
                        .unwrap()
                        .bool()
                        .unwrap()
                        .into_iter()
                        .map(|v| v == context.to_bool())
                        .collect_vec(),
                );
                boolean_chunked = boolean_chunked.bitand(mask);
            }
            Self(self.0.filter(&boolean_chunked).unwrap())
        }
    }

    fn data(&self) -> &DataFrame {
        &self.0
    }
    fn data_mut(&mut self) -> &mut DataFrame {
        &mut self.0
    }
    
    fn first_position(&self) -> Result<GenomicPosition<u64>, PolarsError> {
        let chr_col = self.data().column("chr")?.categorical()?;
        let chr = chr_col.get_rev_map().get(chr_col.physical().first().unwrap());
        let pos = self.data().column("position")?.u32()?.first().unwrap();
        Ok(GenomicPosition::new(chr.to_string(), pos as u64))
    }

    fn last_position(&self) -> Result<GenomicPosition<u64>, PolarsError> {
        let chr_col = self.data().column("chr")?.categorical()?;
        let chr = chr_col.get_rev_map().get(chr_col.physical().last().unwrap());
        let pos = self.data().column("position")?.u32()?.last().unwrap();
        Ok(GenomicPosition::new(chr.to_string(), pos as u64))
    }
}

impl From<EncodedBsxBatch> for DataFrame {
    fn from(bs: EncodedBsxBatch) -> Self {
        bs.0
    }
}

pub trait BsxBatchMethods {
    fn filter(self, context: Option<Context>, strand: Option<Strand>) -> Self where Self: Sized;
    
    fn data(&self) -> &DataFrame;
    
    fn data_mut(&mut self) -> &mut DataFrame;
    
    fn check_chr_unique(&self) -> bool {
        self.data().column("chr").unwrap().as_series().unwrap().unique().unwrap().len() == 1
    }
    
    fn first_position(&self) -> PolarsResult<GenomicPosition<u64>> {
        use crate::io::report::report_batch_utils::first_position;
        first_position(self.data(), BsxBatch::chr_col(), BsxBatch::pos_col())
    }
    fn last_position(&self) -> PolarsResult<GenomicPosition<u64>> {
        use crate::io::report::report_batch_utils::last_position;
        last_position(self.data(), BsxBatch::chr_col(), BsxBatch::pos_col())
    }
    fn extend(&mut self, other: &Self) -> PolarsResult<()> where Self: Sized {
        if !(self.check_chr_unique() && other.check_chr_unique()) {
            return Err(PolarsError::ComputeError("Chromosomes in batches non-unique".into()));
        }
        let self_pos = self.last_position()?;
        let other_pos = self.first_position()?;
        
        if self_pos.chr() != other_pos.chr() {
            return Err(PolarsError::ComputeError("Chromosomes in batches differ".into()))
        }
        if self_pos.position() >= other_pos.position() {
            return Err(PolarsError::ComputeError("First position in other batch must be less than last position in the current".into()))
        }
        
        self.data_mut().extend(other.data())
    }
    
    fn height(&self) -> usize {
        self.data().height()
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::get_categorical_dtype;
    use super::*;
    
    fn dummy_batch() -> BsxBatch {
        BsxBatch::try_from(df![
                "chr" => ["1", "1", "1"],
                "position" => [1, 2, 3],
                "strand" => [".", ".", "."],
                "context" => ["CG", "CHG", "CHH"],
                "count_m" => [None::<u32>, None, None],
                "count_total" => [None::<u32>, None, None],
                "density" => [0, 1, 0]
            ].unwrap().lazy()
            .cast(BsxBatch::hashmap(), true)
            .collect().unwrap()
        ).unwrap()
    }
    
    #[test]
    fn test_chr_unique() {
        let batch = BsxBatch::try_from(df![
                "chr" => ["1", "2", "3"],
                "position" => [1, 2, 3],
                "strand" => [".", ".", "."],
                "context" => [None::<&str>, None, None],
                "count_m" => [None::<u32>, None, None],
                "count_total" => [None::<u32>, None, None],
                "density" => [0, 1, 0]
            ].unwrap().lazy()
                .cast(BsxBatch::hashmap(), true)
                .collect().unwrap()
        ).unwrap();
        
        assert!(!batch.check_chr_unique());

        let batch = dummy_batch();
        
        assert!(batch.check_chr_unique());
        
        let chr_dtype = get_categorical_dtype(vec!["1".into()]);
        let encoded = EncodedBsxBatch::encode(batch, &chr_dtype).unwrap();
        
        assert!(encoded.check_chr_unique())
    }
    
    #[test]
    fn test_position() {
        let batch = dummy_batch();
        let chr_dtype = get_categorical_dtype(vec!["1".into()]);
        let encoded = EncodedBsxBatch::encode(batch.clone(), &chr_dtype).unwrap();
        
        assert_eq!(batch.first_position().unwrap(), encoded.first_position().unwrap());
        assert_eq!(batch.last_position().unwrap(), encoded.last_position().unwrap());
    }
    
    #[test]
    fn test_encode_decode() {
        let batch = dummy_batch();
        let chr_dtype = get_categorical_dtype(vec!["1".into()]);
        let encoded = EncodedBsxBatch::encode(batch.clone(), &chr_dtype).unwrap();
        let decoded = encoded.decode().unwrap();
        assert_eq!(batch, decoded);
    }
    
    #[test]
    fn test_trim() {
        let batch = BsxBatch::try_from(df![
                "chr" => ["1", "1", "1", "1", "1", "1", "1"],
                "position" => [2, 3, 4, 5, 7, 9, 10],
                "strand" => [".", ".", ".", ".", ".", ".", "."],
                "context" => ["CG", "CHG", "CHH", "CG", "CHG", "CHH", "CG"],
                "count_m" => [None::<u32>, None, None, None, None, None, None],
                "count_total" => [None::<u32>, None, None, None, None, None, None],
                "density" => [0, 1, 0, 0, 1, 0, 0]
            ].unwrap().lazy()
            .cast(BsxBatch::hashmap(), true)
            .collect().unwrap()
        ).unwrap();
        let chr_dtype = get_categorical_dtype(vec!["1".into()]);
        let encoded = EncodedBsxBatch::encode(batch.clone(), &chr_dtype).unwrap();
        
        let trimmed = encoded.trim_region::<u64>(
            &RegionCoordinates::new("1".to_string(), 4, 9)
        ).unwrap();
        assert_eq!(trimmed.first_position().unwrap().position(), 4);
        assert_eq!(trimmed.last_position().unwrap().position(), 9);

        let trimmed = encoded.trim_region::<u64>(
            &RegionCoordinates::new("1".to_string(), 1, 9)
        ).unwrap();
        assert_eq!(trimmed.first_position().unwrap().position(), 2);
        assert_eq!(trimmed.last_position().unwrap().position(), 9);

        let trimmed = encoded.trim_region::<u64>(
            &RegionCoordinates::new("1".to_string(), 4, 12)
        ).unwrap();
        assert_eq!(trimmed.first_position().unwrap().position(), 4);
        assert_eq!(trimmed.last_position().unwrap().position(), 10);
    }
}