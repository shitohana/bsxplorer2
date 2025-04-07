/*******************************************************************************
 Copyright (c) 2025
 The Prosperity Public License 3.0.0

 Contributor: [shitohana](https://github.com/shitohana)

 Source Code: https://github.com/shitohana/BSXplorer
 ******************************************************************************/

use crate::data_structs::batch::builder::BsxBatchBuilder;
use crate::data_structs::batch::decoded::BsxBatch;
use crate::data_structs::batch::traits::{BsxBatchMethods, BsxColNames};
use crate::data_structs::region::{GenomicPosition, RegionCoordinates};
use crate::tools::stats::MethylationStats;
use crate::utils::types::{Context, IPCEncodedEnum, Strand};
use anyhow::anyhow;
use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use bio_types::strand::NoStrand;
use itertools::{multiunzip, Itertools};
use once_cell::sync::OnceCell;
use polars::chunked_array::ChunkedArray;
use polars::datatypes::{BooleanChunked, BooleanType, CategoricalChunked, DataType, Float32Type, Int16Type, LogicalType, PlHashMap, PlSmallStr, UInt32Type};
use polars::error::{PolarsError, PolarsResult};
use polars::frame::DataFrame;
use polars::prelude::{col, lit, when, Expr, IntoLazy, NamedFrom, Schema};
use statrs::statistics::Statistics;
use std::cmp::Ordering;
use std::collections::HashMap;

/// Encoded version of [BsxBatch]
///
/// Encodes
/// 1. `context` as [bool]
/// 2. `strand` as [bool]
///
/// Other types are also downcasted
#[derive(Debug, Clone, PartialEq)]
pub struct EncodedBsxBatch {
    data: DataFrame,
    chr: OnceCell<String>,
    start: OnceCell<u32>,
    end: OnceCell<u32>,
}

impl BsxColNames for EncodedBsxBatch {}

impl EncodedBsxBatch {
    pub const POS_DTYPE: DataType = DataType::UInt32;
    pub const STRAND_DTYPE: DataType = DataType::Boolean;
    pub const CONTEXT_DTYPE: DataType = DataType::Boolean;
    pub const COUNT_DTYPE: DataType = DataType::Int16;
    pub const DENSITY_DTYPE: DataType = DataType::Float32;

    pub fn position(&self) -> &ChunkedArray<UInt32Type> {
        self.data
            .column(Self::POS_NAME)
            .unwrap()
            .u32()
            .unwrap()
    }

    pub fn chr(&self) -> &CategoricalChunked {
        self.data
            .column(Self::CHR_NAME)
            .unwrap()
            .categorical()
            .unwrap()
    }

    pub fn strand(&self) -> &ChunkedArray<BooleanType> {
        self.data
            .column(Self::STRAND_NAME)
            .unwrap()
            .bool()
            .unwrap()
    }

    pub fn context(&self) -> &ChunkedArray<BooleanType> {
        self.data
            .column(Self::CONTEXT_NAME)
            .unwrap()
            .bool()
            .unwrap()
    }

    pub fn count_m(&self) -> &ChunkedArray<Int16Type> {
        self.data
            .column(Self::COUNT_M_NAME)
            .unwrap()
            .i16()
            .unwrap()
    }

    pub fn count_total(&self) -> &ChunkedArray<Int16Type> {
        self.data
            .column(Self::COUNT_TOTAL_NAME)
            .unwrap()
            .i16()
            .unwrap()
    }

    pub fn density(&self) -> &ChunkedArray<Float32Type> {
        self.data
            .column(Self::DENSITY_NAME)
            .unwrap()
            .f32()
            .unwrap()
    }

    pub(crate) fn new_with_fields(
        data: DataFrame,
        chr: String,
        start: u32,
        end: u32,
    ) -> EncodedBsxBatch {
        let chr = OnceCell::with_value(chr);
        let start = OnceCell::with_value(start);
        let end = OnceCell::with_value(end);
        Self { data, chr, start, end }
    }
}

impl EncodedBsxBatch {
    /// Creates new [EncodedBsxBatch] without checks
    ///
    /// # Safety
    ///
    /// * input [DataFrame] is expected to follow format
    pub(crate) unsafe fn new_unchecked(data_frame: DataFrame) -> Self {
        EncodedBsxBatch { data: data_frame, chr: OnceCell::new(), start: OnceCell::new(), end: OnceCell::new() }
    }

    /// Trims [EncodedBsxBatch] by `region_coordinates`
    ///
    /// Returns [PolarsError] if
    /// 1. Chromosome names do not match
    /// 2. Batch does not contain region data_structs
    pub fn trim_region(
        &self,
        region_coordinates: &RegionCoordinates<u64>,
    ) -> anyhow::Result<Self>
    where
        Self: Sized, {
        let batch_first = GenomicPosition::new(
            self.chr_val()?.to_string(),
            self.start_pos().ok_or(anyhow!("Empty batch"))? as u64
        );;
        let pos_col = self.data().column("position")?.u32()?;
        let height = self.height();

        match batch_first.partial_cmp(&region_coordinates.end_gpos()) {
            None => {
                Err(PolarsError::ComputeError(
                    "Chromosome does not match".into(),
                )
                .into())
            },
            Some(Ordering::Greater) => {
                Err(PolarsError::ComputeError(
                    "Batch does not contain region information".into(),
                )
                .into())
            },
            _ => {
                let start_shift = pos_col
                    .iter()
                    .position(|v| {
                        v.unwrap() >= region_coordinates.start() as u32
                    })
                    .unwrap_or(height - 1);
                let end_shift = pos_col
                    .iter()
                    .skip(start_shift)
                    .position(|v| v.unwrap() > region_coordinates.end() as u32)
                    .map(|val| val + start_shift)
                    .unwrap_or(height);
                let mask = BooleanChunked::from_iter({
                    let mut start = vec![false; start_shift];
                    start.extend(vec![true; end_shift - start_shift]);
                    start.extend(vec![false; height - end_shift]);
                    start
                });
                Ok(self
                    .data()
                    .filter(&mask)
                    .map(|df| unsafe { Self::new_unchecked(df) })?)
            },
        }
    }

    // TODO: add checks
    /// Vertically stacks with other [EncodedBsxBatch]
    pub fn vstack(
        self,
        other: Self,
    ) -> PolarsResult<Self> {
        let new = self.data.vstack(&other.data)?;
        unsafe { Ok(Self::new_unchecked(new)) }
    }

    /// Get methylation stats for [EncodedBsxBatch]
    pub fn get_methylation_stats(&self) -> PolarsResult<MethylationStats> {
        // TODO move logic out
        let density_col = self.data().column("density")?;
        let mean_methylation = density_col
            .f32()?
            .iter()
            .map(|v| v.unwrap_or(0.0) as f64)
            .filter(|v| !v.is_nan())
            .mean();
        let variance_methylation = density_col
            .f32()?
            .iter()
            .map(|v| v.unwrap_or(0.0) as f64)
            .filter(|v| !v.is_nan())
            .variance();
        let coverage_distribution = self
            .data()
            .column("count_total")?
            .i16()?
            .iter()
            .fold(HashMap::<u16, u32>::new(), |mut counts, value| {
                *counts
                    .entry(value.unwrap() as u16)
                    .or_insert(0) += 1;
                counts
            });
        let context_methylation = {
            let context_col = self.data().column("context")?;

            let hashmap = itertools::izip!(
                context_col.bool()?.iter(),
                density_col.f32()?.iter()
            )
            .fold(
                HashMap::<Option<bool>, (f64, u32)>::new(),
                |mut stats, (context, density)| {
                    if density
                        .map(|v| !v.is_nan())
                        .unwrap_or(true)
                    {
                        let (sum, count) = stats
                            .entry(context)
                            .or_insert((0f64, 0));
                        *sum += density.unwrap_or(0f32) as f64;
                        *count += 1;
                    }
                    stats
                },
            );

            HashMap::<Context, (f64, u32)>::from_iter(
                hashmap
                    .into_iter()
                    .map(|(k, v)| (Context::from_bool(k), v)),
            )
        };
        let strand_methylation = {
            let strand_col = self.data().column("strand")?;

            let hashmap = itertools::izip!(
                strand_col.bool()?.iter(),
                density_col.f32()?.iter()
            )
            .fold(
                HashMap::<Option<bool>, (f64, u32)>::new(),
                |mut stats, (strand, density)| {
                    if density
                        .map(|v| !v.is_nan())
                        .unwrap_or(true)
                    {
                        let (sum, count) =
                            stats.entry(strand).or_insert((0f64, 0));
                        *sum += density.unwrap_or(0f32) as f64;
                        *count += 1;
                    }
                    stats
                },
            );

            HashMap::<Strand, (f64, u32)>::from_iter(
                hashmap
                    .into_iter()
                    .map(|(k, v)| (Strand::from_bool(k), v)),
            )
        };

        Ok(MethylationStats::from_data(
            mean_methylation,
            variance_methylation,
            coverage_distribution,
            context_methylation,
            strand_methylation,
        ))
    }

    /// Creates [Contig] from [EncodedBsxBatch]
    pub fn as_contig(&self) -> anyhow::Result<Contig<&str, NoStrand>> {
        let start = self.start_pos().ok_or(anyhow!("no data"))?;
        let end = self.end_pos().ok_or(anyhow!("no data"))?;
        let chr = self.chr_val()?;

        Ok(Contig::new(
            chr,
            start as isize,
            (end + 1 - start) as usize,
            NoStrand::Unknown,
        ))
    }

    /// Replaces low coverage methylation sites to count_total = 0, density =
    /// NaN
    pub fn mark_low_counts(
        self,
        threshold: i16,
    ) -> PolarsResult<Self> {
        let data = self.data.lazy().with_column(
            when(col("count_total").lt(lit(threshold as u32))).then(lit(0)).otherwise(col("count_total")).alias("count_total"),
        ).with_columns([
            when(col("count_total").eq(lit(0))).then(lit(0)).otherwise(col("count_m")).alias("count_m"),
            when(col("count_total").eq(lit(0))).then(lit(f64::NAN)).otherwise(col("density")).alias("density"),
        ]).collect()?;
        Ok(Self { data: data, chr: self.chr, start: self.start, end: self.end })
    }

    /// Returns position values
    pub fn get_position_vals(&self) -> PolarsResult<Vec<u32>> {
        Ok(self.position()
            .into_iter()
            .map(|x| x.unwrap())
            .collect())
    }

    /// Filter [EncodedBsxBatch] by mask
    pub fn filter_mask(
        self,
        mask: &BooleanChunked,
    ) -> PolarsResult<Self> {
        Ok(unsafe { Self::new_unchecked(self.data.filter(mask)?) })
    }

    /// Split [EncodedBsxBatch] into three batches by context (CG, CHG, CHH)
    pub fn strip_contexts(self) -> PolarsResult<(Self, Self, Self)> {
        let context_col = self.data().column("context")?.bool()?;
        let cg_encoding = Context::CG.to_bool().unwrap();
        let (cg_mask, chg_mask, chh_mask): (Vec<bool>, Vec<bool>, Vec<bool>) =
            multiunzip(context_col.iter().map(|v| {
                let (cg, chh) = if let Some(boolean) = v {
                    (boolean && cg_encoding, false)
                }
                else {
                    (false, true)
                };
                (cg, cg == chh, chh)
            }));
        Ok((
            unsafe {
                Self::new_unchecked(
                    self.data
                        .filter(&BooleanChunked::new("mask".into(), cg_mask))?,
                )
            },
            unsafe {
                Self::new_unchecked(
                    self.data.filter(&BooleanChunked::new(
                        "mask".into(),
                        chg_mask,
                    ))?,
                )
            },
            unsafe {
                Self::new_unchecked(
                    self.data.filter(&BooleanChunked::new(
                        "mask".into(),
                        chh_mask,
                    ))?,
                )
            },
        ))
    }
}

impl From<EncodedBsxBatch> for DataFrame {
    fn from(batch: EncodedBsxBatch) -> Self { batch.data }
}

impl BsxBatchMethods for EncodedBsxBatch {
    fn filter(
        self,
        expr: Expr,
    ) -> anyhow::Result<Self>
      where
          Self: Sized,
    {
        let dtype = self.chr().dtype().clone();
        let df = self.data.lazy().filter(expr).collect()?;
        BsxBatchBuilder::no_checks().build_encoded(df, dtype)
    }
    /// Returns reference to inner [DataFrame]
    fn data(&self) -> &DataFrame { &self.data }
    /// Returns mutable reference to inner [DataFrame]
    fn data_mut(&mut self) -> &mut DataFrame { &mut self.data }

    /// Returns chromosome name from chromosome column
    fn chr_val(&self) -> anyhow::Result<&str> {
        let val = self.chr.get_or_try_init(|| {
            let rev_map = self.chr().get_rev_map();
            let first_value = self.chr().physical().first();
            first_value.map(|v| rev_map.get(v).to_string()).ok_or(anyhow::anyhow!("No data!"))
        })?.as_str();
        Ok(val)
    }

    fn start_pos(&self) -> Option<u32> {
        self.start.get_or_try_init(
            || self.position().first().map(|v| v.into()).ok_or(anyhow::anyhow!("empty")),
        ).ok().cloned()
    }

    fn end_pos(&self) -> Option<u32> {
        self.end.get_or_try_init(
            || self.position().last().map(|v| v.into()).ok_or(anyhow::anyhow!("empty")),
        ).ok().cloned()
    }
}

pub(crate) fn get_encoded_schema(chr_dtype: &DataType) -> Schema {
    Schema::from_iter(get_encoded_hashmap(chr_dtype).into_iter().map(|(k, v)| (PlSmallStr::from(k), v)))
}

pub(crate) fn get_encoded_hashmap(
    chr_dtype: &DataType
) -> PlHashMap<&str, DataType> {
    PlHashMap::from_iter([
        (EncodedBsxBatch::CHR_NAME, chr_dtype.clone()),
        (EncodedBsxBatch::POS_NAME, EncodedBsxBatch::POS_DTYPE.clone()),
        (EncodedBsxBatch::STRAND_NAME, EncodedBsxBatch::STRAND_DTYPE.clone()),
        (EncodedBsxBatch::CONTEXT_NAME, EncodedBsxBatch::CONTEXT_DTYPE.clone()),
        (EncodedBsxBatch::COUNT_M_NAME, EncodedBsxBatch::COUNT_DTYPE.clone()),
        (EncodedBsxBatch::COUNT_TOTAL_NAME, EncodedBsxBatch::COUNT_DTYPE.clone()),
        (EncodedBsxBatch::DENSITY_NAME, EncodedBsxBatch::DENSITY_DTYPE.clone()),
    ])
}