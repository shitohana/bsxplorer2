/*******************************************************************************
Copyright (c) 2025
The Prosperity Public License 3.0.0

Contributor: [shitohana](https://github.com/shitohana)

Source Code: https://github.com/shitohana/BSXplorer
******************************************************************************/

use super::traits::BsxTypeTag;
use crate::data_structs::batch::traits::{BsxBatchMethods, BsxColNames};
use crate::data_structs::region::{GenomicPosition, RegionCoordinates};
use crate::utils::types::{Context, IPCEncodedEnum};
use anyhow::anyhow;
use itertools::{multiunzip, Itertools};
use once_cell::sync::OnceCell;
use polars::prelude::*;
use std::cmp::Ordering;

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
    pub fn chr(&self) -> &CategoricalChunked {
        self.data
            .column(Self::CHR_NAME)
            .unwrap()
            .categorical()
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
        Self {
            data,
            chr,
            start,
            end,
        }
    }
}

impl EncodedBsxBatch {
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
        Self: Sized,
    {
        let batch_first = GenomicPosition::new(
            self.chr_val()?.to_string(),
            self.start_pos()
                .ok_or(anyhow!("Empty batch"))? as u64,
        );
        let pos_col = self.data().column("position")?.u32()?;
        let height = self.height();

        match batch_first.partial_cmp(&region_coordinates.end_gpos()) {
            None => Err(PolarsError::ComputeError(
                "Chromosome does not match".into(),
            )
            .into()),
            Some(Ordering::Greater) => Err(PolarsError::ComputeError(
                "Batch does not contain region information".into(),
            )
            .into()),
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

    /// Replaces low coverage methylation sites to count_total = 0, density =
    /// NaN
    pub fn mark_low_counts(
        self,
        threshold: i16,
    ) -> PolarsResult<Self> {
        let data = self
            .data
            .lazy()
            .with_column(
                when(col("count_total").lt(lit(threshold as u32)))
                    .then(lit(0))
                    .otherwise(col("count_total"))
                    .alias("count_total"),
            )
            .with_columns([
                when(col("count_total").eq(lit(0)))
                    .then(lit(0))
                    .otherwise(col("count_m"))
                    .alias("count_m"),
                when(col("count_total").eq(lit(0)))
                    .then(lit(f64::NAN))
                    .otherwise(col("density"))
                    .alias("density"),
            ])
            .collect()?;
        Ok(Self {
            data,
            chr: self.chr,
            start: self.start,
            end: self.end,
        })
    }

    /// Split [EncodedBsxBatch] into three batches by context (CG, CHG, CHH)
    pub fn strip_contexts(self) -> PolarsResult<(Self, Self, Self)> {
        let context_col = self.data().column("context")?.bool()?;
        let cg_encoding = Context::CG.to_bool().unwrap();
        let (cg_mask, chg_mask, chh_mask): (Vec<bool>, Vec<bool>, Vec<bool>) =
            multiunzip(context_col.iter().map(|v| {
                let (cg, chh) = if let Some(boolean) = v {
                    (boolean && cg_encoding, false)
                } else {
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
    fn from(batch: EncodedBsxBatch) -> Self {
        batch.data
    }
}

impl BsxBatchMethods for EncodedBsxBatch {
    type ChrType = UInt32Type;
    type PosType = UInt32Type;
    type StrandType = BooleanType;
    type ContextType = BooleanType;
    type CountType = Int16Type;
    type DensityType = Float32Type;

    fn chr(&self) -> &ChunkedArray<Self::ChrType> {
        self.data
            .column(Self::CHR_NAME)
            .unwrap()
            .categorical()
            .unwrap()
            .physical()
    }
    fn position(&self) -> &ChunkedArray<Self::PosType> {
        self.data
            .column(Self::POS_NAME)
            .unwrap()
            .u32()
            .unwrap()
    }
    fn strand(&self) -> &ChunkedArray<Self::StrandType> {
        self.data
            .column(Self::STRAND_NAME)
            .unwrap()
            .bool()
            .unwrap()
    }

    fn context(&self) -> &ChunkedArray<Self::ContextType> {
        self.data
            .column(Self::CONTEXT_NAME)
            .unwrap()
            .bool()
            .unwrap()
    }

    fn count_m(&self) -> &ChunkedArray<Self::CountType> {
        self.data
            .column(Self::COUNT_M_NAME)
            .unwrap()
            .i16()
            .unwrap()
    }

    fn count_total(&self) -> &ChunkedArray<Self::CountType> {
        self.data
            .column(Self::COUNT_TOTAL_NAME)
            .unwrap()
            .i16()
            .unwrap()
    }

    fn density(&self) -> &ChunkedArray<Self::DensityType> {
        self.data
            .column(Self::DENSITY_NAME)
            .unwrap()
            .f32()
            .unwrap()
    }

    unsafe fn new_unchecked(data_frame: DataFrame) -> Self {
        EncodedBsxBatch {
            data: data_frame,
            chr: OnceCell::new(),
            start: OnceCell::new(),
            end: OnceCell::new(),
        }
    }

    /// Returns reference to inner [DataFrame]
    fn data(&self) -> &DataFrame {
        &self.data
    }
    /// Returns mutable reference to inner [DataFrame]
    fn data_mut(&mut self) -> &mut DataFrame {
        &mut self.data
    }

    /// Returns chromosome name from chromosome column
    fn chr_val(&self) -> anyhow::Result<&str> {
        let val = self
            .chr
            .get_or_try_init(|| {
                let rev_map = self.chr().get_rev_map();
                let first_value = self.chr().physical().first();
                first_value
                    .map(|v| rev_map.get(v).to_string())
                    .ok_or(anyhow::anyhow!("No data!"))
            })?
            .as_str();
        Ok(val)
    }

    fn start_pos(&self) -> Option<u32> {
        self.start
            .get_or_try_init(|| {
                self.position()
                    .first()
                    .map(|v| v.into())
                    .ok_or(anyhow::anyhow!("empty"))
            })
            .ok()
            .cloned()
    }

    fn end_pos(&self) -> Option<u32> {
        self.end
            .get_or_try_init(|| {
                self.position()
                    .last()
                    .map(|v| v.into())
                    .ok_or(anyhow::anyhow!("empty"))
            })
            .ok()
            .cloned()
    }
}

impl BsxTypeTag for EncodedBsxBatch {
    fn type_name() -> &'static str {
        "encoded"
    }
}
