/*******************************************************************************
 Copyright (c) 2025
 The Prosperity Public License 3.0.0

 Contributor: [shitohana](https://github.com/shitohana)

 Source Code: https://github.com/shitohana/BSXplorer
 ******************************************************************************/

/// ****************************************************************************
/// * Copyright (c) 2025
/// ***************************************************************************

//! Module contains BsxBatches which represent single chromosome methylation
//! data_structs
//!
//! BsxBatch stores data_structs in [DataFrame] with
//! * Non-null chr and position
//! * Single unique chromosome
//! * Ascending sorted positions
//!
//! Module contains two main structs:
//! 1. [BsxBatch] struct for storing non-encoded data_structs
//! 2. [EncodedBsxBatch] struct for storing encoded data_structs (context and
//! strand encoded as boolean)
//!
//! and [BsxBatchMethods] trait
//!
//! Both BsxBatch and EncodedBsxBatch can be filtered using context and/or
//! strand
use std::cmp::Ordering;
use std::collections::HashMap;
use std::ops::BitAnd;

use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use bio_types::annot::pos::Pos;
use bio_types::annot::refids::RefIDSet;
use bio_types::strand::NoStrand;
use itertools::{multiunzip, Itertools};
use log::warn;
use polars::prelude::*;
use statrs::statistics::Statistics;

use crate::data_structs::region::{GenomicPosition, RegionCoordinates};
use crate::tools::stats::MethylationStats;
use crate::utils::types::{Context, IPCEncodedEnum, Strand};
// noinspection ALL
use crate::utils::{decode_context,
                   decode_strand,
                   encode_context,
                   encode_strand,
                   polars_schema};

/// DataFrame with
/// 1. Non-null chr and position
/// 2. Single unique chr
/// 3. Ascending sorted positions
#[derive(Clone, PartialEq, Debug)]
#[cfg_attr(feature = "python", pyclass)]
pub struct BsxBatch(DataFrame);

impl BsxBatch {
    /// Creates new [BsxBatch] without checks
    ///
    /// # Safety
    ///
    /// * input [DataFrame] is expected to follow format
    pub(crate) unsafe fn new_unchecked(data_frame: DataFrame) -> Self {
        BsxBatch(data_frame)
    }

    /// Returns expected column names of dataframe
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

    /// Returns expected types of columns
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

    /// Returns name of chr column
    #[inline(always)]
    pub(crate) const fn chr_col() -> &'static str { "chr" }

    /// Returns name of position column
    #[inline(always)]
    pub(crate) const fn pos_col() -> &'static str { "position" }

    /// Returns expected schema of [BsxBatch]
    pub fn schema() -> Schema {
        use crate::utils::schema_from_arrays;
        schema_from_arrays(Self::col_names(), Self::col_types())
    }

    /// Returns expected schema of [BsxBatch] as [PlHashMap]
    pub fn hashmap() -> PlHashMap<&'static str, DataType> {
        use crate::utils::hashmap_from_arrays;
        hashmap_from_arrays(Self::col_names(), Self::col_types())
    }
}

impl BsxBatchMethods for BsxBatch {
    /// Filters [BsxBatch] by `context` and `strand`
    ///
    /// If both `context` and `strand` are [None], returns [BsxBatch]
    fn filter(
        self,
        context: Option<Context>,
        strand: Option<Strand>,
    ) -> Self where
        Self: Sized,
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

    /// Returns reference to inner [DataFrame]
    #[inline]
    fn data(&self) -> &DataFrame { &self.0 }

    /// Returns mutable reference to inner [DataFrame]
    #[inline]
    fn data_mut(&mut self) -> &mut DataFrame { &mut self.0 }
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

        for colname in [
            Self::chr_col(),
            Self::pos_col(),
            "context",
            "strand",
            "count_m",
            "count_total",
        ] {
            if data_casted.column(colname)?.null_count() != 0
            {
                return Err(PolarsError::ComputeError(
                    format!("Nulls in {} col", colname).into(),
                ));
            }
        }

        // Check sorted
        if !data_casted.column(Self::pos_col())?.as_series().unwrap().is_sorted(SortOptions::default().with_order_descending(false))?
        {
            return Err(PolarsError::ComputeError("Data is not sorted".into()));
        }

        Ok(BsxBatch(data_casted))
    }
}

impl From<BsxBatch> for DataFrame {
    fn from(b: BsxBatch) -> Self { b.0 }
}

/// Encoded version of [BsxBatch]
///
/// Encodes
/// 1. `context` as [bool]
/// 2. `strand` as [bool]
///
/// Other types are also downcasted
#[derive(Debug, Clone, PartialEq)]
pub struct EncodedBsxBatch(DataFrame);

impl EncodedBsxBatch {
    /// Creates new [EncodedBsxBatch] without checks
    ///
    /// # Safety
    ///
    /// * input [DataFrame] is expected to follow format
    pub(crate) unsafe fn new_unchecked(data_frame: DataFrame) -> Self {
        EncodedBsxBatch(data_frame)
    }

    /// Encodes [BsxBatch] into [EncodedBsxBatch]
    ///
    /// Performs
    /// 1. Encoding of `context` and `strand`
    /// 2. Downcasting of column types
    pub fn encode(
        batch: BsxBatch,
        chr_dtype: &DataType,
    ) -> PolarsResult<Self> {
        let batch_data = DataFrame::from(batch);
        let target_schema = Self::get_hashmap(chr_dtype);

        let mut batch_lazy = batch_data.lazy();

        batch_lazy = encode_context(batch_lazy, "context");
        batch_lazy = encode_strand(batch_lazy, "strand");

        batch_lazy = batch_lazy.cast(target_schema, true).select(
            Self::col_names().iter().cloned().map(col).collect_vec(),
        );
        let mut result = batch_lazy.collect()?;
        result.rechunk_mut();

        Ok(Self(result))
    }

    /// Decodes [EncodedBsxBatch] into [BsxBatch]
    ///
    /// Performs
    /// 1. Decoding of `context` and `strand`
    /// 2. Upcasting of column types
    pub fn decode(self) -> PolarsResult<BsxBatch> {
        let batch_data = DataFrame::from(self);
        let target_schema = BsxBatch::hashmap();
        let mut batch_lazy = batch_data.lazy();
        batch_lazy = decode_context(batch_lazy, "context", "context");
        batch_lazy = decode_strand(batch_lazy, "strand", "strand");
        batch_lazy = batch_lazy.cast(target_schema, true).select(
            BsxBatch::col_names().iter().cloned().map(col).collect_vec(),
        );

        let result = batch_lazy.collect()?;
        Ok(BsxBatch(result))
    }

    /// Returns expected schema of [EncodedBsxBatch]
    pub(crate) fn get_schema(chr_dtype: &DataType) -> Schema {
        polars_schema![
            "chr" => chr_dtype.clone(),
            "position" => DataType::UInt32,
            "strand" => DataType::Boolean,
            "context" => DataType::Boolean,
            "count_m" => DataType::Int16,
            "count_total" => DataType::Int16,
            "density" => DataType::Float32
        ]
    }

    /// Returns expected schema of [EncodedBsxBatch] as [PlHashMap]
    pub(crate) fn get_hashmap(
        chr_dtype: &DataType
    ) -> PlHashMap<&str, DataType> {
        PlHashMap::from_iter([
            ("chr", chr_dtype.clone()),
            ("position", DataType::UInt32),
            ("strand", DataType::Boolean),
            ("context", DataType::Boolean),
            ("count_m", DataType::Int16),
            ("count_total", DataType::Int16),
            ("density", DataType::Float32),
        ])
    }

    /// Trims [EncodedBsxBatch] by `region_coordinates`
    ///
    /// Returns [PolarsError] if
    /// 1. Chromosome names do not match
    /// 2. Batch does not contain region data_structs
    pub fn trim_region(
        &self,
        region_coordinates: &RegionCoordinates<u64>,
    ) -> anyhow::Result<Self> where
        Self: Sized,
    {
        let batch_first = self.first_position()?;
        let pos_col = self.data().column("position")?.u32()?;
        let height = self.height();

        match batch_first.partial_cmp(&region_coordinates.end_gpos()) {
            None => {
                Err(PolarsError::ComputeError(
                    "Chromosome does not match".into(),
                ).into())
            }
            Some(Ordering::Greater) => {
                Err(PolarsError::ComputeError(
                    "Batch does not contain region information".into(),
                ).into())
            }
            _ => {
                let start_shift = pos_col.iter().position(|v| {
                    v.unwrap() >= region_coordinates.start() as u32
                }).unwrap_or(height - 1);
                let end_shift = pos_col.iter().skip(start_shift).position(|v| v.unwrap() > region_coordinates.end() as u32).map(|val| val + start_shift).unwrap_or(height);
                let mask = BooleanChunked::from_iter({
                    let mut start = vec![false; start_shift];
                    start.extend(vec![true; end_shift - start_shift]);
                    start.extend(vec![false; height - end_shift]);
                    start
                });
                Ok(self.data().filter(&mask).map(|df| unsafe { Self::new_unchecked(df) })?)
            }
        }
    }

    /// Returns expected schema of inner [DataFrame]
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn schema(&self) -> Schema { self.data().schema() }

    /// Returns expected column names of inner [DataFrame]
    #[inline]
    pub(crate) fn col_names() -> &'static [&'static str] {
        BsxBatch::col_names()
    }

    // TODO: add checks
    /// Vertically stacks with other [EncodedBsxBatch]
    pub fn vstack(
        self,
        other: Self,
    ) -> PolarsResult<Self> {
        let new = self.0.vstack(&other.0)?;
        unsafe { Ok(Self::new_unchecked(new)) }
    }

    /// Get methylation stats for [EncodedBsxBatch]
    pub fn get_methylation_stats(&self) -> PolarsResult<MethylationStats> {
        let density_col = self.data().column("density")?;
        let mean_methylation = density_col.f32()?.iter().map(|v| v.unwrap_or(0.0) as f64).filter(|v| !v.is_nan()).mean();
        let variance_methylation = density_col.f32()?.iter().map(|v| v.unwrap_or(0.0) as f64).filter(|v| !v.is_nan()).variance();
        let coverage_distribution = self.data().column("count_total")?.i16()?.iter().fold(HashMap::<u16, u32>::new(), |mut counts, value| {
            *counts.entry(value.unwrap() as u16).or_insert(0) += 1;
            counts
        });
        let context_methylation = {
            let context_col = self.data().column("context")?;

            let hashmap = itertools::izip!(
                context_col.bool()?.iter(),
                density_col.f32()?.iter()
            ).fold(
                HashMap::<Option<bool>, (f64, u32)>::new(),
                |mut stats, (context, density)| {
                    if density.map(|v| !v.is_nan()).unwrap_or(true)
                    {
                        let (sum, count) = stats.entry(context).or_insert((0f64, 0));
                        *sum += density.unwrap_or(0f32) as f64;
                        *count += 1;
                    }
                    stats
                },
            );

            HashMap::<Context, (f64, u32)>::from_iter(
                hashmap.into_iter().map(|(k, v)| (Context::from_bool(k), v)),
            )
        };
        let strand_methylation = {
            let strand_col = self.data().column("strand")?;

            let hashmap = itertools::izip!(
                strand_col.bool()?.iter(),
                density_col.f32()?.iter()
            ).fold(
                HashMap::<Option<bool>, (f64, u32)>::new(),
                |mut stats, (strand, density)| {
                    if density.map(|v| !v.is_nan()).unwrap_or(true)
                    {
                        let (sum, count) = stats.entry(strand).or_insert((0f64, 0));
                        *sum += density.unwrap_or(0f32) as f64;
                        *count += 1;
                    }
                    stats
                },
            );

            HashMap::<Strand, (f64, u32)>::from_iter(
                hashmap.into_iter().map(|(k, v)| (Strand::from_bool(k), v)),
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

    /// Returns chromosome name from chromosome column
    pub fn chr(&self) -> PolarsResult<String> {
        let chr_col = self.data().column("chr")?.categorical()?;
        let chr = chr_col.get_rev_map().get(chr_col.physical().first().ok_or(
            PolarsError::ComputeError("Could not get first id".into()),
        )?);
        Ok(chr.to_string())
    }

    /// Returns first available position
    pub fn first_seq_pos(
        &self,
        refid_set: &mut Option<RefIDSet<Arc<String>>>,
    ) -> PolarsResult<Pos<Arc<String>, NoStrand>> {
        let pos = self.data().column("position")?.u32()?.first().unwrap();
        let chr = if let Some(refids) = refid_set {
            refids.intern(self.chr()?.as_str())
        } else {
            Arc::new(self.chr()?.to_owned())
        };
        Ok(Pos::new(chr, pos as isize, NoStrand::Unknown))
    }

    /// Returns last available position
    pub fn last_seq_pos(
        &self,
        refid_set: &mut Option<RefIDSet<Arc<String>>>,
    ) -> PolarsResult<Pos<Arc<String>, NoStrand>> {
        let pos = self.data().column("position")?.u32()?.last().unwrap();
        let chr = if let Some(refids) = refid_set {
            refids.intern(self.chr()?.as_str())
        } else {
            Arc::new(self.chr()?.to_owned())
        };
        Ok(Pos::new(chr, pos as isize, NoStrand::Unknown))
    }

    /// Creates [Contig] from [EncodedBsxBatch]
    pub fn as_contig(
        &self,
        refid_set: &mut Option<RefIDSet<Arc<String>>>,
    ) -> PolarsResult<Contig<Arc<String>, NoStrand>> {
        let start = self.first_seq_pos(refid_set)?;
        let end = self.last_seq_pos(refid_set)?;
        if start.refid() != end.refid() {
            Err(PolarsError::ComputeError(
                "Different start and end refererence IDs".into(),
            ))
        } else {
            Ok(Contig::new(
                start.refid().clone(),
                start.pos(),
                (end.pos() + 1 - start.pos()) as usize,
                NoStrand::Unknown,
            ))
        }
    }

    /// Replaces low coverage methylation sites to count_total = 0, density =
    /// NaN
    pub fn mark_low_counts(
        self,
        threshold: i16,
    ) -> PolarsResult<Self> {
        let data = self.0.lazy().with_column(
            when(col("count_total").lt(lit(threshold as u32))).then(lit(0)).otherwise(col("count_total")).alias("count_total"),
        ).with_columns([
            when(col("count_total").eq(lit(0))).then(lit(0)).otherwise(col("count_m")).alias("count_m"),
            when(col("count_total").eq(lit(0))).then(lit(f64::NAN)).otherwise(col("density")).alias("density"),
        ]).collect()?;
        Ok(Self(data))
    }

    /// Returns methylation density values
    pub fn get_density_vals(&self) -> PolarsResult<Vec<f32>> {
        Ok(self.0.column("density")?.f32()?.into_iter().map(|x| x.unwrap_or(f32::NAN)).collect())
    }

    /// Returns position values
    pub fn get_position_vals(&self) -> PolarsResult<Vec<u32>> {
        Ok(self.0.column("position")?.u32()?.into_iter().map(|x| x.unwrap()).collect())
    }

    /// Returns methylated counts
    pub fn get_counts_m(&self) -> PolarsResult<Vec<i16>> {
        Ok(self.0.column("count_m")?.i16()?.into_iter().map(|x| x.unwrap_or(0)).collect())
    }

    /// Returns total counts
    pub fn get_counts_total(&self) -> PolarsResult<Vec<i16>> {
        Ok(self.0.column("count_total")?.i16()?.into_iter().map(|x| x.unwrap_or(0)).collect())
    }

    /// Filter [EncodedBsxBatch] by mask
    pub fn filter_mask(
        self,
        mask: &BooleanChunked,
    ) -> PolarsResult<Self> {
        Ok(unsafe { Self::new_unchecked(self.0.filter(mask)?) })
    }

    /// Split [EncodedBsxBatch] into three batches by context (CG, CHG, CHH)
    pub fn strip_contexts(self) -> PolarsResult<(Self, Self, Self)> {
        let context_col = self.data().column("context")?.bool()?;
        let cg_encoding = Context::CG.to_bool().unwrap();
        let (cg_mask, chg_mask, chh_mask): (Vec<bool>, Vec<bool>, Vec<bool>) = multiunzip(context_col.iter().map(|v| {
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
                    self.0.filter(&BooleanChunked::new("mask".into(), cg_mask))?,
                )
            },
            unsafe {
                Self::new_unchecked(
                    self.0.filter(&BooleanChunked::new(
                        "mask".into(),
                        chg_mask,
                    ))?,
                )
            },
            unsafe {
                Self::new_unchecked(
                    self.0.filter(&BooleanChunked::new(
                        "mask".into(),
                        chh_mask,
                    ))?,
                )
            },
        ))
    }
}

impl From<EncodedBsxBatch> for DataFrame {
    fn from(batch: EncodedBsxBatch) -> Self { batch.0 }
}

impl BsxBatchMethods for EncodedBsxBatch {
    /// Filters [EncodedBsxBatch] by `context` and `strand`
    ///
    /// If both `context` and `strand` are [None], returns [EncodedBsxBatch]
    fn filter(
        self,
        context: Option<Context>,
        strand: Option<Strand>,
    ) -> Self where
        Self: Sized,
    {
        if context.is_none() && strand.is_none() {
            self
        } else {
            let mut boolean_chunked = BooleanChunked::new("mask".into(), vec![true; self.0.height()]);
            if let Some(strand) = strand {
                let mask = BooleanChunked::new(
                    "strand".into(),
                    self.0.column("strand").unwrap().bool().unwrap().into_iter().map(|v| v == strand.to_bool()).collect_vec(),
                );
                boolean_chunked = boolean_chunked.bitand(mask);
            }
            if let Some(context) = context {
                let mask = BooleanChunked::new(
                    "context".into(),
                    self.0.column("context").unwrap().bool().unwrap().into_iter().map(|v| v == context.to_bool()).collect_vec(),
                );
                boolean_chunked = boolean_chunked.bitand(mask);
            }
            Self(self.0.filter(&boolean_chunked).unwrap())
        }
    }

    /// Returns reference to inner [DataFrame]
    fn data(&self) -> &DataFrame { &self.0 }

    /// Returns mutable reference to inner [DataFrame]
    fn data_mut(&mut self) -> &mut DataFrame { &mut self.0 }

    /// Returns first [GenomicPosition]
    fn first_position(&self) -> anyhow::Result<GenomicPosition<u64>> {
        let chr_col = self.data().column("chr")?.categorical()?;
        let chr = chr_col.get_rev_map().get(chr_col.physical().first().unwrap());
        let pos = self.data().column("position")?.u32()?.first().unwrap();
        Ok(GenomicPosition::new(chr.to_string(), pos as u64))
    }

    /// Returns last [GenomicPosition]
    fn last_position(&self) -> anyhow::Result<GenomicPosition<u64>> {
        let chr_col = self.data().column("chr")?.categorical()?;
        let chr = chr_col.get_rev_map().get(chr_col.physical().last().unwrap());
        let pos = self.data().column("position")?.u32()?.last().unwrap();
        Ok(GenomicPosition::new(chr.to_string(), pos as u64))
    }
}

/// Trait for common methods for [BsxBatch] and [EncodedBsxBatch]
pub trait BsxBatchMethods {
    /// Filters [BsxBatch] by `context` and `strand`
    ///
    /// If both `context` and `strand` are [None], returns [BsxBatch]
    fn filter(
        self,
        context: Option<Context>,
        strand: Option<Strand>,
    ) -> Self where
        Self: Sized;

    /// Returns reference to inner [DataFrame]
    fn data(&self) -> &DataFrame;

    /// Returns mutable reference to inner [DataFrame]
    fn data_mut(&mut self) -> &mut DataFrame;

    /// Checks if chromosome column contains only one chromosome
    fn check_chr_unique(&self) -> bool {
        self.data().column("chr").unwrap().unique().unwrap().len() == 1
    }

    /// Returns first [GenomicPosition]
    fn first_position(&self) -> anyhow::Result<GenomicPosition<u64>> {
        use crate::utils::first_position;
        first_position(self.data(), BsxBatch::chr_col(), BsxBatch::pos_col())
    }
    /// Returns last [GenomicPosition]
    fn last_position(&self) -> anyhow::Result<GenomicPosition<u64>> {
        use crate::utils::last_position;
        last_position(self.data(), BsxBatch::chr_col(), BsxBatch::pos_col())
    }

    /// Extends [BsxBatch] by other
    ///
    /// Prints warning if first position of other is not less than last of self,
    /// but still sorts data_structs by position.
    ///
    /// Returns error if
    /// 1. Chromosome names do not match
    /// 2. Chromosome columns non-unique
    /// 3. Resulting data_structs still contains duplicates
    fn extend(
        &mut self,
        other: &Self,
    ) -> anyhow::Result<()> where
        Self: Sized,
    {
        if !(self.check_chr_unique() && other.check_chr_unique()) {
            return Err(PolarsError::ComputeError(
                "Chromosomes in batches non-unique".into(),
            ).into());
        }
        let self_pos = self.last_position()?;
        let other_pos = self.first_position()?;

        if self_pos.chr() != other_pos.chr() {
            return Err(PolarsError::ComputeError(
                "Chromosomes in batches differ".into(),
            ).into());
        }
        if self_pos.position() >= other_pos.position() {
            warn!(
                "First position in other batch ({}) must be less than last \
                 position in the current ({})",
                other_pos, self_pos
            );
            self.data_mut().extend(other.data())?;
            self.data_mut().sort_in_place(
                ["position"],
                SortMultipleOptions::default().with_order_descending(false).with_multithreaded(true),
            )?;

            if self.data().column("position")?.n_unique()? != self.data().height()
            {
                println!(
                    "{} != {}",
                    self.data().column("position")?.n_unique()?,
                    self.data().height()
                );
                println!("{}", self.data());
                return Err(PolarsError::ComputeError(
                    "Position values are duplicated".into(),
                ).into());
            }
        }

        Ok(self.data_mut().extend(other.data())?)
    }

    /// Returns number of rows
    fn height(&self) -> usize { self.data().height() }
}

#[cfg(feature = "python")]
#[pymethods]
impl EncodedBsxBatch {
    #[new]
    pub fn new_py(
        batch: BsxBatch,
        chr_names: Vec<String>,
    ) -> PyResult<EncodedBsxBatch> {
        let chr_dtype = get_categorical_dtype(chr_names);
        wrap_polars_result!(Self::encode(batch, &chr_dtype))
    }

    #[pyo3(name = "decode")]
    pub fn decode_py(&self) -> PyResult<BsxBatch> {
        let clone = self.clone();
        wrap_polars_result!(clone.decode())
    }

    #[pyo3(name = "trim_region")]
    pub fn trim_region_py(
        &self,
        region: PyRegionCoordinates,
    ) -> PyResult<Self> {
        let region_coordinates: RegionCoordinates<_> = region.into();
        wrap_polars_result!(self.trim_region(&region_coordinates))
    }

    #[getter]
    #[pyo3(name = "schema")]
    pub fn schema_py(&self) -> PyResult<PySchema> {
        Ok(PySchema(SchemaRef::new(self.schema())))
    }

    #[pyo3(name = "vstack")]
    pub fn vstack_py(
        &self,
        other: &Self,
    ) -> PyResult<Self> {
        wrap_polars_result!(self.clone().vstack(other.clone()))
    }

    #[getter]
    fn get_data(&self) -> PyDataFrame { PyDataFrame(self.data().clone()) }

    #[pyo3(name="filter", signature = (context=None, strand=None))]
    fn filter_py(
        &self,
        context: Option<Context>,
        strand: Option<Strand>,
    ) -> Self {
        self.clone().filter(context, strand)
    }

    #[getter]
    fn get_first_position(&self) -> PyResult<PyGenomicPosition> {
        wrap_polars_result!(self
            .first_position()
            .map(PyGenomicPosition::from))
    }

    #[getter]
    fn get_last_position(&self) -> PyResult<PyGenomicPosition> {
        wrap_polars_result!(self
            .first_position()
            .map(PyGenomicPosition::from))
    }

    #[getter]
    fn get_height(&self) -> usize { self.height() }

    fn __len__(&self) -> usize { self.height() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::get_categorical_dtype;

    fn dummy_batch() -> BsxBatch {
        BsxBatch::try_from(
            df![
                "chr" => ["1", "1", "1"],
                "position" => [1, 2, 3],
                "strand" => [".", ".", "."],
                "context" => ["CG", "CHG", "CHH"],
                "count_m" => [0,0,0],
                "count_total" => [0,0,0],
                "density" => [0, 1, 0]
            ].unwrap().lazy().cast(BsxBatch::hashmap(), true).collect().unwrap(),
        ).unwrap()
    }

    #[test]
    fn test_chr_unique() {
        let batch = BsxBatch::try_from(
            df![
                "chr" => ["1", "2", "3"],
                "position" => [1, 2, 3],
                "strand" => [".", ".", "."],
                "context" => ["CG", "CHG", "CHH"],
                "count_m" => [0,0,0],
                "count_total" => [0,0,0],
                "density" => [0, 1, 0]
            ].unwrap().lazy().cast(BsxBatch::hashmap(), true).collect().unwrap(),
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

        assert_eq!(
            batch.first_position().unwrap(),
            encoded.first_position().unwrap()
        );
        assert_eq!(
            batch.last_position().unwrap(),
            encoded.last_position().unwrap()
        );
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
        let batch = BsxBatch::try_from(
            df![
                "chr" => ["1", "1", "1", "1", "1", "1", "1"],
                "position" => [2, 3, 4, 5, 7, 9, 10],
                "strand" => [".", ".", ".", ".", ".", ".", "."],
                "context" => ["CG", "CHG", "CHH", "CG", "CHG", "CHH", "CG"],
                "count_m" => [0,0,0,0,0,0,0],
                "count_total" => [0,0,0,0,0,0,0],
                "density" => [0, 1, 0, 0, 1, 0, 0]
            ].unwrap().lazy().cast(BsxBatch::hashmap(), true).collect().unwrap(),
        ).unwrap();
        let chr_dtype = get_categorical_dtype(vec!["1".into()]);
        let encoded = EncodedBsxBatch::encode(batch.clone(), &chr_dtype).unwrap();

        let trimmed = encoded.trim_region(&RegionCoordinates::new(
            "1".to_string(),
            4,
            9,
            Strand::None,
        )).unwrap();
        assert_eq!(
            trimmed.first_position().unwrap().position(),
            4
        );
        assert_eq!(
            trimmed.last_position().unwrap().position(),
            9
        );

        let trimmed = encoded.trim_region(&RegionCoordinates::new(
            "1".to_string(),
            1,
            9,
            Strand::None,
        )).unwrap();
        assert_eq!(
            trimmed.first_position().unwrap().position(),
            2
        );
        assert_eq!(
            trimmed.last_position().unwrap().position(),
            9
        );

        let trimmed = encoded.trim_region(&RegionCoordinates::new(
            "1".to_string(),
            4,
            12,
            Strand::None,
        )).unwrap();
        assert_eq!(
            trimmed.first_position().unwrap().position(),
            4
        );
        assert_eq!(
            trimmed.last_position().unwrap().position(),
            10
        );
    }
}
