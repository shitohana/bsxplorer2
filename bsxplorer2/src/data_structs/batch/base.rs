use std::cmp::Ordering;
use std::ops::Deref;

use anyhow::bail;
use hashbrown::HashMap;
use itertools::{izip, Itertools};
use num::Zero;
use polars::frame::column::ScalarColumn;
use polars::prelude::*;
use statrs::distribution::DiscreteCDF;
use statrs::statistics::Statistics;

use super::lazy::LazyBsxBatch;
use super::{create_empty_categorical_dtype,
            create_empty_series,
            get_col_fn,
            BsxBatchBuilder,
            BsxColumns as BsxCol};
use crate::data_structs::context_data::ContextData;
use crate::data_structs::coords::{Contig, GenomicPosition};
use crate::data_structs::enums::{Context, IPCEncodedEnum, Strand};
use crate::data_structs::methstats::MethylationStats;
use crate::data_structs::typedef::{BsxSmallStr, CountType, DensityType, PosType};
use crate::io::report::ReportType;
use crate::plsmallstr;
#[cfg(feature = "tools")]
use crate::tools::dimred::*;

#[derive(Debug, Clone, PartialEq)]
pub struct BsxBatch {
    data: DataFrame,
}

impl BsxBatch {
    // COLUMN GETTERS
    get_col_fn!(chr, BsxCol::Chr.as_str(), categorical, CategoricalChunked);

    get_col_fn!(position, BsxCol::Position.as_str(), u32, UInt32Chunked);

    get_col_fn!(strand, BsxCol::Strand.as_str(), bool, BooleanChunked);

    get_col_fn!(context, BsxCol::Context.as_str(), bool, BooleanChunked);

    get_col_fn!(count_m, BsxCol::CountM.as_str(), u16, UInt16Chunked);

    get_col_fn!(count_total, BsxCol::CountTotal.as_str(), u16, UInt16Chunked);

    get_col_fn!(density, BsxCol::Density.as_str(), f32, Float32Chunked);

    // CONSTRUCTORS
    /// # Safety
    /// The caller must ensure that the DataFrame is valid.
    #[inline(always)]
    pub unsafe fn new_unchecked(df: DataFrame) -> Self {
        BsxBatch { data: df }
    }

    pub fn try_from_columns(
        chr: &str,
        chr_dtype: Option<DataType>,
        positions: Vec<PosType>,
        strand: Vec<bool>,
        context: Vec<Option<bool>>,
        count_m: Vec<CountType>,
        count_total: Vec<CountType>,
    ) -> PolarsResult<Self> {
        assert!(
            [
                positions.len(),
                strand.len(),
                context.len(),
                count_m.len(),
                count_total.len(),
            ]
            .iter()
            .all_equal(),
            "All input vectors must have the same length"
        );
        let density = count_m
            .iter()
            .zip(count_total.iter())
            .map(|(m, t)| *m as DensityType / *t as DensityType)
            .collect_vec();
        let height = positions.len();
        let df = DataFrame::from_iter([
            {
                let name = plsmallstr!(BsxCol::Chr.as_str());
                let data = Scalar::new(
                    chr_dtype.unwrap_or(create_empty_categorical_dtype()),
                    AnyValue::StringOwned(plsmallstr!(chr)),
                );
                Column::Scalar(ScalarColumn::new(name, data, height))
            },
            BsxCol::Position.create_series(positions)?.into(),
            BsxCol::Strand.create_series(strand)?.into(),
            BsxCol::Context.create_series(context)?.into(),
            BsxCol::CountM.create_series(count_m)?.into(),
            BsxCol::CountTotal.create_series(count_total)?.into(),
            BsxCol::Density.create_series(density)?.into(),
        ]);

        Ok(unsafe { BsxBatch::new_unchecked(df) })
    }

    pub fn empty(chr_dtype: Option<&DataType>) -> Self {
        let df = DataFrame::from_iter([
            Series::new_empty(
                plsmallstr!(BsxCol::Chr.as_str()),
                chr_dtype.unwrap_or(&BsxCol::Chr.dtype()),
            ),
            create_empty_series!(Position),
            create_empty_series!(Strand),
            create_empty_series!(Context),
            create_empty_series!(CountM),
            create_empty_series!(CountTotal),
            create_empty_series!(Density),
        ]);
        unsafe { BsxBatch::new_unchecked(df) }
    }

    // CONVERSION
    #[inline(always)]
    pub fn data(&self) -> &DataFrame {
        &self.data
    }

    #[inline(always)]
    pub fn into_inner(self) -> DataFrame {
        self.data
    }

    #[inline(always)]
    pub fn lazy(self) -> LazyBsxBatch {
        self.into()
    }

    pub fn column(
        &self,
        name: BsxCol,
    ) -> &Series {
        self.data()
            .column(name.as_str())
            .unwrap()
            .as_materialized_series()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    // OPERATIONS
    pub fn split_at(
        &self,
        index: usize,
    ) -> (Self, Self) {
        let (left, right) = self.data().split_at(index as i64);
        #[allow(unsafe_code)]
        unsafe {
            (Self::new_unchecked(left), Self::new_unchecked(right))
        }
    }

    pub fn rechunk(&mut self) {
        self.data.rechunk_mut();
    }

    /// This method modifies the batch in-place. Additional checks for
    /// duplicated positions and the correct positions order are
    /// performed.
    pub fn extend(
        &mut self,
        other: &Self,
    ) -> PolarsResult<()> {
        let mut new = self.data.vstack(other.data())?;
        new.rechunk_mut();
        BsxBatchBuilder::no_checks()
            .with_check_duplicates(true)
            .with_check_sorted(true)
            .with_check_single_chr(true)
            .checks_only(&new)?;
        self.data = new;
        Ok(())
    }

    pub unsafe fn extend_unchecked(
        &mut self,
        other: &Self,
    ) -> PolarsResult<()> {
        self.data.extend(other.data())
    }

    pub fn add_context_data(
        self,
        context_data: ContextData,
    ) -> PolarsResult<Self> {
        let chr = self.seqname().unwrap_or_default().to_owned();
        let context_df = context_data.to_df();

        let self_selected = self
            .data
            .lazy()
            .drop([BsxCol::Strand.col(), BsxCol::Context.col()]);
        let joined = context_df
            .lazy()
            .left_join(
                self_selected,
                BsxCol::Position.col(),
                BsxCol::Position.col(),
            )
            .with_columns([
                BsxCol::Chr
                    .col()
                    .fill_null(lit(chr))
                    .alias(BsxCol::Chr.as_str()),
                BsxCol::CountM
                    .col()
                    .fill_null(lit(0))
                    .alias(BsxCol::CountM.as_str()),
                BsxCol::CountTotal
                    .col()
                    .fill_null(lit(0))
                    .alias(BsxCol::CountTotal.as_str()),
                BsxCol::Density
                    .col()
                    .fill_null(lit(DensityType::NAN))
                    .alias(BsxCol::Density.as_str()),
            ]);
        let res = joined.collect()?.select(BsxCol::colnames()).unwrap();
        Ok(unsafe { Self::new_unchecked(res) })
    }

    #[allow(unsafe_code)]
    pub fn slice(
        &self,
        start: i64,
        length: usize,
    ) -> Self {
        let slice = self.data().slice(start, length);
        unsafe { Self::new_unchecked(slice) }
    }

    pub fn into_report(
        self,
        report_type: ReportType,
    ) -> PolarsResult<DataFrame> {
        let lf = self.data.lazy();
        let mut res = match report_type {
            ReportType::BedGraph => report_type_conversion::bedgraph(lf),
            ReportType::Bismark => report_type_conversion::bismark(lf),
            ReportType::CgMap => report_type_conversion::cgmap(lf),
            ReportType::Coverage => report_type_conversion::coverage(lf),
        }
        .cast(report_type.hashmap(), true)
        .select(
            report_type
                .col_names()
                .iter()
                .map(|s| col(*s))
                .collect_vec(),
        );

        let res_schema = res.collect_schema()?;
        let target_schema = SchemaRef::new(report_type.schema());
        let schemas_equal = res_schema.deref().eq(target_schema.deref());

        if !schemas_equal {
            Err(PolarsError::SchemaMismatch(
                format!("{:?} != {:?}", res_schema, target_schema).into(),
            ))
        }
        else {
            res.collect()
        }
    }

    pub fn get_methylation_stats(&self) -> MethylationStats {
        let nonull = self.density().drop_nulls();
        let mean = nonull
            .mean()
            .map(|v| v as DensityType)
            .unwrap_or(DensityType::NAN);
        let var =
            nonull.into_no_null_iter().map(|x| x as f64).variance() as DensityType;

        MethylationStats::from_data(
            mean,
            var,
            self.get_coverage_dist(),
            self.get_context_stats(),
            self.get_strand_stats(),
        )
    }

    pub fn get_coverage_dist(&self) -> HashMap<CountType, u32> {
        self.count_total()
            .into_iter()
            .filter_map(|v| v.map(|v1| (v1, v1)))
            .into_group_map()
            .into_iter()
            .map(|(k, v)| (k, v.len() as u32))
            .collect()
    }

    /// Returns context -> (sum methylation ratios, total counts)
    pub fn get_context_stats(&self) -> HashMap<Context, (DensityType, u32)> {
        izip!(self.context(), self.density())
            .filter_map(|(k, v)| {
                Option::map(v, |density| (Context::from_bool(k), density))
            })
            .into_group_map()
            .into_iter()
            .map(|(k, v)| (k, (v.iter().sum(), v.len() as u32)))
            .collect()
    }

    /// Returns strand -> (sum methylation ratios, total counts)
    pub fn get_strand_stats(&self) -> HashMap<Strand, (DensityType, u32)> {
        izip!(self.strand(), self.density())
            .filter_map(|(k, v)| {
                Option::map(v, |density| (Strand::from_bool(k), density))
            })
            .into_group_map()
            .into_iter()
            .map(|(k, v)| (k, (v.iter().sum(), v.len() as u32)))
            .collect()
    }

    pub fn as_binom(
        self,
        mean: f64,
        pvalue: f64,
    ) -> PolarsResult<Self> {
        let pvalue_vec = izip!(self.count_m(), self.count_total())
            .map(|(m, n)| (m.unwrap_or(0), n.unwrap_or(0)))
            .map(|(m, n)| {
                if n.is_zero() {
                    f64::NAN
                }
                else if m.is_zero() {
                    1.0
                }
                else if n == m {
                    0.0
                }
                else {
                    let binom =
                        statrs::distribution::Binomial::new(mean, n as u64).unwrap();
                    1.0 - binom.cdf(m as u64)
                }
            })
            .collect_vec();

        let binom_count_m: Column = Series::from_any_values(
            plsmallstr!(BsxCol::CountM.as_str()),
            &pvalue_vec
                .iter()
                .cloned()
                .map(|p| {
                    if 0.0 <= p && p <= pvalue {
                        1
                    }
                    else {
                        0
                    }
                })
                .map(AnyValue::UInt16)
                .collect_vec(),
            true,
        )?
        .into();

        let binom_count_total: Column = Series::from_any_values(
            plsmallstr!(BsxCol::CountTotal.as_str()),
            &pvalue_vec
                .iter()
                .cloned()
                .map(|p| {
                    if !p.is_nan() {
                        1
                    }
                    else {
                        0
                    }
                })
                .map(AnyValue::UInt16)
                .collect_vec(),
            true,
        )?
        .into();

        let density = (&binom_count_m.cast(&DataType::Float32)?
            / &binom_count_total.cast(&DataType::Float32)?)?
            .with_name(plsmallstr!(BsxCol::Density.as_str()));

        let mut new_data = self.data.clone();
        new_data.with_column(binom_count_m)?;
        new_data.with_column(binom_count_total)?;
        new_data.with_column(density)?;

        Ok(unsafe { Self::new_unchecked(new_data) })
    }

    #[cfg(feature = "tools")]
    pub fn shrink(
        &self,
        method: SegmentAlgorithm,
    ) -> anyhow::Result<(Vec<f64>, Vec<f64>)> {
        if self.is_empty() {
            return Ok(Default::default());
        }

        if self.count_m().null_count() != 0 || self.count_total().null_count() != 0 {
            anyhow::bail!("Can not segment region. Missing counts data");
        }
        let (count_m, count_total) = unsafe {
            (
                self.count_m().to_vec_null_aware().left().unwrap_unchecked(),
                self.count_total()
                    .to_vec_null_aware()
                    .left()
                    .unwrap_unchecked(),
            )
        };
        let meth_data = MethDataBinom::new(&count_m, &count_total);

        let segment_boundaries = match method {
            SegmentAlgorithm::Pelt(beta, min_size) => {
                let (segments, _score) = pelt(
                    &meth_data,
                    beta.unwrap_or((meth_data.len() as f64).ln()),
                    min_size,
                );
                segments.iter().map(|v| v + 1).collect_vec()
            },
        };

        self.partition(segment_boundaries, AggMethod::Mean.get_fn())
    }

    pub fn discretise(
        &self,
        n_fragments: usize,
        agg_fn: AggMethod,
    ) -> anyhow::Result<(Vec<f64>, Vec<f64>)> {
        // Check valid number of fragments
        if n_fragments == 0 {
            bail!("Cannot partition batch into 0 fragments")
        }
        if self.is_empty() {
            return Ok(Default::default());
        }
        if n_fragments == 1 {
            return self.partition(vec![self.len()], agg_fn.get_fn());
        }

        let positions_series = self.position();
        let positions: Vec<PosType> = positions_series.into_no_null_iter().collect();
        let size = positions.len();

        let start = self.first_pos().unwrap(); // Guaranteed to exist because not empty
        let end = self.last_pos().unwrap(); // Guaranteed to exist because not empty
        let genomic_length = end + 1 - start; // Genomic range length

        let fragment_genomic_length = genomic_length as f64 / n_fragments as f64;
        let mut target_positions: Vec<PosType> = Vec::with_capacity(n_fragments - 1);
        for i in 1..n_fragments {
            let target_pos_f64 = start as f64 + i as f64 * fragment_genomic_length;
            let target_pos = target_pos_f64.round() as PosType;
            target_positions.push(target_pos);
        }

        let mut breakpoints: Vec<usize> = Vec::with_capacity(n_fragments - 1);
        let mut current_search_start_index = 0; // Optimization: positions are sorted, search from last found index

        for target_pos in target_positions {
            // Search in the remaining slice of positions, starting from the last found
            // index.
            let search_slice = &positions[current_search_start_index..];

            // Find the index within the *slice* where the first position >= target_pos
            // occurs.
            let result_in_slice =
                search_slice.binary_search_by(|pos| pos.cmp(&target_pos));

            let breakpoint_index_in_slice = match result_in_slice {
                Ok(idx) => idx, // Found exact match
                Err(idx) => idx, /* No exact match, idx is the insertion point
                                  * (first element >= target) */
            };

            // Convert the index in the slice back to an index in the original
            // 'positions' vector.
            let breakpoint_index =
                current_search_start_index + breakpoint_index_in_slice;
            let final_index = std::cmp::min(breakpoint_index, size);

            breakpoints.push(final_index);

            current_search_start_index = final_index;
        }

        debug_assert_eq!(breakpoints.len(), n_fragments - 1);
        self.partition(breakpoints, agg_fn.get_fn())
    }

    pub fn partition(
        &self,
        mut breakpoints: Vec<usize>,
        agg_fn: Box<dyn Fn(&[f32]) -> f64>,
    ) -> anyhow::Result<(Vec<f64>, Vec<f64>)> {
        if self.is_empty() || breakpoints.is_empty() {
            return Ok(Default::default());
        }
        let size = self.len();
        if breakpoints.iter().any(|v| v > &size) {
            bail!("Partition index out of bounds")
        }
        if breakpoints.first().unwrap() == &0 {
            bail!("Partition index can not be 0")
        }
        if self.density().null_count() > 0 {
            bail!("Density contains nulls")
        }

        let start = self.first_pos().unwrap();
        let end = self.last_pos().unwrap();
        let length = (end + 1 - start) as f64;

        if breakpoints.last().map(|p| p != &size).unwrap_or(true) {
            breakpoints.push(size);
        }

        let rel_positions = unsafe {
            let positions = self.position();
            let mut boundary_ends = (0..(breakpoints.len() - 1))
                .map(|i| positions.get_unchecked(*&breakpoints[i]).unwrap_unchecked())
                .map(|pos| (pos - start) as f64 / length)
                .collect_vec();
            boundary_ends.push(1.0);
            boundary_ends
        };

        let mut breakpoints_ext = breakpoints;
        breakpoints_ext.insert(0, 0);
        let densities = {
            let density_vals = self.density().iter().map(|v| v.unwrap()).collect_vec();
            breakpoints_ext
                .windows(2)
                .map(|window| {
                    let start_i = window[0];
                    let end_i = window[1];
                    agg_fn(&density_vals[start_i..end_i])
                })
                .collect_vec()
        };
        Ok((rel_positions, densities))
    }

    // POSITION
    #[inline]
    pub fn first_pos(&self) -> Option<u32> {
        if self.data.is_empty() {
            return None;
        }
        self.position().first()
    }

    #[inline]
    pub fn last_pos(&self) -> Option<u32> {
        if self.data.is_empty() {
            return None;
        }
        self.position().last()
    }

    pub fn seqname(&self) -> Option<&str> {
        if self.data.is_empty() {
            return None;
        }
        self.chr().iter_str().next().flatten()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.data().height()
    }

    pub fn first_genomic_pos(&self) -> Option<GenomicPosition> {
        let seqname = self.seqname().map(BsxSmallStr::from);
        let pos = self.first_pos();
        seqname.and_then(|seqname| pos.map(|pos| GenomicPosition::new(seqname, pos)))
    }

    pub fn last_genomic_pos(&self) -> Option<GenomicPosition> {
        let seqname = self.seqname().map(BsxSmallStr::from);
        let pos = self.last_pos();
        seqname.and_then(|seqname| pos.map(|pos| GenomicPosition::new(seqname, pos)))
    }

    pub fn as_contig(&self) -> Option<Contig> {
        let seqname = self.seqname().map(BsxSmallStr::from);
        let first = self.first_pos();
        let last = self.last_pos();
        if let (Some(s), Some(f), Some(l)) = (seqname, first, last) {
            Some(Contig::new(s, f, l, Strand::None))
        }
        else {
            None
        }
    }
}

mod report_type_conversion {
    use super::*;

    pub fn bedgraph(lf: LazyFrame) -> LazyFrame {
        lf.select([
            // chr
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            BsxCol::Position.col().alias("start"),
            (BsxCol::Position.col() + lit(1)).alias("end"),
            BsxCol::Density.col().alias("density"),
        ])
        .drop_nans(Some(vec![BsxCol::Density.col()]))
        .cast(ReportType::BedGraph.hashmap(), true)
    }

    pub fn coverage(lf: LazyFrame) -> LazyFrame {
        lf.select([
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            BsxCol::Position.col().alias("start"),
            (BsxCol::Position.col() + lit(1)).alias("end"),
            BsxCol::Density.col().alias("density"),
            BsxCol::CountM.col().alias("count_m"),
            (BsxCol::CountTotal.col() - BsxCol::CountM.col()).alias("count_um"),
        ])
        .drop_nans(Some(vec![BsxCol::Density.col()]))
        .cast(ReportType::Coverage.hashmap(), true)
    }

    pub fn bismark(lf: LazyFrame) -> LazyFrame {
        lf.with_column(
            when(BsxCol::Context.col().is_null())
                .then(lit("CHH"))
                .when(BsxCol::Context.col().eq(lit(false)))
                .then(lit("CHG"))
                .otherwise(lit("CG"))
                .cast(DataType::String)
                .alias(BsxCol::Context.as_str()),
        )
        .select([
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            BsxCol::Position.col().alias("position"),
            // strand
            when(BsxCol::Strand.col().eq(lit(true)))
                .then(lit("+"))
                .when(BsxCol::Strand.col().eq(lit(false)))
                .then(lit("-"))
                .otherwise(lit("."))
                .cast(DataType::String)
                .alias("strand"),
            BsxCol::CountM.col().alias("count_m"),
            (BsxCol::CountTotal.col() - BsxCol::CountM.col()).alias("count_um"),
            BsxCol::Density.col().alias("density"),
            BsxCol::Context.col().alias("context"),
            BsxCol::Context.col().alias("trinuc"),
        ])
    }

    pub fn cgmap(lf: LazyFrame) -> LazyFrame {
        lf.with_column(
            when(BsxCol::Context.col().is_null())
                .then(lit("CHH"))
                .when(BsxCol::Context.col().eq(lit(false)))
                .then(lit("CHG"))
                .otherwise(lit("CG"))
                .cast(DataType::String)
                .alias(BsxCol::Context.as_str()),
        )
        .select([
            // chr
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            // nuc
            when(BsxCol::Strand.col().eq(lit(true)))
                .then(lit("C"))
                .when(BsxCol::Strand.col().eq(lit(false)))
                .then(lit("G"))
                .otherwise(lit("."))
                .cast(DataType::String)
                .alias("nuc"),
            // position
            BsxCol::Position.col().alias("position"),
            // context
            BsxCol::Context.col().alias("context"),
            // dinuc
            BsxCol::Context.col().str().head(lit(2)).alias("dinuc"),
            // density
            BsxCol::Density.col().alias("density"),
            // count_m
            BsxCol::CountM.col().alias("count_m"),
            // count_total
            BsxCol::CountTotal.col().alias("count_total"),
        ])
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AggMethod {
    Mean,
    GeometricMean,
    Median,
    Max,
    Min,
}

impl AggMethod {
    pub(crate) fn get_fn(&self) -> Box<dyn Fn(&[f32]) -> f64> {
        use statrs::statistics::*;

        Box::new(match self {
            AggMethod::Mean => {
                |arr: &[f32]| {
                    let len = arr.len();
                    arr.iter().sum::<f32>() as f64 / len as f64
                }
            },
            AggMethod::GeometricMean => {
                |arr: &[f32]| {
                    let len = arr.len();
                    (arr.iter().map(|v| (*v as f64).ln()).sum::<f64>() / len as f64)
                        .exp()
                }
            },
            AggMethod::Median => {
                |arr: &[f32]| {
                    Data::new(arr.iter().map(|v| *v as f64).collect_vec())
                        .percentile(50)
                }
            },
            AggMethod::Max => {
                |arr: &[f32]| {
                    *arr.iter()
                        .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                        .unwrap_or(&f32::NAN) as f64
                }
            },
            AggMethod::Min => {
                |arr: &[f32]| {
                    *arr.iter()
                        .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                        .unwrap_or(&f32::NAN) as f64
                }
            },
        })
    }
}
