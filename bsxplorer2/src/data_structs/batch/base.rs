use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::{
    Deref,
    Index,
};
use std::panic::AssertUnwindSafe;
use std::str::FromStr;

use anyhow::bail;
use itertools::{
    izip,
    Itertools,
};
use num::Zero;
use polars::frame::column::ScalarColumn;
use polars::prelude::*;

use super::{
    create_empty_categorical_dtype,
    create_empty_series,
    get_col_fn,
    BsxColumns as BsxCol,
};
use crate::data_structs::typedef::{
    CountType,
    DensityType,
    PosType,
};
use crate::plsmallstr;
use crate::prelude::*;
#[cfg(feature = "tools")]
use crate::tools::dimred::*;

/// Represents a batch of methylation data as a Polars DataFrame.
///
/// This struct wraps a Polars DataFrame with specific columns expected for
/// methylation data (chromosome, position, strand, context, counts, density).
/// It provides methods for accessing columns, manipulating data, and performing
/// specific methylation-related operations.
#[derive(Debug, Clone, PartialEq)]
pub struct BsxBatch {
    data: DataFrame,
}

impl Display for BsxBatch {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        writeln!(f, "BsxBatch with {} rows", self.len())?;
        if !self.is_empty() {
            writeln!(f, "Covers region: {}", self.as_contig().unwrap())?;
            writeln!(f, "Data:")?;
            writeln!(f, "{}", self.data())?;
        }
        Ok(())
    }
}

impl Index<BsxCol> for BsxBatch {
    type Output = Series;

    fn index(
        &self,
        index: BsxCol,
    ) -> &Self::Output {
        self.column(index)
    }
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
    /// Creates a new `BsxBatch` from a `DataFrame` without validation.
    ///
    /// # Safety
    /// The caller must ensure that the DataFrame is valid, containing the
    /// expected columns with correct data types and adhering to internal
    /// invariants (e.g., sorted positions, single chromosome).
    #[inline(always)]
    pub unsafe fn new_unchecked(df: DataFrame) -> Self {
        BsxBatch { data: df }
    }

    /// Tries to create a `BsxBatch` from individual column vectors.
    ///
    /// Validates that all input vectors have the same length and calculates the
    /// density. Creates a DataFrame and wraps it in a `BsxBatch`.
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

    /// Creates an empty `BsxBatch` with the correct schema.
    ///
    /// Allows specifying the dtype for the chromosome column.
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
    /// Returns a reference to the underlying Polars DataFrame.
    #[inline(always)]
    pub fn data(&self) -> &DataFrame {
        &self.data
    }

    /// Consumes the `BsxBatch` and returns the underlying Polars DataFrame.
    #[inline(always)]
    pub fn into_inner(self) -> DataFrame {
        self.data
    }

    /// Converts the `BsxBatch` into a `LazyBsxBatch`.
    #[inline(always)]
    pub fn lazy(self) -> LazyBsxBatch {
        self.into()
    }

    /// Gets a specific column by its enum identifier.
    pub fn column(
        &self,
        name: BsxCol,
    ) -> &Series {
        self.data()
            .column(name.as_str())
            .unwrap()
            .as_materialized_series()
    }

    /// Checks if the `BsxBatch` is empty (contains no rows).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    // OPERATIONS
    /// Splits the `BsxBatch` into two `BsxBatch`es at the given index.
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

    /// Rechunks the underlying DataFrame in place.
    pub fn rechunk(&mut self) {
        self.data.rechunk_mut();
    }

    /// Appends another `BsxBatch` to this one, performing validation.
    ///
    /// This method modifies the batch in-place. Additional checks for
    /// duplicated positions, correct positions order, and single chromosome
    /// are performed after extension.
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

    /// Appends another `BsxBatch` to this one without performing validation.
    ///
    /// # Safety
    /// The caller must ensure that appending `other` results in a valid
    /// `BsxBatch`.
    pub unsafe fn extend_unchecked(
        &mut self,
        other: &Self,
    ) {
        self.data.extend(other.data()).unwrap()
    }

    /// Adds context data (strand and context) to the `BsxBatch`.
    ///
    /// Performs a left join with the context data based on position.
    /// Fills nulls for count_m, count_total, and density for positions
    /// present in context_data but not in the original batch.
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

    /// Creates a slice of the `BsxBatch`.
    #[allow(unsafe_code)]
    pub fn slice(
        &self,
        start: i64,
        length: usize,
    ) -> Self {
        let slice = self.data().slice(start, length);
        unsafe { Self::new_unchecked(slice) }
    }

    /// Converts the `BsxBatch` into a DataFrame formatted according to the
    /// specified `ReportType`.
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

    /// Calculates and returns various methylation statistics for the batch.
    pub fn get_methylation_stats(&self) -> RegionMethAgg {
        izip!(
            self.context(),
            self.strand(),
            self.count_m(),
            self.count_total(),
            self.density()
        )
        .into_iter()
        .filter(|(_, _, _, c, d)| c.is_some() || d.is_some())
        .map(|(ctx, snd, cm, ct, dens)| {
            let count = if cm.is_some() { ct } else { None };
            let sum = cm.map(|v| v as DensityType).unwrap_or(dens.unwrap());
            (ctx, snd, sum, count)
        })
        .fold(RegionMethAgg::new(), |mut acc, (ctx, snd, sum, count)| {
            acc.add_cytosine(sum, count, ctx.into(), snd.into());
            acc
        })
    }

    /// Filters sites based on a binomial test p-value and replaces counts and
    /// density.
    ///
    /// Keeps sites where the methylation p-value is less than or equal to the
    /// given `pvalue`. Sites failing the test have counts set to (0, 0) and
    /// density to NaN. Sites passing the test have counts set to (1, 1) and
    /// density to 1.0. Sites with 0 total reads initially have density NaN
    /// and counts (0, 0).
    #[cfg(feature = "tools")]
    pub fn as_binom(
        self,
        mean: f64,
        pvalue: f64,
    ) -> PolarsResult<Self> {
        use statrs::distribution::{
            Binomial,
            DiscreteCDF,
        };
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
                    let binom = Binomial::new(mean, n as u64).unwrap();
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

    /// Segments the methylation data using a specified algorithm and aggregates
    /// densities.
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

        let mut segment_boundaries = match method {
            SegmentAlgorithm::Pelt(beta, min_size) => {
                let (segments, _score) = pelt(
                    &meth_data,
                    beta.unwrap_or((meth_data.len() as f64).ln()),
                    min_size,
                );
                segments.iter().map(|v| v + 1).collect_vec()
            },
        };
        if segment_boundaries.is_empty() {
            segment_boundaries.push(self.len());
        }

        self.partition(segment_boundaries, AggMethod::Mean.get_fn())
    }

    /// Discretises the batch into a fixed number of fragments based on genomic
    /// position and aggregates densities.
    ///
    /// Divides the genomic range covered by the batch into `n_fragments`
    /// roughly equal-sized bins and calculates an aggregate density for
    /// sites falling into each bin.
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

    /// Partitions the batch based on the provided breakpoints and aggregates
    /// densities within each partition.
    ///
    /// Breakpoints are indices into the sorted list of positions. The function
    /// aggregates the density values within each segment defined by
    /// consecutive breakpoints (and the start/end of the batch).
    /// Returns a tuple of vectors: relative genomic positions (0.0 to 1.0) of
    /// the segment ends, and the aggregated density for each segment.
    #[allow(clippy::type_complexity)]
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

        let start = self.first_pos().unwrap();
        let end = self.last_pos().unwrap();
        let length = (end + 1 - start) as f64;

        if breakpoints.last().unwrap() != &size {
            breakpoints.push(size);
        }

        let rel_positions = unsafe {
            let positions = self
                .position()
                .to_vec_null_aware()
                .left()
                .unwrap_unchecked();
            let mut boundary_ends = (0..(breakpoints.len() - 1))
                .map(|i| positions[breakpoints[i]])
                .map(|pos| (pos - start) as f64 / length)
                .collect_vec();
            boundary_ends.push(1.0);
            boundary_ends
        };

        let densities = self.partition_density(&breakpoints, agg_fn)?;
        Ok((rel_positions, densities))
    }

    #[allow(clippy::type_complexity)]
    pub fn partition_density(
        &self,
        breakpoints: &[usize],
        agg_fn: Box<dyn Fn(&[f32]) -> f64>,
    ) -> anyhow::Result<Vec<f64>> {
        if self.is_empty() {
            return Ok(Default::default());
        }
        let size = self.len();
        if breakpoints.iter().any(|v| v > &size) {
            bail!("Partition index out of bounds")
        }
        if !breakpoints.is_empty() && breakpoints.first().unwrap() == &0 {
            bail!("Partition index can not be 0")
        }
        if self.density().null_count() > 0 {
            bail!("Density contains nulls")
        }

        let mut breakpoints_ext = breakpoints.to_vec();
        breakpoints_ext.insert(0, 0);

        if breakpoints_ext.last().unwrap() != &size {
            breakpoints_ext.push(size);
        }

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
        Ok(densities)
    }

    // POSITION
    /// Gets the position of the first site in the batch.
    #[inline]
    pub fn first_pos(&self) -> Option<u32> {
        if self.data.is_empty() {
            return None;
        }
        let result =
            std::panic::catch_unwind(AssertUnwindSafe(|| self.position().first()));
        if result.is_err() {
            panic!("Polars internal error. Make sure to rechunk the data!")
        }
        else {
            result.unwrap()
        }
    }

    /// Gets the position of the last site in the batch.
    #[inline]
    pub fn last_pos(&self) -> Option<u32> {
        if self.data.is_empty() {
            return None;
        }
        let result =
            std::panic::catch_unwind(AssertUnwindSafe(|| self.position().last()));
        if result.is_err() {
            panic!("Polars internal error. Make sure to rechunk the data!")
        }
        else {
            result.unwrap()
        }
    }

    pub fn positions_vec(&self) -> Vec<PosType> {
        unsafe {
            self.position()
                .to_vec_null_aware()
                .left()
                .unwrap_unchecked()
        }
    }

    /// Gets the sequence name (chromosome) for the batch.
    ///
    /// Assumes all sites in the batch are on the same sequence.
    pub fn seqname(&self) -> Option<&str> {
        if self.data.is_empty() {
            return None;
        }
        self.chr().iter_str().next().flatten()
    }

    /// Gets the number of rows (sites) in the batch.
    #[inline]
    pub fn len(&self) -> usize {
        self.data().height()
    }

    /// Gets the genomic position of the first site in the batch as a
    /// `GenomicPosition`.
    pub fn first_genomic_pos(&self) -> Option<GenomicPosition> {
        let seqname = self.seqname().map(BsxSmallStr::from);
        let pos = self.first_pos();
        seqname.and_then(|seqname| pos.map(|pos| GenomicPosition::new(seqname, pos)))
    }

    /// Gets the genomic position of the last site in the batch as a
    /// `GenomicPosition`.
    pub fn last_genomic_pos(&self) -> Option<GenomicPosition> {
        let seqname = self.seqname().map(BsxSmallStr::from);
        let pos = self.last_pos();
        seqname.and_then(|seqname| pos.map(|pos| GenomicPosition::new(seqname, pos)))
    }

    /// Represents the batch as a `Contig` if possible.
    ///
    /// Returns a `Contig` if the batch is not empty and has a defined
    /// sequence name, first position, and last position.
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

    pub fn set_validity_bitmap(
        &mut self,
        column: BsxCol,
        bitmap: Vec<bool>,
    ) -> PolarsResult<()> {
        if matches!(column, BsxCol::Chr | BsxCol::Position) {
            return Err(PolarsError::InvalidOperation(
                format!("Column {} cannot be set to null", column).into(),
            ));
        }
        if bitmap.len() != self.len() {
            return Err(PolarsError::OutOfBounds(
                "Lengths of data and bitmap do not match".into(),
            ));
        }
        self.data.apply(column.as_ref(), |c| {
            let values = c
                .as_materialized_series()
                .iter()
                .zip(bitmap.iter())
                .map(|(v, valid)| {
                    if *valid {
                        v
                    }
                    else {
                        AnyValue::Null
                    }
                })
                .collect_vec();
            Series::new(column.as_ref().into(), &values)
        })?;
        Ok(())
    }
}

mod report_type_conversion {
    use super::*;

    /// Converts a lazy DataFrame representation of BsxBatch data to BedGraph
    /// format.
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

    /// Converts a lazy DataFrame representation of BsxBatch data to Coverage
    /// format.
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

    /// Converts a lazy DataFrame representation of BsxBatch data to Bismark
    /// format.
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

    /// Converts a lazy DataFrame representation of BsxBatch data to CGmap
    /// format.
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

/// Enumerates different aggregation methods for methylation density values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AggMethod {
    Mean,
    Sum,
    GeometricMean,
    Median,
    Max,
    Min,
}

impl FromStr for AggMethod {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "mean" => Ok(AggMethod::Mean),
            "sum" => Ok(AggMethod::Sum),
            "gmean" => Ok(AggMethod::GeometricMean),
            "median" => Ok(AggMethod::Median),
            "max" => Ok(AggMethod::Max),
            "min" => Ok(AggMethod::Min),
            other => bail!("AggMethod {} not implemented", other),
        }
    }
}

impl AggMethod {
    /// Returns a boxed closure for the selected aggregation method.
    #[allow(clippy::type_complexity)]
    pub fn get_fn(&self) -> Box<dyn Fn(&[f32]) -> f64 + Sync + Send> {
        Box::new(match self {
            AggMethod::Mean => {
                |arr: &[f32]| {
                    let len = arr.len();
                    arr.iter().sum::<f32>() as f64 / len as f64
                }
            },
            AggMethod::Sum => |arr: &[f32]| arr.iter().sum::<f32>() as f64,
            AggMethod::GeometricMean => {
                |arr: &[f32]| {
                    let len = arr.len();
                    if arr.is_empty() {
                        f64::NAN
                    }
                    else {
                        (arr.iter().map(|v| (*v as f64).ln()).sum::<f64>() / len as f64)
                            .exp()
                    }
                }
            },
            #[cfg(feature = "tools")]
            AggMethod::Median => {
                |arr: &[f32]| {
                    use statrs::statistics::*;
                    if arr.is_empty() {
                        f64::NAN
                    }
                    else {
                        Data::new(arr.iter().map(|v| *v as f64).collect_vec())
                            .percentile(50)
                    }
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

    /// Returns a function that can aggregate columns using Polars expressions.
    #[allow(clippy::type_complexity)]
    pub fn get_expr(&self) -> Box<dyn Fn(Vec<&Column>) -> Column> {
        Box::new(match self {
            AggMethod::Mean => {
                |columns: Vec<&Column>| {
                    if columns.is_empty() {
                        return Series::new_empty("".into(), &DataType::Float64).into();
                    }
                    let sum: Column = columns
                        .iter()
                        .map(|col| {
                            col.as_materialized_series()
                                .cast(&DataType::Float64)
                                .unwrap()
                        })
                        .reduce(|acc, series| (&acc + &series).unwrap())
                        .unwrap()
                        .into();
                    let count = columns.len() as f64;
                    &sum / count
                }
            },
            AggMethod::Sum => {
                |columns: Vec<&Column>| {
                    if columns.is_empty() {
                        return Series::new_empty("".into(), &DataType::Float64).into();
                    }
                    columns
                        .iter()
                        .map(|col| {
                            col.as_materialized_series()
                                .cast(&DataType::Float64)
                                .unwrap()
                        })
                        .reduce(|acc, series| (&acc + &series).unwrap())
                        .unwrap()
                        .into()
                }
            },
            AggMethod::GeometricMean => {
                |columns: Vec<&Column>| {
                    if columns.is_empty() {
                        return Series::new_empty("".into(), &DataType::Float64).into();
                    }
                    let log_sum = columns
                        .iter()
                        .map(|col| {
                            col.as_materialized_series()
                                .cast(&DataType::Float64)
                                .unwrap()
                                .log(std::f64::consts::E)
                        })
                        .reduce(|acc, series| (&acc + &series).unwrap())
                        .unwrap();
                    let count = columns.len() as f64;
                    let log_mean = &log_sum / count;
                    log_mean.exp().into()
                }
            },
            AggMethod::Median => {
                |columns: Vec<&Column>| {
                    if columns.is_empty() {
                        return Series::new_empty("".into(), &DataType::Float64).into();
                    }
                    let height = columns[0].len();
                    let mut values = Vec::with_capacity(height);
                    for i in 0..height {
                        let row_values: Vec<f64> = columns
                            .iter()
                            .map(|col| {
                                col.as_materialized_series()
                                    .cast(&DataType::Float64)
                                    .unwrap()
                                    .get(i)
                                    .unwrap()
                                    .extract::<f64>()
                                    .unwrap()
                            })
                            .collect();
                        let median = if row_values.is_empty() {
                            f64::NAN
                        }
                        else {
                            let mut sorted = row_values;
                            sorted.sort_by(|a, b| {
                                a.partial_cmp(b).unwrap_or(Ordering::Equal)
                            });
                            let len = sorted.len();
                            if len % 2 == 0 {
                                (sorted[len / 2 - 1] + sorted[len / 2]) / 2.0
                            }
                            else {
                                sorted[len / 2]
                            }
                        };
                        values.push(median);
                    }
                    Series::from_vec("".into(), values).into()
                }
            },
            AggMethod::Max => {
                todo!()
            },
            AggMethod::Min => {
                todo!()
            },
        })
    }
}
