use std::ops::Deref;

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
use crate::data_structs::typedef::BsxSmallStr;
use crate::io::report::ReportType;
use crate::plsmallstr;


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
    #[inline(always)]
    pub unsafe fn new_unchecked(df: DataFrame) -> Self {
        BsxBatch { data: df }
    }

    pub fn try_from_columns(
        chr: &str,
        chr_dtype: Option<DataType>,
        positions: Vec<u32>,
        strand: Vec<bool>,
        context: Vec<Option<bool>>,
        count_m: Vec<u16>,
        count_total: Vec<u16>,
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
            .map(|(m, t)| *m as f32 / *t as f32)
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
    pub fn data(&self) -> &DataFrame {
        &self.data
    }

    pub fn into_inner(self) -> DataFrame {
        self.data
    }

    #[inline(always)]
    pub fn lazy(self) -> LazyBsxBatch {
        self.into()
    }

    pub fn column(
        &self,
        name: &str,
    ) -> Option<&Series> {
        if BsxCol::has_name(name) {
            Some(self.data().column(name).unwrap().as_materialized_series())
        }
        else {
            None
        }
    }

    pub fn is_empty(&self) -> bool {
        self.data().is_empty()
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

    pub fn extend(
        &mut self,
        other: &Self,
    ) -> PolarsResult<()> {
        let mut new = self.data.vstack(other.data())?;
        new.rechunk_mut();
        BsxBatchBuilder::no_checks()
            .with_check_duplicates(true)
            .with_check_sorted(true)
            .checks_only(&new)?;
        self.data = new;
        Ok(())
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
                    .fill_null(lit(f32::NAN))
                    .alias(BsxCol::Density.as_str()),
            ]);
        let res = joined.collect()?;
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
                .into_iter()
                .map(|s| col(*s))
                .collect_vec(),
        );

        let res_schema = res.collect_schema()?;
        let target_schema = SchemaRef::new(report_type.schema());
        let schemas_equal = res_schema.deref().eq(&target_schema.deref());

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
        let mean = nonull.mean().unwrap_or(f64::NAN);
        let var = nonull.into_no_null_iter().map(|x| x as f64).variance();

        MethylationStats::from_data(
            mean,
            var,
            self.get_coverage_dist(),
            self.get_context_stats(),
            self.get_strand_stats(),
        )
    }

    pub fn get_coverage_dist(&self) -> HashMap<u16, u32> {
        izip!(self.count_m(), self.count_total())
            .filter_map(|(k, v)| Option::zip(k, v))
            .into_group_map()
            .into_iter()
            .map(|(k, v)| (k, v.iter().count() as u32))
            .collect()
    }

    /// Returns context -> (sum methylation ratios, total counts)
    pub fn get_context_stats(&self) -> HashMap<Context, (f64, u32)> {
        izip!(self.context(), self.density())
            .filter_map(|(k, v)| {
                Option::map(v, |density| (Context::from_bool(k), density as f64))
            })
            .into_group_map()
            .into_iter()
            .map(|(k, v)| (k, (v.iter().sum(), v.len() as u32)))
            .collect()
    }

    /// Returns strand -> (sum methylation ratios, total counts)
    pub fn get_strand_stats(&self) -> HashMap<Strand, (f64, u32)> {
        izip!(self.strand(), self.density())
            .filter_map(|(k, v)| {
                Option::map(v, |density| (Strand::from_bool(k), density as f64))
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
                .map(|v| AnyValue::UInt16(v))
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
                .map(|v| AnyValue::UInt16(v))
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

    // POSITION
    #[inline]
    pub fn first_pos(&self) -> Option<u32> {
        self.position().first()
    }

    #[inline]
    pub fn last_pos(&self) -> Option<u32> {
        self.position().last()
    }

    pub fn seqname(&self) -> Option<&str> {
        self.chr().iter_str().next().flatten()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.data().height()
    }

    pub fn first_genomic_pos(&self) -> Option<GenomicPosition<BsxSmallStr, u32>> {
        let seqname = self.seqname().map(BsxSmallStr::from);
        let pos = self.first_pos();
        seqname.and_then(|seqname| pos.map(|pos| GenomicPosition::new(seqname, pos)))
    }

    pub fn last_genomic_pos(&self) -> Option<GenomicPosition<BsxSmallStr, u32>> {
        let seqname = self.seqname().map(BsxSmallStr::from);
        let pos = self.last_pos();
        seqname.and_then(|seqname| pos.map(|pos| GenomicPosition::new(seqname, pos)))
    }

    pub fn as_contig(&self) -> Option<Contig<BsxSmallStr, u32>> {
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
            BsxCol::Position.col().alias("end"),
            BsxCol::Density.col().alias("density"),
        ])
        .drop_nans(Some(vec![BsxCol::Density.col()]))
        .cast(ReportType::BedGraph.hashmap(), true)
    }

    pub fn coverage(lf: LazyFrame) -> LazyFrame {
        lf.select([
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            BsxCol::Position.col().alias("start"),
            BsxCol::Position.col().alias("end"),
            BsxCol::CountM.col().alias("count_m"),
            (BsxCol::CountTotal.col() - BsxCol::CountM.col()).alias("count_um"),
            BsxCol::Density.col().alias("density"),
        ])
        .drop_nans(Some(vec![BsxCol::Density.col()]))
        .cast(ReportType::Coverage.hashmap(), true)
    }

    pub fn bismark(lf: LazyFrame) -> LazyFrame {
        lf.with_column(
            when(BsxCol::Context.col() == lit(NULL))
                .then(lit("CHH"))
                .when(BsxCol::Context.col() == lit(false))
                .then(lit("CHG"))
                .otherwise(lit("CG"))
                .cast(DataType::String)
                .alias(BsxCol::Context.as_str()),
        )
        .select([
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            BsxCol::Position.col().alias("position"),
            // strand
            when(BsxCol::Strand.col() == lit(true))
                .then(lit("+"))
                .when(BsxCol::Strand.col() == lit(false))
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
            when(BsxCol::Context.col() == lit(NULL))
                .then(lit("CHH"))
                .when(BsxCol::Context.col() == lit(false))
                .then(lit("CHG"))
                .otherwise(lit("CG"))
                .cast(DataType::String)
                .alias(BsxCol::Context.as_str()),
        )
        .select([
            // chr
            BsxCol::Chr.col().cast(DataType::String).alias("chr"),
            // nuc
            when(BsxCol::Strand.col() == lit(true))
                .then(lit("C"))
                .when(BsxCol::Strand.col() == lit(false))
                .then(lit("G"))
                .otherwise(lit("."))
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

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};

    use super::*;
    use crate::data_structs::batch::create_caregorical_dtype;
    use crate::utils::get_categorical_dtype;

    #[test]
    fn test_empty_batch() {
        let batch = BsxBatch::empty(None);
        assert!(batch.is_empty());
    }

    #[fixture]
    fn test_batch() -> BsxBatch {
        BsxBatch::try_from_columns(
            "chr1",
            Some(create_caregorical_dtype(vec!["chr1".into()])),
            vec![3, 5, 9, 12, 15],
            vec![true, false, true, true, false],
            vec![Some(true), Some(false), Some(true), None, None],
            vec![5, 10, 15, 10, 5],
            vec![10, 30, 20, 10, 10],
        )
        .unwrap()
    }

    #[rstest]
    fn test_column_getters(test_batch: BsxBatch) {
        assert_eq!(test_batch.len(), 5);

        // Test chr column
        let chr = test_batch.chr();
        assert!(chr.iter_str().all(|c| c.unwrap() == "chr1"));

        // Test position column
        let positions = test_batch.position();
        assert_eq!(positions.into_iter().collect::<Vec<_>>(), vec![
            Some(3),
            Some(5),
            Some(9),
            Some(12),
            Some(15)
        ]);

        // Test strand column
        let strands = test_batch.strand();
        assert_eq!(strands.into_iter().collect::<Vec<_>>(), vec![
            Some(true),
            Some(false),
            Some(true),
            Some(true),
            Some(false)
        ]);

        // Test context column
        let contexts = test_batch.context();
        assert_eq!(contexts.into_iter().collect::<Vec<_>>(), vec![
            Some(true),
            Some(false),
            Some(true),
            None,
            None
        ]);

        // Test count_m column
        let count_m = test_batch.count_m();
        assert_eq!(count_m.into_iter().collect::<Vec<_>>(), vec![
            Some(5),
            Some(10),
            Some(15),
            Some(10),
            Some(5)
        ]);

        // Test count_total column
        let count_total = test_batch.count_total();
        assert_eq!(count_total.into_iter().collect::<Vec<_>>(), vec![
            Some(10),
            Some(30),
            Some(20),
            Some(10),
            Some(10)
        ]);

        // Test density column
        let density = test_batch.density();
        assert_eq!(
            density
                .into_iter()
                .map(|v| v.map(|f| (f * 100.0).round() / 100.0))
                .collect::<Vec<_>>(),
            vec![Some(0.5), Some(0.33), Some(0.75), Some(1.0), Some(0.5)]
        );
    }

    #[rstest]
    fn test_split_at(test_batch: BsxBatch) {
        let (left, right) = test_batch.split_at(2);

        assert_eq!(left.len(), 2);
        assert_eq!(right.len(), 3);

        assert_eq!(left.position().into_iter().collect::<Vec<_>>(), vec![
            Some(3),
            Some(5)
        ]);
        assert_eq!(right.position().into_iter().collect::<Vec<_>>(), vec![
            Some(9),
            Some(12),
            Some(15)
        ]);
    }

    #[rstest]
    fn test_slice(test_batch: BsxBatch) {
        let sliced = test_batch.slice(1, 3);

        assert_eq!(sliced.len(), 3);
        assert_eq!(sliced.position().into_iter().collect::<Vec<_>>(), vec![
            Some(5),
            Some(9),
            Some(12)
        ]);
    }

    #[rstest]
    fn test_position_methods(test_batch: BsxBatch) {
        assert_eq!(test_batch.first_pos(), Some(3));
        assert_eq!(test_batch.last_pos(), Some(15));
        assert_eq!(test_batch.seqname(), Some("chr1"));

        let first_genomic_pos = test_batch.first_genomic_pos();
        assert!(first_genomic_pos.is_some());
        assert_eq!(first_genomic_pos.clone().unwrap().position(), 3);
        assert_eq!(first_genomic_pos.clone().unwrap().seqname(), "chr1");

        let last_genomic_pos = test_batch.last_genomic_pos();
        assert!(last_genomic_pos.is_some());
        assert_eq!(last_genomic_pos.clone().unwrap().position(), 15);
        assert_eq!(last_genomic_pos.clone().unwrap().seqname(), "chr1");

        let contig = test_batch.as_contig();
        assert!(contig.is_some());
        assert_eq!(contig.clone().unwrap().start(), 3);
        assert_eq!(contig.clone().unwrap().end(), 15);
        assert_eq!(contig.clone().unwrap().seqname(), "chr1");
    }

    #[rstest]
    fn test_stats_methods(test_batch: BsxBatch) {
        // Test methylation stats
        let meth_stats = test_batch.get_methylation_stats();
        assert!((meth_stats.mean_methylation() - 0.616).abs() < 0.001);

        // Test coverage distribution
        let cov_dist = test_batch.get_coverage_dist();
        assert_eq!(cov_dist.get(&5), Some(&2));
        assert_eq!(cov_dist.get(&10), Some(&2));
        assert_eq!(cov_dist.get(&15), Some(&1));

        // Test context stats
        let context_stats = test_batch.get_context_stats();
        assert_eq!(context_stats.len(), 3);
        let cg_stats = context_stats.get(&Context::CG);
        assert!(cg_stats.is_some());
        let chg_stats = context_stats.get(&Context::CHG);
        assert!(chg_stats.is_some());

        // Test strand stats
        let strand_stats = test_batch.get_strand_stats();
        assert_eq!(strand_stats.len(), 2);
        let pos_stats = strand_stats.get(&Strand::Forward);
        assert!(pos_stats.is_some());
        let neg_stats = strand_stats.get(&Strand::Reverse);
        assert!(neg_stats.is_some());
    }

    #[rstest]
    fn test_as_binom(test_batch: BsxBatch) {
        let binom_batch = test_batch.as_binom(0.5, 0.05).unwrap();
        assert_eq!(binom_batch.len(), 5);

        // Verify that the binom transform changed count_m and density values
        let count_m = binom_batch.count_m();
        let count_total = binom_batch.count_total();
        let density = binom_batch.density();

        // Check types remain the same
        assert_eq!(count_m.dtype(), &DataType::UInt16);
        assert_eq!(count_total.dtype(), &DataType::UInt16);
        assert_eq!(density.dtype(), &DataType::Float32);
    }

    #[rstest]
    fn test_is_empty(test_batch: BsxBatch) {
        assert!(!test_batch.is_empty());
        assert!(BsxBatch::empty(None).is_empty());
    }

    #[rstest]
    #[case::no_chr(None, None)]
    #[case::both_chr(Some(get_categorical_dtype(vec!["chr1".into()])), Some(get_categorical_dtype(vec!["chr1".into()])))]
    #[should_panic]
    #[case::different_types(None, Some(get_categorical_dtype(vec!["chr1".into()])))]
    fn test_can_extend(
        #[case] first_dtype: Option<DataType>,
        #[case] second_dtype: Option<DataType>,
    ) {
        let batch1 = BsxBatch::empty(first_dtype.as_ref());
        let batch2 = BsxBatch::try_from_columns(
            "chr1",
            second_dtype,
            vec![1, 2, 3],
            vec![true, false, true],
            vec![Some(true), Some(false), None],
            vec![1, 2, 3],
            vec![3, 6, 9],
        )
        .unwrap();

        assert!(matches!(
            batch1.column(BsxCol::Chr.as_str()).unwrap().dtype(),
            DataType::Categorical(_, _) | DataType::Enum(_, _)
        ));
        assert!(matches!(
            batch2.column(BsxCol::Chr.as_str()).unwrap().dtype(),
            DataType::Categorical(_, _) | DataType::Enum(_, _)
        ));

        let vstack = batch1.data().vstack(&batch2.data()).unwrap();
        assert!(matches!(
            vstack.column(BsxCol::Chr.as_str()).unwrap().dtype(),
            DataType::Categorical(_, _) | DataType::Enum(_, _)
        ));
        assert!(vstack
            .column(BsxCol::Chr.as_str())
            .unwrap()
            .categorical()
            .unwrap()
            .get_rev_map()
            .find("chr1")
            .is_some());
    }
}
