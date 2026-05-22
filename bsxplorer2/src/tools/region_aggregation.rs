use std::cmp::Ordering;

use anyhow::{bail, Result};

use crate::data_structs::batch::BsxBatch;
use crate::data_structs::coords::Contig;
use crate::data_structs::typedef::{BsxSmallStr, PosType};
use crate::data_structs::{Context, Strand};
use crate::io::bsx::RegionReader;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionSpec {
    pub region_id: String,
    pub seqname: String,
    pub start: PosType,
    pub end: PosType,
    pub strand: Option<Strand>,
    pub context: Option<Context>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandPolicy {
    Both,
    Plus,
    Minus,
    RegionStrand,
}

impl StrandPolicy {
    pub fn parse(value: &str) -> Result<Self> {
        match value.to_ascii_lowercase().as_str() {
            "both" | "ignore" => Ok(Self::Both),
            "plus" | "+" => Ok(Self::Plus),
            "minus" | "-" => Ok(Self::Minus),
            "region_strand" | "same" => Ok(Self::RegionStrand),
            other => bail!("Unsupported strand policy: {other}"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AggregationOptions {
    pub context: Option<Context>,
    pub strand_policy: StrandPolicy,
    pub min_total: u32,
    pub chunk_size: usize,
    pub include_empty_regions: bool,
}

impl Default for AggregationOptions {
    fn default() -> Self {
        Self {
            context: None,
            strand_policy: StrandPolicy::Both,
            min_total: 0,
            chunk_size: 10_000,
            include_empty_regions: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RegionCountResult {
    pub region_id: String,
    pub seqname: String,
    pub start: PosType,
    pub end: PosType,
    pub strand: String,
    pub context: String,
    pub sample_id: String,
    pub m_c: u64,
    pub u_c: u64,
    pub total: u64,
    pub n_cytosines: u64,
    pub mean_methylation: Option<f64>,
    pub coverage_qc: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CountRecord {
    pub seqname: String,
    pub position: PosType,
    pub strand: Strand,
    pub context: Context,
    pub count_m: u32,
    pub count_total: u32,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct RegionAccumulator {
    m_c: u64,
    u_c: u64,
    total: u64,
    n_cytosines: u64,
}

impl RegionAccumulator {
    fn add(
        &mut self,
        count_m: u32,
        count_total: u32,
    ) {
        let count_u = count_total.saturating_sub(count_m);
        self.m_c += count_m as u64;
        self.u_c += count_u as u64;
        self.total += (count_m + count_u) as u64;
        self.n_cytosines += 1;
    }
}

pub fn parse_context(value: &str) -> Result<Context> {
    match value.to_ascii_uppercase().as_str() {
        "CG" => Ok(Context::CG),
        "CHG" => Ok(Context::CHG),
        "CHH" => Ok(Context::CHH),
        other => bail!("Unsupported methylation context: {other}"),
    }
}

pub fn parse_optional_context(value: Option<&str>) -> Result<Option<Context>> {
    match value {
        Some(context)
            if !context.is_empty() && !context.eq_ignore_ascii_case("all") =>
        {
            parse_context(context).map(Some)
        },
        _ => Ok(None),
    }
}

pub fn normalize_or_filter_context(
    record_context: Context,
    region_context: Option<Context>,
    options: &AggregationOptions,
) -> bool {
    let requested_context = options.context.or(region_context);
    requested_context
        .map(|context| context == record_context)
        .unwrap_or(true)
}

pub fn apply_strand_policy(
    record_strand: Strand,
    region_strand: Option<Strand>,
    options: &AggregationOptions,
) -> bool {
    match options.strand_policy {
        StrandPolicy::Both => true,
        StrandPolicy::Plus => record_strand == Strand::Forward,
        StrandPolicy::Minus => record_strand == Strand::Reverse,
        StrandPolicy::RegionStrand => region_strand
            .filter(|strand| *strand != Strand::None)
            .map(|strand| strand == record_strand)
            .unwrap_or(true),
    }
}

pub fn aggregate_counts_for_regions(
    reader: &mut RegionReader,
    regions: &[RegionSpec],
    sample_id: Option<&str>,
    options: &AggregationOptions,
) -> Result<Vec<RegionCountResult>> {
    validate_regions(regions)?;

    let mut region_order = (0..regions.len()).collect::<Vec<_>>();
    region_order.sort_by(|left, right| {
        let left_region = &regions[*left];
        let right_region = &regions[*right];
        let left_seq = BsxSmallStr::from(left_region.seqname.as_str());
        let right_seq = BsxSmallStr::from(right_region.seqname.as_str());
        let left_chr_idx = reader
            .index()
            .get_chr_index(&left_seq)
            .unwrap_or(usize::MAX);
        let right_chr_idx = reader
            .index()
            .get_chr_index(&right_seq)
            .unwrap_or(usize::MAX);
        left_chr_idx
            .cmp(&right_chr_idx)
            .then(compare_regions(left_region, right_region))
    });

    let mut results = Vec::with_capacity(regions.len());
    for region_idx in region_order {
        let region = &regions[region_idx];
        let result = if !reader
            .index()
            .get_chr_order()
            .contains(&BsxSmallStr::from(region.seqname.as_str()))
        {
            finalize_region_count_result(
                region,
                sample_id,
                options,
                RegionAccumulator::default(),
            )
        } else {
            let contig = Contig::new(
                BsxSmallStr::from(region.seqname.as_str()),
                region.start,
                region.end.saturating_add(1),
                Strand::None,
            );
            let accumulator = match reader.query(contig, None)? {
                Some(batch) => aggregate_batch(region, &batch, options),
                None => RegionAccumulator::default(),
            };
            finalize_region_count_result(region, sample_id, options, accumulator)
        };

        if options.include_empty_regions || result.n_cytosines > 0 {
            results.push((region_idx, result));
        }
    }

    results.sort_by_key(|(idx, _)| *idx);
    Ok(results.into_iter().map(|(_, result)| result).collect())
}

pub fn aggregate_counts_for_region_chunk(
    records: &[CountRecord],
    regions: &[RegionSpec],
    sample_id: Option<&str>,
    options: &AggregationOptions,
) -> Result<Vec<RegionCountResult>> {
    validate_regions(regions)?;
    let mut results = Vec::with_capacity(regions.len());

    for region in regions {
        let mut accumulator = RegionAccumulator::default();
        for record in records {
            if record.seqname != region.seqname
                || record.position < region.start
                || record.position > region.end
                || record.count_total < options.min_total
                || !normalize_or_filter_context(record.context, region.context, options)
                || !apply_strand_policy(record.strand, region.strand, options)
            {
                continue;
            }
            accumulator.add(record.count_m, record.count_total);
        }

        let result =
            finalize_region_count_result(region, sample_id, options, accumulator);
        if options.include_empty_regions || result.n_cytosines > 0 {
            results.push(result);
        }
    }

    Ok(results)
}

fn validate_regions(regions: &[RegionSpec]) -> Result<()> {
    for region in regions {
        if region.start > region.end {
            bail!(
                "Invalid region coordinates for {}: start {} is greater than end {}",
                region.region_id,
                region.start,
                region.end
            );
        }
    }
    Ok(())
}

fn compare_regions(
    left: &RegionSpec,
    right: &RegionSpec,
) -> Ordering {
    left.seqname
        .cmp(&right.seqname)
        .then(left.start.cmp(&right.start))
        .then(left.end.cmp(&right.end))
}

fn aggregate_batch(
    region: &RegionSpec,
    batch: &BsxBatch,
    options: &AggregationOptions,
) -> RegionAccumulator {
    let mut accumulator = RegionAccumulator::default();

    let positions = batch.position();
    let strands = batch.strand();
    let contexts = batch.context();
    let count_m = batch.count_m();
    let count_total = batch.count_total();

    for idx in 0..batch.len() {
        let Some(position) = positions.get(idx) else {
            continue;
        };
        if position < region.start || position > region.end {
            continue;
        }

        let Some(total) = count_total.get(idx) else {
            continue;
        };
        if (total as u32) < options.min_total {
            continue;
        }

        let record_context = Context::from(contexts.get(idx));
        if !normalize_or_filter_context(record_context, region.context, options) {
            continue;
        }

        let record_strand = Strand::from(strands.get(idx));
        if !apply_strand_policy(record_strand, region.strand, options) {
            continue;
        }

        let methylated = count_m.get(idx).unwrap_or(0);
        accumulator.add(methylated as u32, total as u32);
    }

    accumulator
}

pub fn finalize_region_count_result(
    region: &RegionSpec,
    sample_id: Option<&str>,
    options: &AggregationOptions,
    accumulator: RegionAccumulator,
) -> RegionCountResult {
    let mean_methylation = if accumulator.total > 0 {
        Some(accumulator.m_c as f64 / accumulator.total as f64)
    } else {
        None
    };

    let coverage_qc = if accumulator.n_cytosines == 0 {
        "no_records"
    } else if accumulator.total == 0 {
        "zero_coverage"
    } else if accumulator.total < options.min_total as u64 {
        "low_coverage"
    } else {
        "ok"
    };

    RegionCountResult {
        region_id: region.region_id.clone(),
        seqname: region.seqname.clone(),
        start: region.start,
        end: region.end,
        strand: region
            .strand
            .map(|strand| strand.to_string())
            .unwrap_or_else(|| ".".to_string()),
        context: options
            .context
            .or(region.context)
            .map(|context| context.to_string())
            .unwrap_or_else(|| "all".to_string()),
        sample_id: sample_id.unwrap_or("").to_string(),
        m_c: accumulator.m_c,
        u_c: accumulator.u_c,
        total: accumulator.total,
        n_cytosines: accumulator.n_cytosines,
        mean_methylation,
        coverage_qc: coverage_qc.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_options() -> AggregationOptions {
        AggregationOptions {
            context: None,
            strand_policy: StrandPolicy::Both,
            min_total: 0,
            chunk_size: 2,
            include_empty_regions: true,
        }
    }

    fn records() -> Vec<CountRecord> {
        vec![
            CountRecord {
                seqname: "chr1".to_string(),
                position: 10,
                strand: Strand::Forward,
                context: Context::CG,
                count_m: 4,
                count_total: 10,
            },
            CountRecord {
                seqname: "chr1".to_string(),
                position: 20,
                strand: Strand::Reverse,
                context: Context::CHG,
                count_m: 2,
                count_total: 5,
            },
            CountRecord {
                seqname: "chr1".to_string(),
                position: 30,
                strand: Strand::Forward,
                context: Context::CG,
                count_m: 7,
                count_total: 10,
            },
        ]
    }

    fn region(
        start: PosType,
        end: PosType,
    ) -> RegionSpec {
        RegionSpec {
            region_id: "r1".to_string(),
            seqname: "chr1".to_string(),
            start,
            end,
            strand: None,
            context: None,
        }
    }

    #[test]
    fn aggregates_one_region_exact_sums() {
        let results = aggregate_counts_for_region_chunk(
            &records(),
            &[region(1, 25)],
            Some("sample1"),
            &default_options(),
        )
        .unwrap();

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].m_c, 6);
        assert_eq!(results[0].u_c, 9);
        assert_eq!(results[0].total, 15);
        assert_eq!(results[0].n_cytosines, 2);
    }

    #[test]
    fn overlapping_regions_are_independent() {
        let mut r2 = region(15, 35);
        r2.region_id = "r2".to_string();

        let results = aggregate_counts_for_region_chunk(
            &records(),
            &[region(1, 25), r2],
            None,
            &default_options(),
        )
        .unwrap();

        assert_eq!(results[0].m_c, 6);
        assert_eq!(results[1].m_c, 9);
    }

    #[test]
    fn context_filter_works() {
        let mut options = default_options();
        options.context = Some(Context::CG);

        let results = aggregate_counts_for_region_chunk(
            &records(),
            &[region(1, 40)],
            None,
            &options,
        )
        .unwrap();

        assert_eq!(results[0].m_c, 11);
        assert_eq!(results[0].n_cytosines, 2);
        assert_eq!(results[0].context, "CG");
    }

    #[test]
    fn strand_filter_works() {
        let mut options = default_options();
        options.strand_policy = StrandPolicy::Plus;

        let results = aggregate_counts_for_region_chunk(
            &records(),
            &[region(1, 40)],
            None,
            &options,
        )
        .unwrap();

        assert_eq!(results[0].m_c, 11);
        assert_eq!(results[0].n_cytosines, 2);
    }

    #[test]
    fn empty_region_is_returned_when_requested() {
        let results = aggregate_counts_for_region_chunk(
            &records(),
            &[region(100, 200)],
            None,
            &default_options(),
        )
        .unwrap();

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].coverage_qc, "no_records");
        assert_eq!(results[0].total, 0);
    }
}
