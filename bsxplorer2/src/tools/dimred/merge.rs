#![allow(unused)]
use std::collections::BTreeSet;
use std::error::Error;
use std::fmt::Debug;
use std::io::Read;
use std::process::exit;
use std::str::FromStr;
use std::time::Instant;

use anyhow::{
    bail,
    ensure,
};
use arcstr::ArcStr;
use bio::bio_types::annot::pos::Pos;
use hashbrown::{
    HashMap,
    HashSet,
};
use itertools::{
    izip,
    Itertools,
};
use polars_arrow::pushable::Pushable;
use rayon::iter::IndexedParallelIterator;
use rayon::prelude::*;
use rust_lapper::{
    Interval,
    Lapper,
};
use spipe::spipe;
use typed_floats::PositiveFinite;

use super::dbscan::{
    self,
    Classification,
};
use crate::data_structs::coords::ContigIntervalMap;
use crate::data_structs::typedef::{
    DensityType,
    PosType,
};
use crate::utils::THREAD_POOL;
use crate::{
    AggMethod,
    BsxSmallStr,
    Contig,
    Strand,
};

pub type EqFloat = PositiveFinite;

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub enum MergeType {
    Full,
    /// DBSCAN algorithm with parameters (eps, mpt, max_dist)
    Dbscan(PosType, usize, PosType, AggMethod),
}

impl ToString for MergeType {
    fn to_string(&self) -> String {
        format!("{self:?}")
    }
}

impl FromStr for MergeType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let split = s.split(":").collect_vec();
        ensure!(split.len() > 0, "Empty string for merge type");
        match split[0] {
            "full" => Ok(Self::Full),
            "dbscan" => {
                ensure!(
                    split.len() == 5,
                    "There should be exactly 4 ':' delimited parameters for 'dbscan'"
                );
                let eps = split[1].parse::<PosType>()?;
                let mpt = split[2].parse::<usize>()?;
                let max_dist = split[3].parse::<PosType>()?;
                let agg_method = split[4].parse::<AggMethod>()?;
                Ok(Self::Dbscan(eps, mpt, max_dist, agg_method))
            },
            other => bail!("Merge method {} not implemented", other),
        }
    }
}

fn validate_chr(
    intervals: &[ContigIntervalMap<EqFloat>]
) -> anyhow::Result<Vec<BsxSmallStr>> {
    if intervals.is_empty() {
        bail!("No data to merge")
    }

    let all_chr_list = intervals.iter().map(ContigIntervalMap::chr_names).concat();
    if all_chr_list.iter().all_unique() {
        bail!("All chromosome names in samples are different")
    }
    Ok(all_chr_list.into_iter().unique().collect_vec())
}

/// This method gathers __all__ intervals as breakpoints for each chrom
fn gather_breakpoints(
    intervals: &[ContigIntervalMap<EqFloat>]
) -> HashMap<ArcStr, Vec<PosType>> {
    intervals
        .into_par_iter()
        .map(ContigIntervalMap::get_breakpoints)
        .collect::<Vec<_>>()
        .concat()
        .into_iter()
        .into_group_map()
        .into_iter()
        .map(|(k, v)| {
            (
                ArcStr::from(k.as_str()),
                v.concat().into_iter().collect_vec(),
            )
        })
        .map(|(k, mut v)| {
            v.par_sort_unstable();
            (k, v)
        })
        .collect::<HashMap<_, _>>()
}

fn partition_by_dist_chr(
    chr: ArcStr,
    pos: Vec<PosType>,
    max_dist: PosType,
) -> Vec<(ArcStr, Vec<u32>)> {
    let partitioned = pos.into_iter().fold(vec![vec![]], |mut acc, new| {
        let last_vec = acc.last_mut().unwrap();
        if last_vec
            .last()
            .map(|v| (*v as u32).abs_diff(new) <= max_dist)
            .unwrap_or(true)
        {
            last_vec.push(new);
        }
        else {
            acc.push(vec![new]);
        }
        acc
    });

    partitioned
        .into_iter()
        .map(|p| (chr.clone(), p))
        .collect_vec()
}

/// Input vectors should be sorted in ascending order
fn partition_by_dist(
    data: HashMap<ArcStr, Vec<PosType>>,
    max_dist: PosType,
) -> Vec<(ArcStr, Vec<u32>)> {
    data.into_par_iter()
        .map(|(k, v)| partition_by_dist_chr(k, v, max_dist))
        .collect::<Vec<_>>()
        .concat()
}

fn apply_dbscan_to_part(
    pos: Vec<PosType>,
    model: dbscan::Model,
    agg_fn: &(dyn Fn(&[f32]) -> f64 + Sync + Send),
) -> Vec<u32> {
    model
        .run(&pos)
        .iter()
        .enumerate()
        .fold(vec![vec![]], |mut acc, (i, v)| {
            match *v {
                Classification::Noise => acc.push(vec![]),
                _ => acc.last_mut().unwrap().push(i),
            }
            acc
        })
        .into_par_iter()
        .filter(|v| !v.is_empty())
        .map(|indices| indices.iter().map(|i| pos[*i] as DensityType).collect_vec())
        .map(|res| agg_fn(&res).round() as PosType)
        .collect::<Vec<_>>()
}

fn apply_dbscan(
    partitioned: Vec<(ArcStr, Vec<u32>)>,
    eps: u32,
    min_points: usize,
    agg: AggMethod,
) -> HashMap<ArcStr, Vec<PosType>> {
    let db = dbscan::Model::new(eps, min_points);
    let agg_fn = agg.get_fn();

    partitioned
        .into_par_iter()
        .map(|(k, v)| (k, apply_dbscan_to_part(v, db.clone(), &agg_fn)))
        .filter(|(_, v)| !v.is_empty())
        .collect::<Vec<_>>()
        .into_iter()
        .into_group_map()
        .into_iter()
        .map(|(k, v)| (k, v.concat()))
        .collect()
}

fn merge_signal(
    new: Vec<PosType>,
    reference: Vec<PosType>,
    reference_values: Vec<EqFloat>,
) {
    let zipped_values = vec![
        vec![-1isize; new.len()],
        (0isize..reference_values.len() as isize).collect_vec(),
    ]
    .concat();
    let zipped = vec![new, reference].concat();

    let sorted_indices = {
        let mut indices = (0..zipped.len()).collect_vec();
        indices.par_sort_unstable_by_key(|i| zipped[*i]);
        indices
    };
}

pub fn preprocess_all_bpoints(
    bpoints: HashMap<ArcStr, Vec<PosType>>
) -> HashMap<ArcStr, Vec<PosType>> {
    bpoints
        .into_iter()
        .map(|(k, mut v)| {
            v.par_sort_unstable();
            v.dedup();
            (k, v)
        })
        .collect()
}

fn interleave_single(
    intervals: &Lapper<PosType, EqFloat>,
    reference: &[PosType],
) -> Vec<(u32, f64)> {
    if intervals.is_empty() {
        return Default::default();
    }
    if reference.len() == 1 {
        let mut prev_end = 0;
        // Calculate mean of all intervals (as there is only one reference point
        // we need just a single value over whole region)
        let acc_value = intervals.iter().fold(0.0, |acc, interval| {
            let out =
                acc + (f64::from(interval.val) * (interval.stop - prev_end) as f64);
            prev_end = interval.stop;
            out
        });
        return vec![(0, acc_value / prev_end as f64)];
    }

    let intervals = intervals.intervals.as_slice();

    let mut cached_pos = 0;
    let mut cached_oft = 0;
    let mut cached_val = 0.0;
    let mut cached_agg = 0.0;
    let mut cached_is_ref = true;

    let mut i_idx = 0;
    let mut r_idx = 0;

    let i_len = intervals.len();
    let r_len = reference.len();

    let mut result = Vec::with_capacity(reference.len() - 1);

    loop {
        if r_idx == r_len {
            break;
        }
        else if i_idx != i_len {
            let ref_pos = unsafe { *reference.get_unchecked(r_idx) };
            let int_pos = unsafe { intervals.get_unchecked(i_idx) }.start;

            let this_is_ref = ref_pos <= int_pos;

            // If both are zero - skip this point, as it is already set
            // Only set the value of the first interval
            if ref_pos == int_pos && ref_pos == 0 {
                cached_val = unsafe { intervals.get_unchecked(i_idx) }.val.into();
                r_idx += 1;
                i_idx += 1;
                continue;
            // Else if positions are equal and not zero reference point
            // is main (i.e. this_is_ref = true). But we update value
            // as well
            }
            else if ref_pos == int_pos {
                cached_val = unsafe { intervals.get_unchecked(i_idx) }.val.into();
            }

            // ref -> ref
            // As there were no interval boundaries, value stays the same,
            // no need to update agg or cached_was ref as well
            if cached_is_ref && this_is_ref {
                result.push((cached_pos, cached_val));

                cached_pos = ref_pos;
                cached_oft = ref_pos;
            // int -> ref
            // End of aggregation - we encountered reference point. Add the
            // remaining points to agg and finalize it.
            // Agg is reset
            }
            else if this_is_ref {
                cached_agg += cached_val * (ref_pos - cached_oft) as f64;
                result.push((cached_pos, cached_agg / (ref_pos - cached_pos) as f64));

                cached_agg = 0.0;
                cached_pos = ref_pos;
                cached_oft = ref_pos;
                cached_is_ref = true;
            // int -> int
            // Another interval boundary, just update cached value and increment
            // agg.
            }
            else if cached_is_ref {
                cached_agg = cached_val * (int_pos - cached_oft) as f64;
                cached_val = unsafe { intervals.get_unchecked(i_idx) }.val.into();
                cached_oft = int_pos;
                cached_is_ref = false;
            // ref -> int
            // Starting an aggregation, update current value (cached_val)
            // and cached_oft;
            }
            else {
                cached_agg += cached_val * (int_pos - cached_oft) as f64;
                cached_val = unsafe { intervals.get_unchecked(i_idx) }.val.into();
                cached_oft = int_pos
            }

            // If reference was greater this time - increment intervals
            if this_is_ref {
                r_idx += 1;
                if int_pos == ref_pos {
                    i_idx += 1;
                }
            // Otherwise increment reference
            }
            else {
                i_idx += 1
            }
        // If there is no more intervals - just save rest of reference points
        // with the same value (because we assume both reference and sample
        // cover the same region)
        }
        else {
            result.push((unsafe { *reference.get_unchecked(r_idx) }, cached_val));
            r_idx += 1;
        }
    }

    debug_assert_eq!(result.len(), reference.len() - 1);
    debug_assert!(result.iter().all(|(_, v)| &0.0 <= v && v <= &1.0));

    result
}

fn interleave_sample(
    sample: &ContigIntervalMap<EqFloat>,
    bpoints: &HashMap<ArcStr, Vec<PosType>>,
) -> HashMap<ArcStr, Vec<(u32, f64)>> {
    sample
        .inner()
        .par_iter()
        .map(|(k, v)| (ArcStr::from(k.as_str()), v))
        .filter_map(|(chr, lapper)| bpoints.get(&chr).map(|r| (chr, lapper, r)))
        .map(|(chr, lapper, merged)| {
            (
                ArcStr::from(chr.as_str()),
                interleave_single(lapper, merged),
            )
        })
        .collect::<HashMap<_, _>>()
}

pub fn merge_breakpoints(
    intervals: &[ContigIntervalMap<EqFloat>],
    merge_type: MergeType,
) -> Vec<ContigIntervalMap<EqFloat>> {
    let all_breakpoints = gather_breakpoints(intervals);

    let merged_breakpoints = THREAD_POOL.install(move || {
        match merge_type {
            MergeType::Full => preprocess_all_bpoints(all_breakpoints),
            MergeType::Dbscan(eps, mpt, max_dist, agg_method) => {
                let partitioned = partition_by_dist(all_breakpoints, max_dist);
                apply_dbscan(partitioned, eps, mpt, agg_method)
            },
        }
    });

    let result = THREAD_POOL.install(|| {
        let interleaved = intervals
            .par_iter()
            .map(|int| interleave_sample(int, &merged_breakpoints))
            .collect::<Vec<_>>();

        let all_chr = intervals
            .iter()
            .map(|imap| imap.chr_names())
            .flatten()
            .collect::<HashSet<_>>();

        interleaved
            .into_par_iter()
            .map(|hmap| {
                let mut res = hmap
                    .into_par_iter()
                    .map(|(chr, values)| {
                        let mut res_intervals = Vec::with_capacity(values.len());
                        let mut prev_value = 0;
                        for (v, score) in values {
                            res_intervals.push(Interval {
                                start: prev_value,
                                stop:  v,
                                val:   EqFloat::try_from(score).unwrap(),
                            });
                            prev_value = v;
                        }
                        (BsxSmallStr::from(chr.as_str()), Lapper::new(res_intervals))
                    })
                    .collect::<HashMap<_, _>>();
                for chr in all_chr.iter() {
                    if !res.contains_key(chr) {
                        res.insert(chr.to_owned(), Lapper::new(vec![]));
                    }
                }

                res
            })
            .map(ContigIntervalMap::from)
            .collect::<Vec<_>>()
    });

    result
}
