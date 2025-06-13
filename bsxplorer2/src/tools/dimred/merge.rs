#![allow(unused)]
use std::collections::BTreeSet;
use std::error::Error;
use std::io::Read;

use anyhow::bail;
use arcstr::ArcStr;
use bio::bio_types::annot::pos::Pos;
use hashbrown::HashMap;
use itertools::{
    izip,
    Itertools,
};
use rayon::iter::IndexedParallelIterator;
use rayon::prelude::*;
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

pub fn read_segments<R: Read>(handle: R) -> anyhow::Result<ContigIntervalMap<EqFloat>> {
    bio::io::bed::Reader::new(handle)
        .records()
        .into_iter()
        .map(|record_res| -> anyhow::Result<(Contig, EqFloat)> {
            let record = record_res?;
            let score = record
                .score()
                .ok_or(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "File should contain a valid score column",
                ))?
                .to_string()
                .parse::<EqFloat>()?;
            let contig = Contig::from(record);
            Ok((contig, score))
        })
        .collect::<anyhow::Result<ContigIntervalMap<EqFloat>>>()
}

pub enum MergeType {
    Full,
    /// DBSCAN algorithm with parameters (eps, mpt, max_dist)
    Dbscan(f64, usize, PosType),
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
        .map(|(k, v)| (ArcStr::from(k.as_str()), v.concat().into_iter().collect_vec()))
        .map(|(k, mut v)| {v.par_sort_unstable(); (k, v)})
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
    agg_fn: &(dyn Fn(&[f32]) -> f64 + Sync + Send)
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
        .map(|indices| {
            indices.iter().map(|i| pos[*i] as DensityType).collect_vec()
        })
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

fn merge_signal(new: Vec<PosType>, reference: Vec<PosType>, reference_values: Vec<EqFloat>) {
    let zipped_values = vec![vec![-1isize; new.len()], (0isize..reference_values.len() as isize).collect_vec()].concat();
    let zipped = vec![new, reference].concat();

    let sorted_indices = {
        let mut indices = (0..zipped.len())
            .collect_vec();
        indices.par_sort_unstable_by_key(|i| zipped[*i]);
        indices
    };
}

pub fn merge_breakpoints(
    intervals: &[ContigIntervalMap<EqFloat>],
    merge_type: MergeType,
) -> anyhow::Result<ContigIntervalMap<EqFloat>> {
    THREAD_POOL.install(|| -> anyhow::Result<()> {
        let gathered = gather_breakpoints(intervals);
        println!("Before {}", gathered.values().map(Vec::len).sum::<usize>());
        println!("Mean {}", gathered.values().map(Vec::len).sum::<usize>() as f64 / intervals.len() as f64);
        let partitioned = partition_by_dist(gathered, 1000);
        let res = apply_dbscan(partitioned, 20, 3, AggMethod::Median);
        println!("After {:?}", res.values().find_or_first(|_| true).unwrap());
        Ok(())
    })?;
    // let merged = match merge_type {
    //     MergeType::Full => gathered,
    //     MergeType::Dbscan(eps, min_points) => {
    //         let db = dbscan::Model::new(eps, min_points);
    //         gathered.into_iter()
    //             .map(|(k, v)| {
    //                 let labels = db.clone().run(&v);
    //                 let res = v.into_iter()
    //                     .zip(labels.into_iter())
    //                     .filter(|(_, l)| !matches!(l,
    // dbscan::Classification::Noise))                     .map(|(v, _)| v)
    //                     .collect_vec();
    //                 (k, res)
    //             })
    //             .collect()
    //     },
    // };

    // let merged = merged.into_iter()
    //     .map(|(k, v)| (k, v.into_iter().collect::<BTreeSet<_>>()))
    //     .collect::<HashMap<_, _>>();

    // let val = EqFloat::new(0.0).unwrap();
    // let imap = merged.into_iter()
    //     .map(|(k, v)| {
    //        let v_len = v.len();
    //        izip!(v.iter().take(v_len - 1), v.iter().skip(1))
    //            .map(|(start, end)| Contig::new(k.clone(), *start, *end,
    // Strand::None))            .collect_vec()
    //     })
    //     .concat()
    //     .into_iter()
    //     .map(|c| (c, val))
    //     .collect::<ContigIntervalMap<EqFloat>>();

    // Ok(imap)
    Ok(Default::default())
}
