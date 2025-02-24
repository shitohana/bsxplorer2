use core::sync;
use itertools::Itertools;
use rayon::prelude::*;
use statrs::statistics::Statistics;
use std::cmp::{Ordering, PartialEq};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::ops::Not;
use std::sync::atomic::AtomicUsize;
use std::sync::{Arc, Mutex, RwLock};

mod dmr_model;
mod ks2d;
pub use dmr_model::*;
pub use ks2d::ks2d_2sample;

/// Numerator    = |S(a, b)| - (b - a) * |S(start, end)| / (t - s)
///
/// Denominator  = (b - a) * \[1 - (b - a) / (end - start)\]
///
/// Z = Numerator^2 / Denominator
fn scoring_function(mds_cumsum: &[f64], start: usize, end: usize, a: usize, b: usize) -> f64 {
    let numerator = (mds_cumsum[b] - mds_cumsum[a]).abs()
        - (b - a) as f64 / (end - start) as f64 * (mds_cumsum[end] - mds_cumsum[start]).abs();

    let denominator = (b - a) as f64 * (1f64 - (b - a) as f64 / (end - start) as f64);

    numerator.powi(2) / denominator
}

/// ### Returns
///
/// Tuple (a, b, score), where a, b - interval [a, b) and score - maximum score value
fn maximize_scoring(mds_cumsum: &[f64], start: usize, end: usize) -> (usize, usize, f64) {
    (start..end)
        .combinations(2)
        .par_bridge()
        .map(|ab_vec| {
            (
                ab_vec[0],
                ab_vec[1],
                scoring_function(mds_cumsum, start, end, ab_vec[0], ab_vec[1]),
            )
        })
        .max_by(|x, y| x.2.total_cmp(&y.2))
        .unwrap_or((0, 0, 0f64))
}

#[derive(Clone, Debug)]
pub struct PreSegmentView<'a> {
    pvalue: Option<f64>,
    rel_start: usize,
    rel_end: usize,
    parent: Arc<&'a PreSegmentOwned>,
}

impl<'a> PreSegmentView<'a> {
    pub fn init_pvalue(&mut self) {
        if self.pvalue.is_none() {
            let (_, prob) = self.ks_test();
            self.pvalue = Some(prob)
        }
    }

    pub fn mds_cumsum(&self) -> &[f64] {
        &self.parent.mds_cumsum[self.rel_start..self.rel_end]
    }
    pub fn mds_orig(&self) -> &[f64] {
        &self.parent.mds_orig[self.rel_start..self.rel_end]
    }
    pub fn group_a(&self) -> &[f64] {
        &self.parent.group_a[self.rel_start..self.rel_end]
    }
    pub fn group_b(&self) -> &[f64] {
        &self.parent.group_b[self.rel_start..self.rel_end]
    }
    pub fn mean_diff(&self) -> f64 {
        self.group_a().mean() - self.group_b().mean()
    }
    pub fn start_pos(&self) -> u64 {
        self.parent.positions[self.rel_start]
    }
    pub fn end_pos(&self) -> u64 {
        self.parent.positions[self.rel_end - 1]
    }
    pub fn positions(&self) -> &[u64] {
        &self.parent.positions[self.rel_start..self.rel_end]
    }
    pub fn size(&self) -> usize {
        self.rel_end - self.rel_start
    }
    pub fn iter_points_a(&self) -> impl Iterator<Item = (&u64, &f64)> {
        self.positions().iter().zip(self.group_a().iter())
    }
    pub fn iter_points_b(&self) -> impl Iterator<Item = (&u64, &f64)> {
        self.positions().iter().zip(self.group_b().iter())
    }

    pub fn slice(self, start: usize, end: usize) -> PreSegmentView<'a> {
        PreSegmentView {
            pvalue: None,
            rel_start: self.rel_start + start,
            rel_end: if self.rel_start + end <= self.rel_end {
                self.rel_start + end
            } else {
                self.rel_end
            },
            parent: Arc::clone(&self.parent),
        }
    }

    pub fn ks_test(&self) -> (f64, f64) {
        ks2d::ks2d_2sample(
            self.positions(),
            self.group_a(),
            self.positions(),
            self.group_b(),
        )
    }
    pub fn calc_trend_abs(&self) -> f64 {
        if self.mds_cumsum().is_empty() {
            return 0.0;
        }
        (self.mds_cumsum()[self.size() - 1] - self.mds_cumsum()[0]).abs() / (self.size() as f64)
    }
    pub fn has_significant_trend(&self, trend_threshold: f64) -> bool {
        let trend = self.calc_trend_abs();
        trend > trend_threshold
    }

    pub fn no_valley(&self, tmin: usize, valley_threshold: f64) -> bool {
        // If there are fewer CpGs than the minimum, we say “no valley”
        if self.size() < tmin {
            return true;
        }
        // Compute the overall average absolute MDS in the segment.
        let overall_avg = self.calc_trend_abs();
        // Slide a window of length tmin over the segment.
        (0..(self.size() - tmin))
            .par_bridge()
            .map(|start_idx| {
                (self.mds_cumsum()[start_idx + tmin] - self.mds_cumsum()[start_idx]).abs()
                    / tmin as f64
            })
            .map(|window_avg| window_avg < valley_threshold * overall_avg)
            .any(|x| x)
            .not()
    }

    pub fn to_owned(&self) -> PreSegmentOwned {
        PreSegmentOwned {
            positions: self.positions().to_vec(),
            group_a: self.group_a().to_vec(),
            group_b: self.group_b().to_vec(),
            mds_cumsum: self.mds_cumsum().to_vec(),
            mds_orig: self.mds_orig().to_vec(),
        }
    }
}

impl PartialEq for PreSegmentOwned {
    fn eq(&self, other: &Self) -> bool {
        self.size() == other.size()
            && self.positions.first() == other.positions.first()
            && self.positions.last() == other.positions.last()
    }
}

impl PartialEq for PreSegmentView<'_> {
    fn eq(&self, other: &PreSegmentView) -> bool {
        if self.parent == other.parent {
            self.rel_start == other.rel_start && self.rel_end == other.rel_end
        } else {
            false
        }
    }
}

impl PartialOrd for PreSegmentView<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.parent == other.parent {
            self.rel_start.partial_cmp(&other.rel_start)
        } else {
            None
        }
    }
}

impl Eq for PreSegmentView<'_> {}

impl Ord for PreSegmentView<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other)
            .unwrap_or_else(|| panic!("Different parent segments"))
    }
}

#[derive(Debug)]
pub struct PreSegmentOwned {
    group_a: Vec<f64>,
    group_b: Vec<f64>,
    mds_cumsum: Vec<f64>,
    mds_orig: Vec<f64>,
    positions: Vec<u64>,
}

impl PreSegmentOwned {
    pub fn new(positions: Vec<u64>, group_a: Vec<f64>, group_b: Vec<f64>) -> Self {
        let mds_orig = group_a
            .iter()
            .zip(group_b.iter())
            .map(|(a, b)| a - b)
            .collect_vec();
        assert_eq!(positions.len(), mds_orig.len());
        let mut cum = Vec::with_capacity(mds_orig.len());
        for &val in mds_orig.iter() {
            cum.push(cum.last().unwrap_or(&0.0) + val);
        }

        Self {
            group_a,
            group_b,
            mds_cumsum: cum,
            mds_orig,
            positions,
        }
    }

    pub fn split_by_dist(mut self, max_dist: u64) -> Vec<Self> {
        let mut split_idxs = arg_split_segment(&self.positions, max_dist);
        split_idxs.reverse();
        if split_idxs.is_empty() {
            vec![self]
        } else {
            let mut res = split_idxs.into_iter().fold(Vec::new(), |mut acc, idx| {
                let positions = self.positions.drain(idx..).collect_vec();
                let group_a = self.group_a.drain(idx..).collect_vec();
                let group_b = self.group_b.drain(idx..).collect_vec();
                let mds_orig = self.mds_orig.drain(idx..).collect_vec();
                let mds_cumsum = self.mds_cumsum.drain(idx..).collect_vec();

                acc.push(Self {
                    positions,
                    group_a,
                    group_b,
                    mds_cumsum,
                    mds_orig,
                });
                acc
            });
            res.reverse();
            res
        }
    }

    pub fn concat(mut self, mut other: PreSegmentOwned) -> PreSegmentOwned {
        self.group_a.append(&mut other.group_a);
        self.group_b.append(&mut other.group_b);
        self.mds_orig.append(&mut other.mds_orig);
        self.positions.append(&mut other.positions);
        self.mds_cumsum.append(&mut other.mds_cumsum);
        self
    }

    pub fn to_view(&self) -> PreSegmentView {
        PreSegmentView {
            pvalue: None,
            rel_start: 0,
            rel_end: self.mds_orig.len(),
            parent: Arc::new(self),
        }
    }

    pub fn slice(&self, start: usize, end: usize) -> PreSegmentView {
        PreSegmentView {
            pvalue: None,
            rel_start: start,
            rel_end: end,
            parent: Arc::new(self),
        }
    }

    pub fn length(&self) -> u64 {
        self.positions
            .first()
            .cloned()
            .unwrap_or(0)
            .saturating_sub(self.positions.last().cloned().unwrap_or(0))
    }

    pub fn size(&self) -> usize {
        self.positions.len()
    }

    pub fn group_a(&self) -> &Vec<f64> {
        &self.group_a
    }

    pub fn group_b(&self) -> &Vec<f64> {
        &self.group_b
    }

    pub fn mds_cumsum(&self) -> &Vec<f64> {
        &self.mds_cumsum
    }

    pub fn mds_orig(&self) -> &Vec<f64> {
        &self.mds_orig
    }

    pub fn positions(&self) -> &Vec<u64> {
        &self.positions
    }
}

#[derive(Debug)]
pub struct PotentialDMR {
    segment: PreSegmentOwned,
    p_value: f64,
}

impl PotentialDMR {
    pub fn new(segment: PreSegmentOwned, p_value: f64) -> Self {
        Self { segment, p_value }
    }

    pub fn from_view(view: &PreSegmentView) -> Self {
        let mut view = view.clone();
        view.init_pvalue();
        Self {
            p_value: view.pvalue.unwrap(),
            segment: view.to_owned(),
        }
    }

    pub fn segment(&self) -> &PreSegmentOwned {
        &self.segment
    }

    pub fn p_value(&self) -> f64 {
        self.p_value
    }
}

impl Eq for PotentialDMR {}

impl PartialEq<Self> for PotentialDMR {
    fn eq(&self, other: &Self) -> bool {
        self.segment.positions.first() == other.segment.positions.first()
            && self.segment.positions.last() == other.segment.positions.last()
    }
}

impl PartialOrd<Self> for PotentialDMR {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PotentialDMR {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.segment.positions.last().unwrap() < other.segment.positions.first().unwrap() {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    }
}

// Currently this just calculates mean
// Maybe change it to random number from distribution
pub fn estimate_missing(non_missing_density: &[f64]) -> f64 {
    // Todo: mean and variance can be calculated by one pass
    let mean = non_missing_density.mean();
    let variance = non_missing_density.variance();

    let c = mean * (1f64 - mean) / variance - 1f64;
    let alpha = c * mean;
    let beta = (1f64 - mean) * c;

    alpha / (alpha + beta)
}

pub fn run_segmentation(
    segment: PreSegmentView,
    min_cpg: usize,
    trend_threshold: f64,
    valley: f64,
) -> Vec<Option<PreSegmentView>> {
    use crossbeam::channel;

    // Each segment gets its ID (just by count), to keep track of all previous segments.
    // Value represents [PreSegmentView] and boolean flag, indicating whether it is a potential
    // DMR or not.
    let seg_mapping = Arc::new(RwLock::new(BTreeMap::<usize, (PreSegmentView, bool)>::new()));
    // Each segment can have only one parent. There is no need to find children of searching
    // for children nodes. Instead, we need only parent nodes, that is why we just use BTreeMap.
    let parent_graph = Arc::new(RwLock::new(BTreeMap::<usize, Option<usize>>::new()));
    // Insert root node
    let root_segment = segment.clone();
    let root_id = update_mappings(
        root_segment,
        false,
        None,
        seg_mapping.clone(),
        parent_graph.clone(),
    );

    // Sender and receiver channels are used to avoid recursion, both  sender and receiver take
    // segment ID from seg_mapping, which means, that before sending new segment it needs first
    // to be added into the mapping.
    let (sender, reciever): (channel::Sender<usize>, channel::Receiver<usize>) =
        channel::bounded(rayon::max_num_threads());
    // Initial segmentation step
    let (initial_a, initial_b, _) = maximize_scoring(segment.mds_cumsum(), 0, segment.size() - 1);

    let presegments = [
        segment.clone().slice(0, initial_a),
        segment.clone().slice(initial_a, initial_b),
        segment.clone().slice(initial_b, segment.size()),
    ]
    .into_iter()
    .filter(|s| s.size() > 0)
    .collect_vec();

    for ps in presegments {
        let seg_id = update_mappings(
            ps.clone(),
            false,
            Some(root_id),
            seg_mapping.clone(),
            parent_graph.clone(),
        );
        sender.send(seg_id).expect("Failed to send segment");
    }
    let sender_arc = Arc::new(sender);

    let running_tasks = Arc::new(AtomicUsize::new(0));
    // Segmentation loop
    while !reciever.is_empty() {
        let seg_id = reciever.recv().unwrap();
        let seg_mapping_clone = Arc::clone(&seg_mapping);
        let parent_graph_clone = Arc::clone(&parent_graph);
        let sender_arc_clone = Arc::clone(&sender_arc);
        let running_tasks_clone = running_tasks.clone();

        rayon::scope(|scope| {
            scope.spawn(move |_| {
                running_tasks_clone.fetch_add(1, sync::atomic::Ordering::SeqCst);
                // Get segment from mapping
                let (mut cur_seg, _) = seg_mapping_clone
                    .read()
                    .expect("Segment mapping read error")
                    .get(&seg_id)
                    .cloned()
                    .unwrap_or_else(|| panic!("Failed to get segment {:?}", seg_id));

                let mut do_segmentation = true;
                // Check is there enough cytosines
                if cur_seg.size() <= min_cpg {
                    cur_seg.init_pvalue();
                    seg_mapping_clone
                        .write()
                        .unwrap()
                        .insert(seg_id, (cur_seg.clone(), true));
                    do_segmentation = false;
                } else if !cur_seg.has_significant_trend(trend_threshold)
                    && cur_seg.no_valley(min_cpg, valley)
                {
                    // Assign 2D KS statistics
                    cur_seg.init_pvalue();

                    let parent_id = parent_graph_clone
                        .read()
                        .expect("Parent mapping read error")
                        .get(&seg_id)
                        .cloned()
                        .unwrap_or_else(|| panic!("Failed to get segment {:?}", seg_id))
                        .unwrap_or_else(|| {
                            unreachable!("Root node should not be addressed in channel")
                        });

                    let (parent_seg, _) = seg_mapping_clone
                        .read()
                        .expect("Segment mapping read error")
                        .get(&parent_id)
                        .unwrap_or_else(|| unreachable!("Parent node should exist in channel"))
                        .clone();
                    let parent_pvalue = parent_seg.pvalue.clone();

                    // If parent ID exists - mark as DMR
                    let mark_dmr = if let Some(parent_pvalue) = parent_pvalue {
                        cur_seg.pvalue.unwrap() < parent_pvalue
                    } else {
                        false
                    };

                    // Update mapping with p-value, and possible DMR
                    seg_mapping_clone
                        .write()
                        .unwrap()
                        .insert(seg_id, (cur_seg.clone(), mark_dmr));
                    // If is DMR - return
                    if mark_dmr {
                        do_segmentation = false;
                    }
                } else {
                    // TODO maybe remove. Test
                    cur_seg.init_pvalue();
                    seg_mapping_clone
                        .write()
                        .unwrap()
                        .insert(seg_id, (cur_seg.clone(), false));
                }

                if do_segmentation {
                    // Segment again
                    let (a, b, _) = maximize_scoring(cur_seg.mds_cumsum(), 0, cur_seg.size() - 1);
                    let presegments = [
                        cur_seg.clone().slice(0, a),
                        cur_seg.clone().slice(a, b + 1),
                        cur_seg.clone().slice(b + 1, cur_seg.size()),
                    ]
                    .into_iter()
                    .filter(|s| s.size() > 0)
                    .collect_vec();
                    // If no segmentation occured (i.e. maximum score is the full window)
                    // stop recursion
                    if presegments.len() == 1 {
                        cur_seg.init_pvalue();
                        seg_mapping_clone
                            .write()
                            .unwrap()
                            .insert(seg_id, (cur_seg.clone(), true));
                        do_segmentation = false;
                    }
                    // Otherwise continue
                    for ps in presegments {
                        let seg_id = update_mappings(
                            ps.clone(),
                            false,
                            Some(seg_id),
                            seg_mapping_clone.clone(),
                            parent_graph_clone.clone(),
                        );
                        sender_arc_clone
                            .send(seg_id)
                            .expect("Failed to send segment");
                    }
                }
                running_tasks_clone.fetch_sub(1, sync::atomic::Ordering::SeqCst);
            });
        });
    }

    let c: HashSet<usize> = HashSet::from_iter(parent_graph.read().unwrap().keys().copied());
    let most_prob_res = c
        .iter()
        .map(|id| {
            let seg_mapping_handle = seg_mapping.read().expect("Segment mapping read error");
            let (segment, is_dmr) = seg_mapping_handle.get(id).unwrap();
            (id, segment.pvalue, *is_dmr)
        })
        // .filter(|(_, pvalue, is_dmr)| *is_dmr && pvalue.is_some())
        .filter(|(_, pvalue, is_dmr)| pvalue.is_some())
        .map(|(idx, pvalue, dmr)| (idx, pvalue.unwrap(), dmr))
        .filter(|(_, pvalue, _)| !pvalue.is_nan())
        .min_by(|(_, pvalue, _), (_, pvalue2, _)| pvalue.partial_cmp(pvalue2).unwrap());

    let mut result: Vec<Option<PreSegmentView>> = Vec::new();
    if let Some((most_prob_id, _, _)) = most_prob_res {
        let seg_mapping_handle = seg_mapping.read().expect("Segment mapping read error");
        let (dmr_segment, is_dmr) = seg_mapping_handle.get(&most_prob_id).unwrap();
        let (a, b) = (dmr_segment.rel_start, dmr_segment.rel_end);

        // Before DMR
        if a != segment.rel_start {
            let mut new_seg = segment.clone();
            new_seg.rel_end = a;
            new_seg.pvalue = None;
            result.push(Some(new_seg));
        } else {
            result.push(None);
        }
        // DMR
        if a != b {
            result.push(Some(dmr_segment.clone()));
        } else {
            result.push(None);
        }

        if b + 1 < segment.rel_end {
            let mut new_seg = segment.clone();
            new_seg.rel_start = b + 1;
            new_seg.pvalue = None;
            result.push(Some(new_seg));
        } else {
            result.push(None);
        }
    } else {
        result.push(Some(segment));
    }
    result
}

fn recurse_segmentation(
    segment: PreSegmentView,
    min_cpg: usize,
    trend_threshold: f64,
    valley: f64,
    diff_threshold: f64,
) -> BTreeSet<PreSegmentView> {
    let mut pending_segments = vec![segment];
    let mut result_segments = BTreeSet::new();

    loop {
        if pending_segments.len() != 0 {
            let mut segmentation_result = run_segmentation(
                pending_segments.pop().unwrap(),
                min_cpg,
                trend_threshold,
                valley,
            );

            if segmentation_result.len() == 3 {
                if let Some(segment) = segmentation_result.pop().unwrap() {
                    pending_segments.push(segment);
                }

                // DMR
                if let Some(dmr_seg) = segmentation_result.pop().unwrap() {
                    result_segments.insert(dmr_seg);
                }

                if let Some(segment) = segmentation_result.pop().unwrap() {
                    pending_segments.push(segment);
                }
            } else if segmentation_result.len() == 1 {
                result_segments.insert(segmentation_result.pop().unwrap().unwrap());
            } else {
                unreachable!();
            }
        } else {
            break;
        }
    }

    result_segments
}

fn update_mappings<'a>(
    node: PreSegmentView<'a>,
    potential_dmr: bool,
    parent: Option<usize>,
    seg_mapping: Arc<RwLock<BTreeMap<usize, (PreSegmentView<'a>, bool)>>>,
    parent_graph: Arc<RwLock<BTreeMap<usize, Option<usize>>>>,
) -> usize {
    let new_id = parent_graph
        .read()
        .expect("Failed to get lock on parent graph")
        .len();

    parent_graph
        .write()
        .expect("Failed to get lock on parent graph")
        .insert(new_id, parent);
    seg_mapping
        .write()
        .expect("Failed to get lock on segment mapping")
        .insert(new_id, (node, potential_dmr));
    new_id
}

fn pre_segment<'a>(segment: &'a PreSegmentView<'a>) -> Vec<PreSegmentView<'a>> {
    let (initial_a, initial_b, _) = maximize_scoring(segment.mds_cumsum(), 0, segment.size() - 1);

    let pre_segments = [
        segment.clone().slice(0, initial_a),
        segment.clone().slice(initial_a, initial_b),
        segment.clone().slice(initial_b, segment.size()),
    ]
    .into_iter()
    .filter(|s| s.size() > 0)
    .collect_vec();

    pre_segments
}

pub fn segment_recursive(
    cur_segment: PreSegmentView,
    min_cpg: usize,
    trend_threshold: f64,
    valley: f64,
    potential_dmrs: Arc<Mutex<BTreeSet<PotentialDMR>>>,
) {
    let (initial_a, initial_b, _) =
        maximize_scoring(cur_segment.mds_cumsum(), 0, cur_segment.size() - 1);

    let pre_segments = vec![
        cur_segment.clone().slice(0, initial_a),
        cur_segment.clone().slice(initial_a, initial_b),
        cur_segment.clone().slice(initial_b, cur_segment.size()),
    ];

    pre_segments.par_iter().for_each(|seg| {
        if seg.size() == 0 {
            return;
        }
        if seg.size() <= min_cpg {
            let mut result = seg.clone();
            result.init_pvalue();

            let dmr = PotentialDMR {
                p_value: result.pvalue.unwrap(),
                segment: result.to_owned(),
            };

            potential_dmrs.clone().lock().unwrap().insert(dmr);
        } else if !seg.has_significant_trend(trend_threshold) {
            if seg.no_valley(min_cpg, valley) {
                let mut parent_clone = cur_segment.clone();
                parent_clone.init_pvalue();
                let mut self_clone = seg.clone();
                self_clone.init_pvalue();

                if self_clone.pvalue.unwrap() > parent_clone.pvalue.unwrap() {
                    let dmr = PotentialDMR {
                        p_value: self_clone.pvalue.unwrap(),
                        segment: self_clone.to_owned(),
                    };
                    potential_dmrs.clone().lock().unwrap().insert(dmr);
                } else {
                    segment_recursive(
                        seg.clone(),
                        min_cpg,
                        trend_threshold,
                        valley,
                        potential_dmrs.clone(),
                    )
                }
            } else {
                segment_recursive(
                    seg.clone(),
                    min_cpg,
                    trend_threshold,
                    valley,
                    potential_dmrs.clone(),
                )
            }
        } else {
            segment_recursive(
                seg.clone(),
                min_cpg,
                trend_threshold,
                valley,
                potential_dmrs.clone(),
            )
        }
    });
}

pub fn arg_split_segment(positions: &[u64], max_dist: u64) -> Vec<usize> {
    (1..positions.len())
        .into_iter()
        .filter_map(|i| {
            if positions[i] - positions[i - 1] >= max_dist {
                Some(i)
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use polars::prelude::*;
    use rand::prelude::*;
    use rand::thread_rng;
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::Write;
    use std::ops::Div;

    #[test]
    fn estimate_missing_test() {
        let non_missing_density = vec![0.3, 0.5];
        assert_eq!(estimate_missing(&non_missing_density), 0.4);
    }

    #[test]
    fn pearson_r_test() {
        let x = vec![1, 2, 3, 4, 5, 6];
        let y = vec![6, 7, 8, 9, 10, 11];
        assert_eq!(pearson_r(&x, &y), 1f64);
    }

    #[test]
    fn ks_test() {
        let x1 = vec![1, 2, 3, 4, 5, 6];
        let y1 = vec![1, 4, 9, 16, 25, 36];

        let x2 = vec![1, 2, 3, 4, 5, 6];
        let y2 = vec![1, 4, 9, 16, 25, 36];
        let (prob, _) = ks2d_2sample(&x1, &y1, &x2, &y2);
        assert_approx_eq!(prob, 0f64);

        let x1 = vec![1, 2, 3, 4, 5, 6];
        let y1 = statrs::distribution::Normal::standard()
            .sample_iter(thread_rng())
            .take(x1.len())
            .collect_vec();

        let x2 = vec![1, 2, 3, 4, 5, 6];
        let y2 = statrs::distribution::Beta::new(1f64, 3f64)
            .unwrap()
            .sample_iter(thread_rng())
            .take(x2.len())
            .collect_vec();
        let (prob, d) = ks2d_2sample(&x1, &y1, &x2, &y2);
        println!("Normal vs Beta(1, 3): prob {prob:.3e}, d: {d:.3e}")
    }

    use crate::tools::metilene::ks2d::{ks2d_2sample, pearson_r};
    use assert_approx_eq::assert_approx_eq;
    use itertools::Itertools;
    use polars::prelude::{CsvParseOptions, SerReader};
    use std::sync::{Arc, Mutex};

    #[test]
    fn test_segment_recursive() {
        // Create a synthetic pre-segment with 8 CpGs having a very strong difference,
        // so that the MDS (group_a - group_b) is high and relatively constant.
        let positions = vec![90, 91, 92, 100, 101, 102, 103, 104, 105, 106, 107];
        let group_a = vec![0.3, 0.3, 0.3, 0.4, 0.9, 0.9, 0.9, 0.9, 0.3, 0.3, 0.3];
        let group_b = vec![0.3; 11];
        let pre_seg = PreSegmentOwned::new(positions, group_a, group_b);
        let view = pre_seg.to_view();

        // Set parameters:
        // - min_cpg: at least 3 CpGs required
        // - trend_threshold: e.g., 0.5 (a strong trend is expected when the cumulative change is high)
        // - valley: 0.8 (if any window’s average MDS is less than 80% of overall average, that indicates a valley)
        let min_cpg = 10;
        let trend_threshold = 0.5;
        let valley = 0.3;

        let potential_dmrs = Arc::new(Mutex::new(BTreeSet::<PotentialDMR>::new()));
        segment_recursive(
            view,
            min_cpg,
            trend_threshold,
            valley,
            potential_dmrs.clone(),
        );

        let results = potential_dmrs.lock().unwrap();
        // We expect that at least one potential DMR is produced.
        assert!(
            !results.is_empty(),
            "No potential DMRs were produced by segment_recursive"
        );
        // For debugging, print the results.
        for dmr in results.iter() {
            println!(
                "Potential DMR: p_value: {}, segment size: {:?}, {:?}",
                dmr.p_value, dmr.segment.mds_orig, dmr.segment.positions
            );
        }
    }
}
