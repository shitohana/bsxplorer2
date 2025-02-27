use crate::tools::dmr::segmentation;
use crate::utils::ks2d_2sample;
use itertools::Itertools;
use serde::{Deserialize, Serialize, Serializer};
use statrs::statistics::Statistics;
use std::cmp::Ordering;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct SegmentView<'a> {
    pub(crate) pvalue: Option<f64>,
    pub(crate) rel_start: usize,
    pub(crate) rel_end: usize,
    parent: Arc<&'a SegmentOwned>,
    mds_cumsum: Vec<f64>,
}

impl<'a> SegmentView<'a> {
    pub fn new(
        pvalue: Option<f64>,
        rel_start: usize,
        rel_end: usize,
        parent: Arc<&'a SegmentOwned>,
    ) -> SegmentView<'a> {
        let iter = &parent.mds_cumsum[rel_start..rel_end];
        let first = iter.first().cloned().unwrap_or(0.0);
        let mds_cumsum = iter.into_iter().map(|x| *x - first).collect::<Vec<_>>();

        SegmentView {
            pvalue,
            rel_start,
            rel_end,
            mds_cumsum,
            parent,
        }
    }

    pub fn init_pvalue(&mut self) {
        if self.pvalue.is_none() {
            let (_, prob) = self.ks_test();
            self.pvalue = Some(prob)
        }
    }

    pub fn mds_cumsum(&self) -> &[f64] {
        &self.mds_cumsum
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

    pub fn slice(&self, start: usize, end: usize) -> SegmentView<'a> {
        SegmentView::new(
            None,
            self.rel_start + start,
            if self.rel_start + end <= self.rel_end {
                self.rel_start + end
            } else {
                self.rel_end
            },
            self.parent.clone(),
        )
    }

    pub fn merge(&self, other: &Self) -> SegmentView<'a> {
        debug_assert!(other.rel_start.saturating_sub(self.rel_end) < 2);
        SegmentView::new(None, self.rel_start, other.rel_end, self.parent.clone())
    }

    pub fn ks_test(&self) -> (f64, f64) {
        ks2d_2sample(
            self.positions(),
            self.group_a(),
            self.positions(),
            self.group_b(),
        )
    }

    pub fn to_owned(&self) -> SegmentOwned {
        SegmentOwned {
            positions: self.positions().to_vec(),
            group_a: self.group_a().to_vec(),
            group_b: self.group_b().to_vec(),
            mds_cumsum: self.mds_cumsum().to_vec(),
            mds_orig: self.mds_orig().to_vec(),
        }
    }
}

impl PartialEq for SegmentOwned {
    fn eq(&self, other: &Self) -> bool {
        self.size() == other.size()
            && self.positions.first() == other.positions.first()
            && self.positions.last() == other.positions.last()
    }
}

impl PartialEq for SegmentView<'_> {
    fn eq(&self, other: &SegmentView) -> bool {
        if self.parent == other.parent {
            self.rel_start == other.rel_start && self.rel_end == other.rel_end
        } else {
            false
        }
    }
}

impl PartialOrd for SegmentView<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.parent == other.parent {
            self.rel_start.partial_cmp(&other.rel_start)
        } else {
            None
        }
    }
}

impl Eq for SegmentView<'_> {}

impl Ord for SegmentView<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other)
            .unwrap_or_else(|| panic!("Different parent segments"))
    }
}

#[derive(Debug)]
pub struct SegmentOwned {
    group_a: Vec<f64>,
    group_b: Vec<f64>,
    mds_cumsum: Vec<f64>,
    mds_orig: Vec<f64>,
    positions: Vec<u64>,
}

impl SegmentOwned {
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
        let mut split_idxs = segmentation::arg_split_segment(&self.positions, max_dist);
        split_idxs.reverse();
        if split_idxs.is_empty() {
            vec![self]
        } else {
            let mut res = split_idxs.into_iter().fold(Vec::new(), |mut acc, idx| {
                // TODO not sure here
                let first_cumsum = self.mds_cumsum.get(idx).cloned().unwrap_or(0.0);
                let positions = self.positions.drain(idx..).collect_vec();
                let group_a = self.group_a.drain(idx..).collect_vec();
                let group_b = self.group_b.drain(idx..).collect_vec();
                let mds_orig = self.mds_orig.drain(idx..).collect_vec();
                let mds_cumsum = self
                    .mds_cumsum
                    .drain(idx..)
                    .map(|v| v - first_cumsum)
                    .collect_vec();

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

    pub fn concat(mut self, mut other: SegmentOwned) -> SegmentOwned {
        self.group_a.append(&mut other.group_a);
        self.group_b.append(&mut other.group_b);
        self.mds_orig.append(&mut other.mds_orig);
        self.positions.append(&mut other.positions);
        self.mds_cumsum.append(&mut other.mds_cumsum);
        self
    }

    pub fn to_view(&self) -> SegmentView {
        SegmentView::new(None, 0, self.mds_orig.len(), Arc::new(self))
    }

    pub fn slice(&self, start: usize, end: usize) -> SegmentView {
        SegmentView::new(None, start, end, Arc::new(self))
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

pub struct ReaderMetadata {
    pub(crate) blocks_total: usize,
    pub(crate) current_block: usize,
}

impl ReaderMetadata {
    pub(crate) fn new(blocks_total: usize) -> Self {
        Self {
            blocks_total,
            current_block: 0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DMRegion {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    #[serde(serialize_with = "serialize_scientific")]
    pub p_value: f64,
    pub meth_left: f64,
    pub meth_right: f64,
    pub n_cytosines: usize,
    pub meth_diff: f64,
    pub meth_mean: f64,
}

impl DMRegion {
    pub(crate) fn new(
        chr: Arc<String>,
        start: u32,
        end: u32,
        p_value: f64,
        meth_left: f64,
        meth_right: f64,
        n_cytosines: usize,
    ) -> Self {
        DMRegion {
            chr: chr.to_string(),
            start,
            end,
            p_value,
            meth_left,
            meth_right,
            n_cytosines,
            meth_diff: meth_left - meth_right,
            meth_mean: (meth_left + meth_right) / 2.0,
        }
    }

    pub fn meth_diff(&self) -> f64 {
        self.meth_left - self.meth_right
    }

    fn meth_mean(&self) -> f64 {
        (self.meth_left + self.meth_right) / 2.0
    }

    pub fn length(&self) -> u32 {
        self.end - self.start + 1
    }
}

fn serialize_scientific<S>(x: &f64, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(&format!("{:e}", x))
}
