use once_cell::sync::OnceCell;
use std::cmp::Ordering;
use std::sync::Arc;

use crate::tools::dmr::segmentation;
use crate::utils::mann_whitney_u;
use itertools::Itertools;
use serde::{Deserialize, Serialize, Serializer};

#[derive(Clone, Debug)]
pub struct SegmentView<'a> {
    pub(crate) pvalue:    OnceCell<f64>,
    pub(crate) rel_start: usize,
    pub(crate) rel_end:   usize,
    parent:               Arc<&'a SegmentOwned>,
}

impl<'a> SegmentView<'a> {
    pub fn new(
        rel_start: usize,
        rel_end: usize,
        parent: Arc<&'a SegmentOwned>,
    ) -> SegmentView<'a> {
        let pvalue = OnceCell::new();
        SegmentView {
            pvalue,
            rel_start,
            rel_end,
            parent,
        }
    }

    pub fn mds_orig(&self) -> &[f32] {
        &self.parent.mds_orig[self.rel_start..self.rel_end]
    }

    pub fn group_a(&self) -> &[f32] {
        &self.parent.group_a[self.rel_start..self.rel_end]
    }

    pub fn group_b(&self) -> &[f32] {
        &self.parent.group_b[self.rel_start..self.rel_end]
    }

    pub fn start_pos(&self) -> u64 { self.parent.positions[self.rel_start] }

    pub fn end_pos(&self) -> u64 { self.parent.positions[self.rel_end - 1] }

    #[allow(dead_code)]
    pub fn positions(&self) -> &[u64] {
        &self.parent.positions[self.rel_start..self.rel_end]
    }

    pub fn size(&self) -> usize { self.rel_end - self.rel_start }

    pub fn slice(
        &self,
        start: usize,
        end: usize,
    ) -> SegmentView<'a> {
        SegmentView::new(
            self.rel_start + start,
            if self.rel_start + end <= self.rel_end {
                self.rel_start + end
            }
            else {
                self.rel_end
            },
            self.parent.clone(),
        )
    }

    pub fn merge(
        &self,
        other: &Self,
    ) -> SegmentView<'a> {
        debug_assert!(
            other
                .rel_start
                .saturating_sub(self.rel_end)
                < 2
        );
        SegmentView::new(
            self.rel_start,
            other.rel_end,
            self.parent.clone(),
        )
    }

    #[allow(dead_code)]
    pub fn to_owned(&self) -> SegmentOwned {
        SegmentOwned {
            positions: self.positions().to_vec(),
            group_a:   self.group_a().to_vec(),
            group_b:   self.group_b().to_vec(),
            mds_orig:  self.mds_orig().to_vec(),
        }
    }

    pub fn get_pvalue(&self) -> f64 {
        *self.pvalue.get_or_init(|| {
            let (_, prob) = mann_whitney_u(
                self.group_a(),
                self.group_b()
            );
            prob
        })
    }
}

impl PartialEq for SegmentOwned {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.size() == other.size()
            && self.positions.first() == other.positions.first()
            && self.positions.last() == other.positions.last()
    }
}

impl PartialEq for SegmentView<'_> {
    fn eq(
        &self,
        other: &SegmentView,
    ) -> bool {
        if self.parent == other.parent {
            self.rel_start == other.rel_start && self.rel_end == other.rel_end
        }
        else {
            false
        }
    }
}

impl PartialOrd for SegmentView<'_> {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {
        if self.parent == other.parent {
            self.rel_start
                .partial_cmp(&other.rel_start)
        }
        else {
            None
        }
    }
}

impl Eq for SegmentView<'_> {}

impl Ord for SegmentView<'_> {
    fn cmp(
        &self,
        other: &Self,
    ) -> Ordering {
        self.partial_cmp(other)
            .unwrap_or_else(|| panic!("Different parent segments"))
    }
}

#[derive(Debug)]
pub struct SegmentOwned {
    group_a:   Vec<f32>,
    group_b:   Vec<f32>,
    mds_orig:  Vec<f32>,
    positions: Vec<u64>,
}

impl SegmentOwned {
    pub fn new(
        positions: Vec<u64>,
        group_a: Vec<f32>,
        group_b: Vec<f32>,
    ) -> Self {
        let mds_orig = group_a
            .iter()
            .zip(group_b.iter())
            .map(|(a, b)| a - b)
            .collect_vec();
        assert_eq!(positions.len(), mds_orig.len());

        Self {
            group_a,
            group_b,
            mds_orig,
            positions,
        }
    }

    pub fn split_by_dist(
        mut self,
        max_dist: u64,
    ) -> Vec<Self> {
        let mut split_idxs =
            segmentation::arg_split_segment(&self.positions, max_dist);
        split_idxs.reverse();
        if split_idxs.is_empty() {
            vec![self]
        }
        else {
            let mut res =
                split_idxs
                    .into_iter()
                    .fold(Vec::new(), |mut acc, idx| {
                        let positions = self
                            .positions
                            .drain(idx..)
                            .collect_vec();
                        let group_a = self.group_a.drain(idx..).collect_vec();
                        let group_b = self.group_b.drain(idx..).collect_vec();
                        let mds_orig = self.mds_orig.drain(idx..).collect_vec();

                        acc.push(Self {
                            positions,
                            group_a,
                            group_b,
                            mds_orig,
                        });
                        acc
                    });
            res.reverse();
            res
        }
    }

    pub fn concat(
        mut self,
        mut other: SegmentOwned,
    ) -> SegmentOwned {
        self.group_a.append(&mut other.group_a);
        self.group_b.append(&mut other.group_b);
        self.mds_orig
            .append(&mut other.mds_orig);
        self.positions
            .append(&mut other.positions);
        self
    }

    pub fn to_view(&self) -> SegmentView {
        SegmentView::new(0, self.mds_orig.len(), Arc::new(self))
    }

    pub fn size(&self) -> usize { self.positions.len() }
}

pub struct ReaderMetadata {
    pub(crate) blocks_total:  usize,
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
    pub chr:         String,
    pub start:       u32,
    pub end:         u32,
    #[serde(serialize_with = "serialize_scientific")]
    pub p_value:     f64,
    pub meth_left:   f32,
    pub meth_right:  f32,
    pub n_cytosines: usize,
    pub meth_diff:   f32,
    pub meth_mean:   f32,
}

impl DMRegion {
    pub(crate) fn from_segment_view(segment: SegmentView, chr: String) -> Self {
        let a_mean = segment.group_a().iter().sum::<f32>() / segment.size() as f32;
        let b_mean = segment.group_b().iter().sum::<f32>() / segment.size() as f32;

        DMRegion {
            chr,
            start: segment.start_pos() as u32,
            end: segment.end_pos() as u32,
            p_value: segment.get_pvalue(),
            meth_left: a_mean,
            meth_right: b_mean,
            n_cytosines: segment.size(),
            meth_diff: a_mean - b_mean,
            meth_mean: (a_mean + b_mean) / 2.0,
        }
    }

    pub fn meth_diff(&self) -> f32 { self.meth_left - self.meth_right }

    #[allow(dead_code)]
    fn meth_mean(&self) -> f32 { (self.meth_left + self.meth_right) / 2.0 }

    pub fn length(&self) -> u32 { self.end - self.start + 1 }
}

fn serialize_scientific<S>(
    x: &f64,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer, {
    serializer.serialize_str(&format!("{:e}", x))
}
