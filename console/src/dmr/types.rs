use std::cmp::Ordering;
use std::sync::Arc;

use bsxplorer2::data_structs::typedef::{
    DensityType,
    PosType,
};
use bsxplorer2::prelude::{
    BsxBatch,
    Contig,
    Strand,
};
use bsxplorer2::utils::mann_whitney_u;
use itertools::Itertools;
use once_cell::sync::OnceCell;
use serde::{
    Deserialize,
    Serialize,
    Serializer,
};

use super::segmentation;

#[derive(Clone, Debug)]
pub struct SegmentView<'a> {
    pub(crate) pvalue:    OnceCell<f64>,
    pub(crate) rel_start: usize,
    pub(crate) rel_end:   usize,
    #[allow(clippy::redundant_allocation)]
    parent:               Arc<&'a SegmentOwned>,
}

#[allow(clippy::redundant_allocation)]
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

    pub fn mds_orig(&self) -> &[DensityType] {
        &self.parent.mds_orig[self.rel_start..self.rel_end]
    }

    pub fn group_a(&self) -> &[DensityType] {
        &self.parent.group_a[self.rel_start..self.rel_end]
    }

    pub fn group_b(&self) -> &[DensityType] {
        &self.parent.group_b[self.rel_start..self.rel_end]
    }

    pub fn start_pos(&self) -> PosType {
        self.parent.positions[self.rel_start]
    }

    pub fn end_pos(&self) -> PosType {
        self.parent.positions[self.rel_end - 1]
    }

    #[allow(dead_code)]
    pub fn positions(&self) -> &[PosType] {
        &self.parent.positions[self.rel_start..self.rel_end]
    }

    pub fn size(&self) -> usize {
        self.rel_end - self.rel_start
    }

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
        debug_assert!(other.rel_start.saturating_sub(self.rel_end) < 2);
        SegmentView::new(self.rel_start, other.rel_end, self.parent.clone())
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
            let (_, prob) = mann_whitney_u(self.group_a(), self.group_b());
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

#[allow(clippy::non_canonical_partial_ord_impl)]
impl PartialOrd for SegmentView<'_> {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {
        if self.parent == other.parent {
            self.rel_start.partial_cmp(&other.rel_start)
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
    group_a:   Vec<DensityType>,
    group_b:   Vec<DensityType>,
    mds_orig:  Vec<DensityType>,
    positions: Vec<PosType>,
}

impl TryFrom<(&BsxBatch, &BsxBatch)> for SegmentOwned {
    type Error = anyhow::Error;

    fn try_from(value: (&BsxBatch, &BsxBatch)) -> Result<Self, Self::Error> {
        let (left, right) = value;

        if right.len() != left.len() {
            anyhow::bail!("Batches have different lengths")
        }
        if right.as_contig().unwrap() != left.as_contig().unwrap() {
            anyhow::bail!("Batches cover different regions")
        }

        let positions = left
            .position()
            .into_no_null_iter()
            .map(|v| v as PosType)
            .collect_vec();

        let left_density = left
            .density()
            .into_iter()
            .map(|v| v.unwrap_or(DensityType::NAN))
            .collect_vec();
        let right_density = right
            .density()
            .into_iter()
            .map(|v| v.unwrap_or(DensityType::NAN))
            .collect_vec();

        Ok(Self::new(positions, left_density, right_density))
    }
}

impl SegmentOwned {
    pub fn new(
        positions: Vec<PosType>,
        group_a: Vec<DensityType>,
        group_b: Vec<DensityType>,
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
        max_dist: PosType,
    ) -> Vec<Self> {
        let mut split_idxs = segmentation::arg_split_segment(&self.positions, max_dist);
        split_idxs.reverse();
        if split_idxs.is_empty() {
            vec![self]
        }
        else {
            let mut res = split_idxs.into_iter().fold(Vec::new(), |mut acc, idx| {
                let positions = self.positions.drain(idx..).collect_vec();
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
        self.mds_orig.append(&mut other.mds_orig);
        self.positions.append(&mut other.positions);
        self
    }

    pub fn to_view(&self) -> SegmentView {
        SegmentView::new(0, self.mds_orig.len(), Arc::new(self))
    }

    pub fn size(&self) -> usize {
        self.positions.len()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DMRegion {
    pub chr:         String,
    pub start:       PosType,
    pub end:         PosType,
    #[serde(serialize_with = "serialize_scientific")]
    pub p_value:     f64,
    pub meth_left:   DensityType,
    pub meth_right:  DensityType,
    pub n_cytosines: usize,
    pub meth_diff:   DensityType,
    pub meth_mean:   DensityType,
}

impl DMRegion {
    pub(crate) fn from_segment_view(
        segment: SegmentView,
        chr: String,
    ) -> Self {
        let a_mean = segment.group_a().iter().sum::<f32>() / segment.size() as f32;
        let b_mean = segment.group_b().iter().sum::<f32>() / segment.size() as f32;

        DMRegion {
            chr,
            start: segment.start_pos(),
            end: segment.end_pos(),
            p_value: segment.get_pvalue(),
            meth_left: a_mean,
            meth_right: b_mean,
            n_cytosines: segment.size(),
            meth_diff: a_mean - b_mean,
            meth_mean: (a_mean + b_mean) / 2.0,
        }
    }

    #[allow(dead_code)]
    fn meth_mean(&self) -> DensityType {
        (self.meth_left + self.meth_right) / 2.0
    }

    pub fn as_contig(&self) -> Contig {
        Contig::new(self.chr.as_str().into(), self.start, self.end, Strand::None)
    }
}

fn serialize_scientific<S>(
    x: &f64,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer, {
    serializer.serialize_str(&format!("{:e}", x))
}
