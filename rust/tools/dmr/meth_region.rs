use crate::tools::dmr::beta_binom::BetaBinomObservation;
use crate::tools::dmr::penalty_segment::PenaltySegmentModel;
use crate::tools::dmr::sure_segment;
use crate::tools::dmr::sure_segment::SureSegmentModelConfig;
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// An owned methylation region. This struct owns its data.
#[derive(Clone, Debug)]
pub struct MethylatedRegionOwned {
    pub positions: Vec<u32>,
    pub density: Vec<f64>,
    pub count_m: Vec<u32>,
    pub count_total: Vec<u32>,
}

impl MethylatedRegionOwned {
    /// Creates a new region, validating that all vectors have the same nonzero length.
    pub fn new(
        positions: Vec<u32>,
        density: Vec<f64>,
        count_m: Vec<u32>,
        count_total: Vec<u32>,
    ) -> Self {
        assert!(!positions.is_empty(), "Region cannot be empty");
        assert_eq!(positions.len(), density.len());
        assert_eq!(positions.len(), count_m.len());
        assert_eq!(positions.len(), count_total.len());
        Self {
            positions,
            density,
            count_m,
            count_total,
        }
    }

    /// Returns a read-only view of the entire region.
    pub fn as_view(&self) -> MethylatedRegionView {
        MethylatedRegionView {
            positions: &self.positions,
            density: &self.density,
            count_m: &self.count_m,
            count_total: &self.count_total,
        }
    }

    /// Returns a read-only view for the slice [start, end) without copying.
    pub fn slice(&self, start: usize, end: usize) -> MethylatedRegionView {
        assert!(start <= end, "Start must be <= end");
        assert!(
            end <= self.positions.len(),
            "Slice bounds exceed region size"
        );
        MethylatedRegionView {
            positions: &self.positions[start..end],
            density: &self.density[start..end],
            count_m: &self.count_m[start..end],
            count_total: &self.count_total[start..end],
        }
    }

    pub fn slice_multiple(&self, indices: &[(usize, usize)]) -> Vec<MethylatedRegionView> {
        indices
            .into_par_iter()
            .map(|(start, end)| self.slice(*start, end + 1))
            .collect()
    }
}

/// A lightweight view into a methylation region.
/// This struct borrows from an underlying owned region without extra allocations.
#[derive(Clone, Debug)]
pub struct MethylatedRegionView<'a> {
    pub positions: &'a [u32],
    pub density: &'a [f64],
    pub count_m: &'a [u32],
    pub count_total: &'a [u32],
}

impl<'a> MethylatedRegionView<'a> {
    pub fn to_owned(self) -> MethylatedRegionOwned {
        MethylatedRegionOwned {
            positions: self.positions.to_vec(),
            density: self.density.to_vec(),
            count_m: self.count_m.to_vec(),
            count_total: self.count_total.to_vec(),
        }
    }
    pub fn size(&self) -> usize {
        self.positions.len()
    }

    /// Returns the number of sites in this view.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn segment_indices_sure(&self, config: &SureSegmentModelConfig) -> Vec<(usize, usize)> {
        let segment_indices = sure_segment::segment_signal_sure(&self.density, config);
        segment_indices
            .into_par_iter()
            .map(|(start, end)| {
                let mut result = Vec::new();
                let mut current = start;
                for i in start..end {
                    if self.positions[i + 1] - self.positions[i] > config.max_dist {
                        result.push((current, i));
                        current = i + 1
                    }
                }

                result.push((current, end));
                result
            })
            .flatten()
            .collect()
    }

    pub fn segment_indices_penalty(&self, config: &PenaltySegmentModel) -> Vec<(usize, usize)> {
        let (_, denoised) = config.denoise(self.density);
        let segment_indices = config.get_segments(&denoised);
        let segments_merged = config
            .merge_segments(self.density, &segment_indices)
            .into_iter()
            .map(|(a, b, _)| (a, b))
            .collect_vec();

        segments_merged
            .into_par_iter()
            .map(|(start, end)| {
                let mut result = Vec::new();
                let mut current = start;
                for i in start..end {
                    if self.positions[i + 1] - self.positions[i] > config.max_dist {
                        result.push((current, i));
                        current = i + 1
                    }
                }

                result.push((current, end));
                result
            })
            .flatten()
            .collect()
    }

    /// Returns a read-only view for the slice [start, end) without copying.
    pub fn slice(&self, start: usize, end: usize) -> MethylatedRegionView {
        assert!(start <= end, "Start must be <= end");
        assert!(
            end <= self.positions.len(),
            "Slice bounds exceed region size"
        );
        MethylatedRegionView {
            positions: &self.positions[start..end],
            density: &self.density[start..end],
            count_m: &self.count_m[start..end],
            count_total: &self.count_total[start..end],
        }
    }

    pub fn slice_multiple(&self, indices: &[(usize, usize)]) -> Vec<MethylatedRegionView> {
        indices
            .into_par_iter()
            .map(|(start, end)| self.slice(*start, end + 1))
            .collect()
    }

    /// Splits the view at the provided index into two sub-views.
    pub fn split_at(&self, index: usize) -> (MethylatedRegionView<'a>, MethylatedRegionView<'a>) {
        let (pos_left, pos_right) = self.positions.split_at(index);
        let (dens_left, dens_right) = self.density.split_at(index);
        let (cm_left, cm_right) = self.count_m.split_at(index);
        let (ct_left, ct_right) = self.count_total.split_at(index);
        (
            MethylatedRegionView {
                positions: pos_left,
                density: dens_left,
                count_m: cm_left,
                count_total: ct_left,
            },
            MethylatedRegionView {
                positions: pos_right,
                density: dens_right,
                count_m: cm_right,
                count_total: ct_right,
            },
        )
    }

    /// Concatenates two views into a new owned region.
    /// Note: This operation does incur an allocation.
    pub fn concat(
        left: MethylatedRegionView<'a>,
        right: MethylatedRegionView<'a>,
    ) -> MethylatedRegionOwned {
        let mut positions = Vec::with_capacity(left.len() + right.len());
        positions.extend_from_slice(left.positions);
        positions.extend_from_slice(right.positions);

        let mut density = Vec::with_capacity(left.len() + right.len());
        density.extend_from_slice(left.density);
        density.extend_from_slice(right.density);

        let mut count_m = Vec::with_capacity(left.len() + right.len());
        count_m.extend_from_slice(left.count_m);
        count_m.extend_from_slice(right.count_m);

        let mut count_total = Vec::with_capacity(left.len() + right.len());
        count_total.extend_from_slice(left.count_total);
        count_total.extend_from_slice(right.count_total);

        MethylatedRegionOwned {
            positions,
            density,
            count_m,
            count_total,
        }
    }

    pub fn to_beta_binom_observations(&self) -> Vec<BetaBinomObservation> {
        self.count_m
            .iter()
            .zip(self.count_total.iter())
            .map(|(&count_m, &count_total)| {
                BetaBinomObservation::new(count_m as f64, count_total as f64).unwrap()
            })
            .collect()
    }
}
