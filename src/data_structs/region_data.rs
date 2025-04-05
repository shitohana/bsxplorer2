use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use crate::data_structs::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use crate::data_structs::region::RegionCoordinates;
use crate::utils::types::{Data, PosNum, RefId, Strand};

/// Represents data_structs associated with a genomic region.
#[derive(Clone, Eq, PartialEq)]
pub struct RegionData<R, N, T>
where
    R: RefId,
    N: PosNum,
    T: Data, {
    /// The chromosome ID.
    chr:        R,
    /// The start position of the region.
    start:      N,
    /// The end position of the region.
    end:        N,
    /// The strand of the region.
    strand:     Strand,
    /// The associated data_structs.
    data:       T,
    /// Optional attributes associated with the region.
    attributes: HashMap<String, String>,
}

impl<R: RefId, N: PosNum, T: Data> Hash for RegionData<R, N, T> {
    /// Hashes the region data_structs based on chromosome, start, end, and
    /// strand.
    fn hash<H: Hasher>(
        &self,
        state: &mut H,
    ) {
        (self.chr.clone(), self.start, self.end, self.strand).hash(state)
    }
}

impl<R, N> From<RegionCoordinates<N>> for RegionData<R, N, ()>
where
    R: RefId,
    N: PosNum,
{
    /// Creates a new `RegionData` from region coordinates with no data_structs.
    fn from(region: RegionCoordinates<N>) -> Self {
        RegionData {
            chr:        region.chr.into(),
            start:      region.start,
            end:        region.end,
            strand:     region.strand,
            data:       (),
            attributes: HashMap::default(),
        }
    }
}

impl<R, N, T> RegionData<R, N, T>
where
    R: RefId,
    N: PosNum,
    T: Data,
{
    /// Returns a reference to the chromosome ID.
    pub fn chr(&self) -> &R { &self.chr }

    /// Returns the start position of the region.
    pub fn start(&self) -> N { self.start }

    /// Returns the end position of the region.
    pub fn end(&self) -> N { self.end }

    /// Returns the strand of the region.
    pub fn strand(&self) -> Strand { self.strand }

    /// Returns a reference to the associated data_structs.
    pub fn data(&self) -> &T { &self.data }

    /// Sets the chromosome ID.
    pub fn set_chr(
        &mut self,
        chr: R,
    ) {
        self.chr = chr;
    }

    /// Sets the start position.
    pub fn set_start(
        &mut self,
        start: N,
    ) {
        self.start = start;
    }

    /// Sets the end position.
    pub fn set_end(
        &mut self,
        end: N,
    ) {
        self.end = end;
    }

    /// Sets the strand.
    pub fn set_strand(
        &mut self,
        strand: Strand,
    ) {
        self.strand = strand;
    }

    /// Sets the associated data_structs.
    pub fn set_data(
        &mut self,
        data: T,
    ) {
        self.data = data;
    }

    /// Returns the length of the region.
    pub fn length(&self) -> N { self.end.clone() - self.start.clone() }

    /// Creates a new `RegionData` with the same region information but
    /// different data_structs.
    pub fn with_data<T2: Data>(
        self,
        data: T2,
    ) -> RegionData<R, N, T2> {
        RegionData {
            chr: self.chr,
            start: self.start,
            end: self.end,
            strand: self.strand,
            data,
            attributes: self.attributes,
        }
    }

    /// Creates a new `RegionData` instance.
    pub fn new(
        chr: R,
        start: N,
        end: N,
        strand: Strand,
        data: T,
        attributes: HashMap<String, String>,
    ) -> Self {
        RegionData {
            chr,
            start,
            end,
            strand,
            data,
            attributes,
        }
    }

    /// Creates a new `RegionData` from region coordinates and a value.
    pub fn from_region_and_value(
        region: RegionCoordinates<N>,
        value: T,
        attributes: HashMap<String, String>,
    ) -> RegionData<String, N, T> {
        RegionData {
            chr: region.chr,
            start: region.start,
            end: region.end,
            strand: region.strand,
            data: value,
            attributes,
        }
    }

    /// Tries to create a `RegionData` instance from an `EncodedBsxBatch`.
    pub fn try_from_batch(
        batch: &EncodedBsxBatch
    ) -> anyhow::Result<RegionData<String, u64, &EncodedBsxBatch>> {
        let first_position = batch.first_position()?;

        Ok(RegionData {
            chr:        first_position.chr().to_string(),
            start:      first_position.position(),
            end:        batch.last_position()?.position(),
            strand:     Strand::None,
            data:       batch,
            attributes: HashMap::default(),
        })
    }

    /// Consumes the `RegionData` and returns the associated data_structs.
    pub fn into_data(self) -> T { self.data }

    /// Separates the region data_structs from the associated data_structs,
    /// returning both.
    pub fn drain_data(self) -> (RegionData<R, N, ()>, T) {
        let data = self.data.clone();
        (self.with_data(()), data)
    }
}

impl<R, N, T> RegionData<R, N, T>
where
    R: RefId,
    N: PosNum,
    T: Data,
    String: From<R>,
{
    /// Converts the `RegionData` into a `RegionCoordinates` instance.
    pub fn into_region(self) -> RegionCoordinates<N> {
        RegionCoordinates {
            chr:    self.chr.into(),
            start:  self.start,
            end:    self.end,
            strand: self.strand,
        }
    }

    /// Returns the region as a `RegionCoordinates` instance.
    pub fn as_region(&self) -> RegionCoordinates<N> {
        RegionCoordinates::new(
            self.chr.clone().into(),
            self.start,
            self.end,
            self.strand,
        )
    }
}

// /// ## Notes
// /// Region strand (i.e. reversion) handling is the responsibility of the caller
// impl<R, N> RegionData<R, N, EncodedBsxBatch>
// where
//     R: RefId,
//     N: PosNum,
// {
//     // /// Returns the methylation density values for the region.
//     // TODO move to ChunkedArray
//     // pub fn meth_density(&self) -> PolarsResult<Vec<f32>> {
//     //     self.data.get_density_vals()
//     // }
//
//     /// Returns the positions within the region.
//     pub fn get_positions(&self) -> PolarsResult<Vec<u32>> {
//         self.data.get_position_vals()
//     }
//
//     /// Returns the number of data_structs points in the associated
//     /// data_structs.
//     pub fn size(&self) -> usize { self.data.height() }
//
//     /// Discretizes the region into a specified number of bins.
//     pub fn discretize(
//         &self,
//         resolution: usize,
//     ) -> PolarsResult<Vec<usize>> {
//         let length = self.length();
//         let positions_norm = self
//             .get_positions()?
//             .into_par_iter()
//             .map(|x| (x as f64) / length.to_f64().unwrap() * resolution as f64)
//             .collect::<Vec<f64>>();
//         let mut resolution_iter = (0..resolution).into_iter().peekable();
//         let sort_indices = positions_norm
//             .into_iter()
//             .enumerate()
//             .filter_map(|(idx, x)| {
//                 match resolution_iter.peek() {
//                     Some(&y) if x > y as f64 => {
//                         resolution_iter.next().unwrap();
//                         Some(idx)
//                     },
//                     _ => None,
//                 }
//             })
//             .collect_vec();
//         Ok(sort_indices)
//     }
//
//     // /// Calculates the discrete methylation density for a specified resolution.
//     // TODO move to ChunkedArray
//     // pub fn discrete_density(
//     //     &self,
//     //     resolution: usize,
//     // ) -> PolarsResult<Vec<f32>> {
//     //     let sort_indices = self.discretize(resolution)?;
//     //     let cum_meth_density = self
//     //         .meth_density()?
//     //         .into_iter()
//     //         .enumerate()
//     //         .fold(vec![0.0; self.size()], |mut acc, (idx, new)| {
//     //             acc[idx] = new;
//     //             if idx > 0 {
//     //                 acc[idx] += acc[idx - 1]
//     //             };
//     //             acc
//     //         });
//     //     let mut density = vec![0.0; resolution];
//     //     for (idx, bound_idx) in sort_indices.into_iter().enumerate() {
//     //         density[idx] = cum_meth_density[bound_idx];
//     //     }
//     //     Ok(density)
//     // }
//
//     /// Calculates the discrete methylation counts for a specified resolution.
//     pub fn discrete_counts(
//         &self,
//         resolution: usize,
//     ) -> PolarsResult<Vec<(i16, i16)>> {
//         let sort_indices = self.discretize(resolution)?;
//         let cum_counts_meth = self
//             .counts_meth()?
//             .into_iter()
//             .enumerate()
//             .fold(vec![0; self.size()], |mut acc, (idx, new)| {
//                 acc[idx] = new;
//                 if idx > 0 {
//                     acc[idx] += acc[idx - 1]
//                 };
//                 acc
//             });
//         let cum_counts_total = self
//             .counts_total()?
//             .into_iter()
//             .enumerate()
//             .fold(vec![0; self.size()], |mut acc, (idx, new)| {
//                 acc[idx] = new;
//                 if idx > 0 {
//                     acc[idx] += acc[idx - 1]
//                 };
//                 acc
//             });
//         let mut counts = vec![(0, 0); resolution];
//         for (idx, (meth, total)) in sort_indices.into_iter().zip(
//             cum_counts_meth
//                 .into_iter()
//                 .zip(cum_counts_total.into_iter()),
//         ) {
//             counts[idx] = (meth, total);
//         }
//         Ok(counts)
//     }
//
//     /// Returns the methylated counts for each position in the region.
//     pub fn counts_meth(&self) -> PolarsResult<Vec<i16>> {
//         self.data.get_counts_m()
//     }
//
//     /// Returns the total counts for each position in the region.
//     pub fn counts_total(&self) -> PolarsResult<Vec<i16>> {
//         self.data.get_counts_total()
//     }
//
//     /// Filters the data_structs based on the provided context.
//     pub fn filter(
//         mut self,
//         context: Context,
//     ) -> Self {
//         self.data = self.data.filter(Some(context), None);
//         self
//     }
//
//     /// Separates the data_structs by context (CG, CHG, CHH).
//     pub fn strip_contexts(self) -> anyhow::Result<(Self, Self, Self)> {
//         let (drained, data) = self.drain_data();
//         let (cg, chg, chh) = data.strip_contexts()?;
//         Ok((
//             drained.clone().with_data(cg),
//             drained.clone().with_data(chg),
//             drained.clone().with_data(chh),
//         ))
//     }
// }
