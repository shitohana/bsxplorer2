use std::collections::{BTreeMap, BTreeSet};
use std::io::{Read, Seek};

use anyhow::{anyhow, bail};
use itertools::Itertools;
use polars::prelude::search_sorted::binary_search_ca;
use polars::prelude::SearchSortedSide;

use super::BsxFileReader;
use crate::data_structs::batch::{BsxBatchBuilder,
                                 BsxBatchMethods,
                                 EncodedBsxBatch};
use crate::data_structs::coords::Contig;
use crate::data_structs::typedef::{SeqNameStr, SeqPosNum};
use crate::io::bsx::BatchIndex;

/// RegionReader is a reader for BSX files that operates on a specific region of
/// the genome.
pub struct RegionReader<R, S, P>
where
    R: Read + Seek,
    S: SeqNameStr,
    P: SeqPosNum, {
    /// Cache of encoded BSX batches.
    cache:         BTreeMap<usize, EncodedBsxBatch>,
    /// Inner reader for the BSX file.
    inner:         BsxFileReader<R>,
    /// Index of the BSX file.
    index:         BatchIndex<S, P>,
    /// Preprocessing function to be applied to each batch before it is cached.
    preprocess_fn:
        Option<Box<dyn Fn(EncodedBsxBatch) -> anyhow::Result<EncodedBsxBatch>>>,
}

// PRIVATE METHODS
impl<R, S, P> RegionReader<R, S, P>
where
    R: Read + Seek,
    S: SeqNameStr,
    P: SeqPosNum,
{
    /// Finds the batches that overlap the given contig.
    fn find(
        &self,
        contig: &Contig<S, P>,
    ) -> Option<Vec<usize>> {
        let required_batches = self.index().find(&contig.clone());

        let batches_found = required_batches
            .as_ref()
            .map(|v| v.len() > 0)
            .unwrap_or(false);

        if batches_found {
            required_batches
        }
        else {
            None
        }
    }

    /// Updates the cache with the given batches.
    fn update_cache(
        &mut self,
        batches: &[usize],
    ) -> anyhow::Result<()> {
        let required_batches = BTreeSet::from_iter(batches.iter().cloned());
        let existing_batches = BTreeSet::from_iter(self.cache.keys().cloned());

        let min_required = required_batches.first().unwrap();

        let to_discard = existing_batches
            .difference(&required_batches)
            .cloned()
            .filter(|idx| idx < min_required)
            .collect_vec();
        let to_read = required_batches
            .difference(&existing_batches)
            .cloned()
            .collect_vec();

        for idx in to_discard {
            self.cache.remove(&idx);
        }

        for idx in to_read {
            let mut data = self
                .inner
                .get_batch(idx)
                .transpose()?
                .ok_or_else(|| {
                    anyhow!(
                        "Batch index {} reported by find() but not found by \
                         get_batch()",
                        idx
                    )
                })?;

            if let Some(postprocess_fn) = self.preprocess_fn.as_ref() {
                data = postprocess_fn(data)?;
            }

            self.cache.insert(idx, data);
        }
        Ok(())
    }

    /// Gets the batches for the given contig from the cache.
    fn get_batches_for_contig(
        &self,
        contig: &Contig<S, P>,
    ) -> anyhow::Result<Vec<&EncodedBsxBatch>> {
        self.index()
            .find(&contig.clone())
            .unwrap()
            .into_iter()
            .sorted()
            .map(|idx| self.cache.get(&idx))
            .collect::<Option<Vec<&EncodedBsxBatch>>>()
            .ok_or(anyhow!("Batch data missing"))
    }

    /// Assembles the region from the given batches.
    ///
    /// This method assumes that the cache is up-to-date and that the required
    /// batches are present.
    fn assemble_region(
        &self,
        contig: &Contig<S, P>,
        batches: Vec<&EncodedBsxBatch>,
    ) -> anyhow::Result<Option<EncodedBsxBatch>> {
        let batches_total = batches.len();

        let mut res = vec![];
        for (read_batches, batch) in batches.into_iter().enumerate() {
            let is_first = read_batches == 0;
            let is_last = read_batches == batches_total - 1;
            let slice_start = if is_first {
                binary_search_ca(
                    batch.position(),
                    [Some(contig.start().to_u32().unwrap())].into_iter(),
                    SearchSortedSide::Left,
                    false,
                )[0]
            }
            else {
                0
            };

            let slice_end = if is_last {
                binary_search_ca(
                    batch.position(),
                    [Some(contig.end().to_u32().unwrap())].into_iter(),
                    SearchSortedSide::Left,
                    false,
                )[0]
            }
            else {
                batch.height() as u32
            };

            let slice = batch.slice(slice_start, slice_end - slice_start);
            res.push(slice);
        }
        res.sort_by_key(|b| b.start_pos());
        if res.is_empty() {
            Ok(None)
        }
        else if res.len() == 1 {
            Ok(Some(res.pop().unwrap()))
        }
        else {
            Some(BsxBatchBuilder::concat(res)).transpose()
        }
    }

    /// Determines the intersection kind between the cached region and the query
    /// contig.
    fn determine_intersection(
        &mut self,
        contig: &Contig<S, P>,
    ) -> anyhow::Result<IntersectionKind> {
        let min_cached_pos = self
            .cache
            .first_key_value()
            .map(|(_k, v)| v.start_gpos::<String, u32>());
        let max_cached_pos = self
            .cache
            .last_key_value()
            .map(|(_k, v)| v.end_gpos::<String, u32>());

        let intersection_kind = if let Some((min_pos, max_pos)) =
            min_cached_pos.zip(max_cached_pos)
        {
            assert_eq!(min_pos.seqname(), max_pos.seqname());
            let seqname = min_pos.seqname();
            let min_pos_val = P::from(min_pos.position()).unwrap();
            let max_pos_val = P::from(max_pos.position()).unwrap();

            // |-------------<cache>------------|
            //        |------<contig>----|
            if seqname == contig.seqname().as_ref()
                && contig.start() >= min_pos_val
                && contig.end() <= max_pos_val
            {
                IntersectionKind::Full
            }
            //        |-------------<cache>------------|
            // |--------<contig>----|
            // But this will qualify as PartialRight too, even though it
            // is not a partial intersection
            //        |-------------<cache>------------|
            // tig>-|
            else if min_pos.seqname() == contig.seqname().as_ref()
                && contig.start() <= min_pos_val
                && contig.end() <= max_pos_val
            {
                IntersectionKind::PartialRight
            }
            // |-------------<cache>------------|
            //                  |----<region>------|
            else if min_pos.seqname() == contig.seqname().as_ref()
                && contig.start() > min_pos_val
                && contig.start() < max_pos_val
                && contig.end() >= max_pos_val
            // |----<cache>----|
            //                   |----<contig>-----|
            {
                IntersectionKind::PartialLeft
            }
            else {
                IntersectionKind::None
            }
        }
        else {
            IntersectionKind::None
        };

        Ok(intersection_kind)
    }
}

// PUBLIC METHODS
impl<R> RegionReader<R, String, u32>
where
    R: Read + Seek,
{
    pub fn from_reader(reader: BsxFileReader<R>) -> anyhow::Result<Self> {
        let mut reader = reader;
        let index = reader.index()?.clone();
        Ok(Self::new(reader, index, None))
    }
}

impl<R, S, P> RegionReader<R, S, P>
where
    R: Read + Seek,
    S: SeqNameStr,
    P: SeqPosNum,
{
    /// Sets the preprocessing function.
    pub fn set_preprocess_fn(
        &mut self,
        preprocess_fn: Option<
            Box<dyn Fn(EncodedBsxBatch) -> anyhow::Result<EncodedBsxBatch>>,
        >,
    ) {
        self.preprocess_fn = preprocess_fn;
    }

    /// Creates a new RegionReader.
    pub fn new(
        inner: BsxFileReader<R>,
        index: BatchIndex<S, P>,
        preprocess_fn: Option<
            fn(EncodedBsxBatch) -> anyhow::Result<EncodedBsxBatch>,
        >,
    ) -> Self {
        Self {
            cache: BTreeMap::new(),
            inner,
            index,
            preprocess_fn: preprocess_fn.map(|c| {
                Box::new(c)
                    as Box<
                        dyn Fn(
                            EncodedBsxBatch,
                        )
                            -> anyhow::Result<EncodedBsxBatch>,
                    >
            }),
        }
    }

    /// Returns the index of the BSX file.
    pub fn index(&self) -> &BatchIndex<S, P> {
        &self.index
    }

    /// Resets the cache.
    pub fn reset(&mut self) {
        self.cache.clear();
    }

    /// Queries the BSX file for the given contig.
    pub fn query(
        &mut self,
        contig: Contig<S, P>,
        postprocess_fn: Option<
            Box<fn(EncodedBsxBatch) -> anyhow::Result<EncodedBsxBatch>>,
        >,
    ) -> anyhow::Result<Option<EncodedBsxBatch>> {
        const MAX_ITERATION_LIMIT: usize = 10;
        let mut depth = 0usize;

        if !self
            .index()
            .get_chr_order()
            .contains(contig.seqname())
        {
            bail!("Contig seqname not found in index")
        }

        loop {
            if depth > MAX_ITERATION_LIMIT {
                bail!("Maximum recursion depth exceeded");
            }
            depth += 1;

            match self.determine_intersection(&contig)? {
                IntersectionKind::PartialLeft => {
                    bail!(
                        "Required batch has already been processed. Make sure \
                         the regions are sorted with BatchIndex.sort"
                    )
                },
                IntersectionKind::Full => {
                    let batches = self.get_batches_for_contig(&contig)?;
                    let res = self.assemble_region(&contig, batches)?;
                    if let Some(data) = res {
                        if let Some(postprocess_fn) = postprocess_fn {
                            return Some(postprocess_fn(data)).transpose();
                        }
                        else {
                            return Some(Ok(data)).transpose();
                        }
                    }
                    else {
                        return Ok(None);
                    }
                },
                IntersectionKind::PartialRight | IntersectionKind::None => {
                    if let Some(required_batches) = self.find(&contig) {
                        self.update_cache(&required_batches)?;
                        continue;
                    }
                    else {
                        return Ok(None);
                    }
                },
            }
        }
    }
}

/// Enum representing the intersection kind between the cached region and the
/// query contig.
enum IntersectionKind {
    Full,
    PartialRight,
    PartialLeft,
    None,
}
