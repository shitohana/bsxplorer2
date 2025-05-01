use std::{
    collections::{BTreeMap, BTreeSet},
    io::{Read, Seek},
};

use anyhow::{anyhow, bail};
use itertools::Itertools;
use polars::prelude::{search_sorted::binary_search_ca, SearchSortedSide};

use crate::data_structs::{
    batch::{BsxBatchBuilder, BsxBatchMethods, EncodedBsxBatch},
    coords::Contig,
};

use super::{read::BatchIndex, BsxFileReader};

pub struct RegionReader<R: Read + Seek> {
    cache: BTreeMap<usize, EncodedBsxBatch>,
    inner: BsxFileReader<R>,
}

impl<R: Read + Seek> RegionReader<R> {
    pub fn new(reader: BsxFileReader<R>) -> Self {
        Self { inner: reader, cache: BTreeMap::new() }
    }

    pub fn index(&mut self) -> anyhow::Result<&BatchIndex<String, u32>> {
        self.inner.index()
    }

    pub fn reset(&mut self) {
        self.cache.clear();
    }

    pub fn query<S: AsRef<str> + Clone + Into<String>>(
        &mut self,
        contig: Contig<S, u32>,
    ) -> anyhow::Result<Option<EncodedBsxBatch>> {
        const MAX_ITERATION_LIMIT: usize = 10;
        let mut depth = 0usize;

        if !self.inner.index()?.chr_order().contains(&contig.seqname().as_ref().to_string()) {
            bail!("Contig seqname not found in index")
        }

        loop {
            if depth > MAX_ITERATION_LIMIT {
                bail!("Maximum recursion depth exceeded");
            }
            depth += 1;

            let min_cached_pos = self
                .cache
                .first_key_value()
                .map(|(_k, v)| v.start_gpos())
                .transpose()?;
            let max_cached_pos = self
                .cache
                .last_key_value()
                .map(|(_k, v)| v.end_gpos())
                .transpose()?;

            let intersection_kind = if let Some((min_pos, max_pos)) =
                min_cached_pos.zip(max_cached_pos)
            {
                assert_eq!(min_pos.seqname(), max_pos.seqname());
                if min_pos.seqname() == contig.seqname().as_ref()
                    && min_pos.position() <= contig.start()
                    && max_pos.position() >= contig.end()
                {
                    IntersectionKind::Full
                } else if min_pos.seqname() == contig.seqname().as_ref()
                    && min_pos.position() <= contig.start()
                    && max_pos.position() < contig.end()
                {
                    IntersectionKind::PartialRight
                } else if min_pos.seqname() == contig.seqname().as_ref()
                    && min_pos.position() > contig.start()
                    && max_pos.position() >= contig.end()
                {
                    IntersectionKind::PartialLeft
                } else {
                    IntersectionKind::None
                }
            } else {
                IntersectionKind::None
            };

            match intersection_kind {
                IntersectionKind::PartialLeft => {
                    bail!("Required batch has already been processed. Make sure the regions are sorted with BatchIndex.sort")
                },
                IntersectionKind::Full => {
                    let batches = self
                        .inner
                        .index()?
                        .find(&contig.clone().cast(
                            |seqname| seqname.as_ref().to_string(),
                            |pos| pos,
                        ))
                        .unwrap()
                        .into_iter()
                        .sorted()
                        .map(|idx| self.cache.get(&idx))
                        .collect::<Option<Vec<&EncodedBsxBatch>>>()
                        .ok_or(anyhow!("Batch data missing"))?;

                    let batches_total = batches.len();
                    let mut read_batches = 0usize;

                    let mut res = vec![];
                    for batch in batches {
                        let is_first = read_batches == 0;
                        let is_last = read_batches == batches_total - 1;
                        let slice_start = if is_first {
                            binary_search_ca(
                                batch.position(),
                                [Some(contig.start())].into_iter(),
                                SearchSortedSide::Left,
                                false,
                            )[0]
                        } else {
                            0
                        };

                        let slice_end = if is_last {
                            binary_search_ca(
                                batch.position(),
                                [Some(contig.end())].into_iter(),
                                SearchSortedSide::Left,
                                false,
                            )[0]
                        } else {
                            batch.height() as u32
                        };

                        let slice = batch.slice(slice_start, slice_end - slice_start);
                        res.push(slice);
                        read_batches += 1;
                    }
                    res.sort_by_key(|b| b.start_pos());
                    if res.is_empty() {
                        return Ok(None);
                    } else if res.len() == 1 {
                        return Ok(Some(res.pop().unwrap()));
                    } else {
                        return Some(BsxBatchBuilder::concat(res)).transpose();
                    }
                },
                IntersectionKind::PartialRight | IntersectionKind::None => {
                    let required_batches =
                        self.inner.index()?.find(&contig.clone().cast(
                            |seqname| seqname.as_ref().to_string(),
                            |pos| pos,
                        ));
                    if required_batches.is_none()
                        || (required_batches.is_some()
                            && required_batches
                                .as_ref()
                                .unwrap()
                                .is_empty())
                    {
                        return Ok(None);
                    }
                    let required_batches =
                        BTreeSet::from_iter(required_batches.unwrap());
                    let existing_batches =
                        BTreeSet::from_iter(self.cache.keys().cloned());

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
                        self.cache.insert(
                            idx,
                            self.inner
                                .get_batch(idx)
                                .transpose()?
                                .ok_or_else(|| anyhow!("Batch index {} reported by find() but not found by get_batch()", idx))?,
                        );
                    }
                    continue;
                },
            }
        }
    }
}

enum IntersectionKind {
    Full,
    PartialRight,
    PartialLeft,
    None,
}
