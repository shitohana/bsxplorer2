use std::collections::{
    BTreeMap,
    BTreeSet,
    VecDeque,
};

use anyhow::{
    anyhow,
    bail,
};
use itertools::Itertools;
use polars::error::PolarsResult;
use polars::prelude::search_sorted::binary_search_ca;
use polars::prelude::SearchSortedSide;

use super::BsxFileReader;
use crate::data_structs::batch::{
    BsxBatch,
    BsxBatchBuilder,
};
use crate::data_structs::coords::Contig;
use crate::io::bsx::BatchIndex;

/// RegionReader is a reader for BSX files that operates on a specific region of
/// the genome.
type PreprocessFn =
    Box<dyn Fn(BsxBatch) -> anyhow::Result<BsxBatch> + Send + Sync + 'static>;

pub struct RegionReader {
    /// Cache of encoded BSX batches.
    cache:         BTreeMap<usize, BsxBatch>,
    /// Inner reader for the BSX file.
    inner:         BsxFileReader,
    /// Index of the BSX file.
    index:         BatchIndex,
    /// Preprocessing function to be applied to each batch before it is cached.
    preprocess_fn: Option<PreprocessFn>,
}

// PRIVATE METHODS
impl RegionReader {
    /// Finds the batches that overlap the given contig.
    fn find(
        &self,
        contig: &Contig,
    ) -> Option<Vec<usize>> {
        let required_batches = self.index().find(&contig.clone());

        let batches_found = required_batches
            .as_ref()
            .map(|v| !v.is_empty())
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

        let batches = self
            .inner
            .get_batches(&to_read)
            .into_iter()
            .collect::<Option<PolarsResult<Vec<_>>>>();
        if batches.is_none() {
            return Err(anyhow!("Some batches not found in file"));
        }
        let batches = batches.unwrap()?;

        for (idx, batch) in to_read.into_iter().zip(batches.into_iter()) {
            self.cache.insert(idx, batch);
        }

        Ok(())
    }

    /// Gets the batches for the given contig from the cache.
    fn get_batches_for_contig(
        &self,
        contig: &Contig,
    ) -> anyhow::Result<Vec<&BsxBatch>> {
        self.index()
            .find(&contig.clone())
            .unwrap()
            .into_iter()
            .sorted()
            .map(|idx| self.cache.get(&idx))
            .collect::<Option<Vec<&BsxBatch>>>()
            .ok_or(anyhow!("Batch data missing"))
    }

    /// Assembles the region from the given batches.
    ///
    /// This method assumes that the cache is up-to-date and that the required
    /// batches are present.
    fn assemble_region(
        &self,
        contig: &Contig,
        batches: Vec<&BsxBatch>,
    ) -> anyhow::Result<Option<BsxBatch>> {
        let batches_total = batches.len();

        let mut res = vec![];
        for (read_batches, batch) in batches.into_iter().enumerate() {
            let is_first = read_batches == 0;
            let is_last = read_batches == batches_total - 1;
            let slice_start = if is_first {
                binary_search_ca(
                    batch.position(),
                    [Some(contig.start())].into_iter(),
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
                    [Some(contig.end())].into_iter(),
                    SearchSortedSide::Left,
                    false,
                )[0]
            }
            else {
                batch.len() as u32
            };

            let slice =
                batch.slice(slice_start as i64, (slice_end - slice_start) as usize);
            res.push(slice);
        }
        res.sort_by_key(|b| b.first_pos());
        if res.is_empty() {
            Ok(None)
        }
        else if res.len() == 1 {
            Ok(Some(res.pop().unwrap()))
        }
        else {
            Some(BsxBatchBuilder::concat(res))
                .transpose()
                .map_err(|e| anyhow::anyhow!(e))
        }
    }

    /// Determines the intersection kind between the cached region and the query
    /// contig.
    fn determine_intersection(
        &mut self,
        contig: &Contig,
    ) -> anyhow::Result<IntersectionKind> {
        let min_cached_pos = self
            .cache
            .first_key_value()
            .map(|(_k, v)| v.first_genomic_pos());
        let max_cached_pos = self
            .cache
            .last_key_value()
            .map(|(_k, v)| v.last_genomic_pos());

        if min_cached_pos.is_none() || max_cached_pos.is_none() {
            return Ok(IntersectionKind::None);
        }
        let min_cached_pos = min_cached_pos.unwrap();
        let max_cached_pos = max_cached_pos.unwrap();

        let intersection_kind =
            if let Some((min_pos, max_pos)) = min_cached_pos.zip(max_cached_pos) {
                assert_eq!(min_pos.seqname(), max_pos.seqname());
                let seqname = min_pos.seqname();
                let min_pos_val = min_pos.position();
                let max_pos_val = max_pos.position();

                // |-------------<cache>------------|
                //        |------<contig>----|
                if seqname == contig.seqname()
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
                else if min_pos.seqname() == contig.seqname()
                    && contig.start() <= min_pos_val
                    && contig.end() <= max_pos_val
                {
                    IntersectionKind::PartialLeft
                }
                // |-------------<cache>------------|
                //                  |----<region>------|
                else if min_pos.seqname() == contig.seqname()
                    && contig.start() > min_pos_val
                    && contig.start() < max_pos_val
                    && contig.end() >= max_pos_val
                // |----<cache>----|
                //                   |----<contig>-----|
                {
                    IntersectionKind::PartialRight
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
impl RegionReader {
    /// Creates a new RegionReader from an existing BsxFileReader.
    ///
    /// Reads the index from the reader during initialization.
    pub fn from_reader(reader: BsxFileReader) -> anyhow::Result<Self> {
        let mut reader = reader;
        let index = BatchIndex::from_reader(&mut reader)?;
        Ok(Self::new(reader, index, None))
    }

    /// Sets the preprocessing function.
    pub fn set_preprocess_fn(
        &mut self,
        preprocess_fn: Option<PreprocessFn>,
    ) {
        self.preprocess_fn = preprocess_fn;
    }

    /// Creates a new RegionReader instance.
    ///
    /// # Arguments
    ///
    /// * `inner`: The inner BsxFileReader.
    /// * `index`: The BatchIndex for the BSX file.
    /// * `preprocess_fn`: An optional function to apply to each batch after
    ///   reading but before caching.
    pub fn new(
        inner: BsxFileReader,
        index: BatchIndex,
        preprocess_fn: Option<fn(BsxBatch) -> anyhow::Result<BsxBatch>>,
    ) -> Self {
        Self {
            cache: BTreeMap::new(),
            inner,
            index,
            preprocess_fn: preprocess_fn.map(|c| Box::new(c) as PreprocessFn),
        }
    }

    /// Returns the index of the BSX file.
    pub fn index(&self) -> &BatchIndex {
        &self.index
    }

    /// Resets the cache.
    pub fn reset(&mut self) {
        self.cache.clear();
    }

    /// Creates an iterator over a list of contigs.
    ///
    /// # Arguments
    ///
    /// * `contigs`: A slice of Contig regions to iterate over.
    pub fn iter_contigs(
        &mut self,
        contigs: &[Contig],
    ) -> RegionReaderIterator {
        RegionReaderIterator {
            reader:          self,
            pending_contigs: VecDeque::from(contigs.to_vec()),
            cached_contigs:  VecDeque::new(),
        }
    }

    /// Queries the BSX file for the given contig region.
    ///
    /// Attempts to retrieve the data for the specified contig, potentially
    /// reading batches from the file and caching them.
    ///
    /// # Arguments
    ///
    /// * `contig`: The Contig region to query.
    /// * `postprocess_fn`: An optional function to apply to the assembled batch
    ///   after querying.
    pub fn query(
        &mut self,
        contig: Contig,
        postprocess_fn: Option<PreprocessFn>,
    ) -> anyhow::Result<Option<BsxBatch>> {
        const MAX_ITERATION_LIMIT: usize = 10;
        let mut depth = 0usize;

        if !self.index().get_chr_order().contains(contig.seqname()) {
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
                        "Required batch has already been processed. Make sure the \
                         regions are sorted with BatchIndex.sort"
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

/// Iterator for reading BSX data for multiple contigs using a RegionReader.
pub struct RegionReaderIterator<'a> {
    reader:          &'a mut RegionReader,
    pending_contigs: VecDeque<Contig>,
    cached_contigs:  VecDeque<Contig>,
}

impl Iterator for RegionReaderIterator<'_> {
    type Item = anyhow::Result<BsxBatch>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(contig) = self.cached_contigs.pop_front() {
            self.reader.query(contig, None::<PreprocessFn>).transpose()
        }
        else {
            let fill_cache_res = self.fill_cache();
            if let Some(Ok(_)) = fill_cache_res {
                self.next()
            }
            else if let Some(Err(e)) = fill_cache_res {
                Some(Err(e))
            }
            else {
                None
            }
        }
    }
}

impl RegionReaderIterator<'_> {
    fn fill_cache(&mut self) -> Option<anyhow::Result<()>> {
        let mut required_batches = BTreeSet::new();
        let n_threads = self.reader.inner.n_threads();
        let mut cur_chr = None;

        while required_batches.len() < n_threads && !self.pending_contigs.is_empty() {
            let contig = self.pending_contigs.pop_front().unwrap();
            if cur_chr.is_some() && cur_chr.as_ref() != Some(contig.seqname()) {
                self.pending_contigs.push_front(contig);
                break;
            }
            else {
                cur_chr = Some(contig.seqname().clone());
            }

            if let Some(batches) = self.reader.find(&contig) {
                required_batches.append(&mut BTreeSet::from_iter(batches));
                self.cached_contigs.push_back(contig);
            }
            else {
                continue;
            };
        }

        if required_batches.is_empty() {
            None
        }
        else {
            let res = self
                .reader
                .update_cache(&required_batches.into_iter().collect_vec());
            Some(res)
        }
    }
}

/// Enum representing the intersection kind between the cached region and the
/// query contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum IntersectionKind {
    Full,
    PartialRight,
    PartialLeft,
    None,
}
