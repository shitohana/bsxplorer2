use crate::data_structs::batch::EncodedBsxBatch;
use crate::data_structs::enums::Context;
use crate::io::bsx::BsxFileReader;
use crate::tools::dmr::data_structs::ReaderMetadata;
use crate::tools::dmr::segmentation::FilterConfig;
use crate::tools::dmr::{segment_reading, DmrIterator};

use bio_types::annot::refids::RefIDSet;
use itertools::Itertools;
use polars::io::mmap::MmapBytesReader;
use std::collections::BTreeMap;
use std::fmt::Display;
use std::hash::Hash;
use std::io::{Read, Seek};
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct DmrConfig {
    pub context: Context,
    pub n_missing: usize,
    pub min_coverage: i16,
    pub diff_threshold: f32,
    pub min_cpgs: usize,
    pub max_dist: u64,
    pub initial_l: f64,
    pub l_min: f64,
    pub l_coef: f64,
    pub seg_tolerance: f32,
    pub merge_pvalue: f64,
    pub seg_pvalue: f64,
}

impl DmrConfig {
    pub fn filter_config(&self) -> FilterConfig {
        FilterConfig {
            context: self.context,
            n_missing: self.n_missing,
            min_coverage: self.min_coverage,
        }
    }

    pub fn try_finish<F, R>(
        &self,
        readers: Vec<(R, F)>,
    ) -> anyhow::Result<DmrIterator>
    where
        F: Read + Seek + Send + Sync + MmapBytesReader + 'static,
        R: Display
            + Eq
            + Hash
            + Clone
            + Default
            + std::fmt::Debug
            + Send
            + Ord
            + 'static,
    {
        // Use BTreeMap to sort labels and get consistent left/right order
        let mut readers_pair = BTreeMap::from_iter(
            readers
                .into_iter()
                .into_group_map_by(|(label, reader)| label.clone()),
        )
        .into_values()
        .collect_vec();

        match readers_pair.len() {
            0 => anyhow::bail!("No readers supplied"),
            2 => {},
            n => anyhow::bail!("Too many readers groups supplied: {}", n),
        }

        let (sender, receiver) = crossbeam::channel::bounded(10);
        let mut block_count = 0;
        let join_handle = std::thread::spawn(move || {
            let (right_n, right_readers) = init_bsx_readers(
                readers_pair
                    .pop()
                    .unwrap()
                    .into_iter()
                    .map(|(_, r)| r)
                    .collect_vec(),
            );
            let (left_n, left_readers) = init_bsx_readers(
                readers_pair
                    .pop()
                    .unwrap()
                    .into_iter()
                    .map(|(_, r)| r)
                    .collect_vec(),
            );
            assert_eq!(left_n, right_n, "Number of batches does not match");
            block_count += left_n;
            segment_reading(left_readers, right_readers, sender)
        });

        let out = DmrIterator {
            config: self.clone(),
            ref_idset: RefIDSet::new(),
            leftover: None,
            regions_cache: Vec::new(),
            last_chr: Arc::new(String::new()),
            receiver,
            reader_stat: ReaderMetadata::new(block_count),
            _join_handle: join_handle,
        };

        Ok(out)
    }

    pub fn new(
        context: Context,
        n_missing: usize,
        min_coverage: i16,
        diff_threshold: f32,
        min_cpgs: usize,
        max_dist: u64,
        initial_l: f64,
        l_min: f64,
        l_coef: f64,
        seg_tolerance: f32,
        merge_pvalue: f64,
        seg_pvalue: f64,
    ) -> Self {
        Self {
            context,
            n_missing,
            min_coverage,
            diff_threshold,
            min_cpgs,
            max_dist,
            initial_l,
            l_min,
            l_coef,
            seg_tolerance,
            merge_pvalue,
            seg_pvalue,
        }
    }
}

impl Default for DmrConfig {
    fn default() -> Self {
        Self {
            context: Context::CG,
            n_missing: 0,
            min_coverage: 5,
            diff_threshold: 0.1,
            min_cpgs: 10,
            max_dist: 100,
            initial_l: 2.0,
            l_min: 1e-3,
            l_coef: 1.5,
            seg_tolerance: 1e-6,
            merge_pvalue: 1e-3,
            seg_pvalue: 1e-2,
        }
    }
}

fn init_bsx_readers<F: Read + Seek + 'static>(
    handles: Vec<F>
) -> (usize, Vec<Box<dyn Iterator<Item = EncodedBsxBatch>>>) {
    assert!(!handles.is_empty(), "No readers supplied");
    let bsx_readers = handles
        .into_iter()
        .map(|reader| BsxFileReader::new(reader))
        .collect_vec();
    assert!(
        bsx_readers
            .iter()
            .map(|r| r.blocks_total())
            .all_equal(),
        "Number of blocks not equal in files"
    );
    let n_batches = bsx_readers[0].blocks_total();
    let iterators = bsx_readers
        .into_iter()
        .map(|reader| {
            reader.map(|batch_res| batch_res.expect("could not read batch"))
        })
        .map(|reader| {
            Box::new(reader) as Box<dyn Iterator<Item = EncodedBsxBatch>>
        })
        .collect_vec();
    (n_batches, iterators)
}
