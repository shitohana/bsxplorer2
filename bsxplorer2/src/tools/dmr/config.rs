use std::collections::BTreeMap;
use std::fmt::Display;
use std::hash::Hash;
use std::os::unix::prelude::AsRawFd;
use std::sync::atomic::AtomicUsize;
use std::sync::Arc;

use bio::bio_types::annot::refids::RefIDSet;
use itertools::Itertools;

use crate::data_structs::typedef::*;
use crate::{prelude::*, with_field_fn};
use crate::tools::dmr::data_structs::ReaderMetadata;
use crate::tools::dmr::segmentation::FilterConfig;
use crate::tools::dmr::{
    segment_reading,
    DmrIterator,
};

#[derive(Debug, Clone)]
pub struct DmrConfig {
    pub context:        Context,
    pub n_missing:      usize,
    pub min_coverage:   CountType,
    pub diff_threshold: DensityType,
    pub min_cpgs:       usize,
    pub max_dist:       PosType,
    pub initial_l:      f64,
    pub l_min:          f64,
    pub l_coef:         f64,
    pub seg_tolerance:  DensityType,
    pub merge_pvalue:   f64,
    pub seg_pvalue:     f64,
}

impl DmrConfig {
    with_field_fn!(context, Context);
    with_field_fn!(n_missing, usize);
    with_field_fn!(min_coverage, CountType);
    with_field_fn!(diff_threshold, DensityType);
    with_field_fn!(min_cpgs, usize);
    with_field_fn!(max_dist, PosType);
    with_field_fn!(initial_l, f64);
    with_field_fn!(l_min, f64);
    with_field_fn!(l_coef, f64);
    with_field_fn!(seg_tolerance, DensityType);
    with_field_fn!(merge_pvalue, f64);
    with_field_fn!(seg_pvalue, f64);

    pub fn filter_config(&self) -> FilterConfig {
        FilterConfig {
            context:      self.context,
            n_missing:    self.n_missing,
            min_coverage: self.min_coverage,
        }
    }

    pub fn try_finish<F, R>(
        &self,
        readers: Vec<(R, F)>,
    ) -> anyhow::Result<DmrIterator>
    where
        F: AsRawFd + 'static + Sync + Send,
        R: Display
            + Eq
            + Hash
            + Clone
            + Default
            + std::fmt::Debug
            + Send
            + Ord
            + 'static, {
        // Use BTreeMap to sort labels and get consistent left/right order
        let mut readers_pair = BTreeMap::from_iter(
            readers
                .into_iter()
                .into_group_map_by(|(label, _reader)| label.clone()),
        )
        .into_values()
        .collect_vec();

        match readers_pair.len() {
            0 => anyhow::bail!("No readers supplied"),
            2 => {},
            n => anyhow::bail!("Too many readers groups supplied: {}", n),
        }

        let (sender, receiver) = crossbeam::channel::bounded(10);
        let block_count = Arc::new(AtomicUsize::new(0));
        let local_block_count = block_count.clone();

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
            local_block_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            segment_reading(left_readers, right_readers, sender)
        });

        let out = DmrIterator {
            config: self.clone(),
            ref_idset: RefIDSet::new(),
            leftover: None,
            regions_cache: Vec::new(),
            last_chr: Arc::new(String::new()),
            receiver,
            reader_stat: ReaderMetadata::new(
                block_count.load(std::sync::atomic::Ordering::Relaxed),
            ),
            _join_handle: join_handle,
        };

        Ok(out)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new(
        context: Context,
        n_missing: usize,
        min_coverage: CountType,
        diff_threshold: DensityType,
        min_cpgs: usize,
        max_dist: PosType,
        initial_l: f64,
        l_min: f64,
        l_coef: f64,
        seg_tolerance: DensityType,
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
            context:        Context::CG,
            n_missing:      0,
            min_coverage:   5,
            diff_threshold: 0.1,
            min_cpgs:       10,
            max_dist:       100,
            initial_l:      2.0,
            l_min:          1e-3,
            l_coef:         1.5,
            seg_tolerance:  1e-6,
            merge_pvalue:   1e-3,
            seg_pvalue:     1e-2,
        }
    }
}

fn init_bsx_readers<F: AsRawFd + 'static>(
    handles: Vec<F>
) -> (usize, Vec<Box<dyn Iterator<Item = BsxBatch>>>) {
    assert!(!handles.is_empty(), "No readers supplied");
    let bsx_readers = handles
        .into_iter()
        .map(|reader| BsxFileReader::try_new(reader))
        .collect::<anyhow::Result<Vec<_>>>()
        .expect("Failed to initialize readers");
    assert!(
        bsx_readers.iter().map(|r| r.blocks_total()).all_equal(),
        "Number of blocks not equal in files"
    );
    #[allow(unsafe_code)]
    let n_batches = unsafe { bsx_readers.get_unchecked(0).blocks_total() };
    let iterators = bsx_readers
        .into_iter()
        .map(|reader| {
            Box::new(
                reader
                    .into_iter()
                    .map(|batch_res| batch_res.expect("could not read batch")),
            ) as Box<dyn Iterator<Item = BsxBatch>>
        })
        .collect_vec();
    (n_batches, iterators)
}
