use std::collections::{BTreeMap, HashMap};
use std::fmt::Display;
use std::hash::Hash;
use std::io::{Read, Seek};
use std::sync::Arc;

use anyhow::anyhow;
use bio_types::annot::refids::RefIDSet;
use itertools::Itertools;

use crate::data_structs::bsx_batch_group::EncodedBsxBatchGroup;
use crate::io::bsx::multiple_reader::MultiBsxFileReader;
use crate::io::bsx::read::BsxFileReader;
use crate::tools::dmr::data_structs::ReaderMetadata;
use crate::tools::dmr::segmentation::FilterConfig;
use crate::tools::dmr::DmrIterator;
use crate::utils::types::Context;

#[derive(Debug, Clone)]
pub struct DmrConfig {
    pub context:        Context,
    pub n_missing:      usize,
    pub min_coverage:   i16,
    pub diff_threshold: f64,
    pub min_cpgs:       usize,
    pub max_dist:       u64,
    pub initial_l:      f64,
    pub l_min:          f64,
    pub l_coef:         f64,
    pub seg_tolerance:  f64,
    pub merge_pvalue:   f64,
    pub seg_pvalue:     f64,
}

impl DmrConfig {
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
    ) -> anyhow::Result<DmrIterator<R>>
    where
        F: Read + Seek + Send + Sync + 'static,
        R: Display
            + Eq
            + Hash
            + Clone
            + Default
            + std::fmt::Debug
            + Send
            + 'static, {
        let sample_mapping: HashMap<uuid::Uuid, (R, F)> = HashMap::from_iter(
            readers
                .into_iter()
                .map(|(group, handle)| (uuid::Uuid::new_v4(), (group, handle))),
        );

        let group_mapping: HashMap<uuid::Uuid, R> = HashMap::from_iter(
            sample_mapping
                .iter()
                .map(|(id, (group, _))| (id.clone(), group.clone())),
        );
        let group_order = {
            let group_set = BTreeMap::from_iter(
                group_mapping
                    .values()
                    .map(|x| (x.to_string(), x)),
            );
            if group_set.len() != 2 {
                return Err(anyhow!(
                    "There should be only two groups of samples! ({:?})",
                    group_set
                ));
            }
            let groups_vec = group_set.into_iter().collect_vec();
            (
                groups_vec[0].clone().1.clone(),
                groups_vec[1].clone().1.clone(),
            )
        };
        let readers_mapping: HashMap<uuid::Uuid, BsxFileReader<F>> =
            HashMap::from_iter(
                sample_mapping
                    .into_iter()
                    .map(|(id, (_, handle))| (id, BsxFileReader::new(handle))),
            );

        let config_copy = self.clone();

        let (sender, receiver) = crossbeam::channel::bounded(10);
        let last_chr = Arc::new(String::new());
        let ref_idset = RefIDSet::new();

        let group_mapping_clone = group_mapping.clone();
        let mut multi_reader = MultiBsxFileReader::try_new(readers_mapping)
            .map_err(|e| anyhow!("{:?}", e))?;
        let reader_stat = ReaderMetadata::new(multi_reader.blocks_total());

        let join_handle = std::thread::spawn(move || {
            let group_mapping = group_mapping_clone;

            while let Some(batches) = multi_reader.next() {
                let labels = batches
                    .iter()
                    .map(|(id, _)| group_mapping.get(id).unwrap().clone())
                    .collect();
                let data = batches
                    .into_iter()
                    .map(|(_, batch)| batch)
                    .collect();
                let group =
                    EncodedBsxBatchGroup::try_new(data, Some(labels)).unwrap();
                sender.send(group).unwrap();
            }
        });

        let out = DmrIterator {
            config: config_copy,
            ref_idset,
            group_pair: group_order,
            leftover: None,
            regions_cache: Vec::new(),
            last_chr,
            receiver,
            reader_stat,
            _join_handle: join_handle,
        };

        Ok(out)
    }

    pub fn new(
        context: Context,
        n_missing: usize,
        min_coverage: i16,
        diff_threshold: f64,
        min_cpgs: usize,
        max_dist: u64,
        initial_l: f64,
        l_min: f64,
        l_coef: f64,
        seg_tolerance: f64,
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
