use log::info;
use polars::prelude::*;
use polars_arrow::io::ipc::read::FileMetadata;
use std::collections::btree_map::Cursor;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{Seek, SeekFrom};
use std::ops::Bound::Included;

use crate::io::bsx::ipc::IpcFileReader;
use crate::ubatch::UBatch;

type BSXIndex = HashMap<String, BTreeMap<u32, usize>>;

/// Structure, that holds chromosome-wise index of batches from BSXplorer
/// IPC file.
#[derive(Clone)]
pub struct BSXBTree(BSXIndex);

impl BSXBTree {
    pub(crate) fn new() -> Self {
        Self(HashMap::new())
    }

    /// Insert node to a [BTreeMap].
    pub(crate) fn insert(&mut self, chr: String, start: u32, idx: usize) {
        match self.0.get_mut(&chr) {
            Some(btree) => {
                btree.insert(start, idx);
            }
            None => {
                let mut btree = BTreeMap::new();
                btree.insert(start, idx);
                self.0.insert(chr, btree);
            }
        };
    }
    /// Get batch, which starts exactly on specified position.
    pub fn get_exact(&self, chr: String, start: u32) -> Option<&usize> {
        self.0.get(&chr)?.get(&start)
    }
    /// Get cursor to the first batch, which has lower start position.
    pub fn get_lower_bound(&self, chr: String, start: u32) -> Option<Cursor<'_, u32, usize>> {
        Some(self.0.get(&chr)?.lower_bound(Included(&start)))
    }
    /// Get all batch indexes, which contain data about specified region
    pub fn get_region(&self, chr: &str, start: u32, end: u32) -> Option<Vec<usize>> {
        let mut lower_bound = self.get_lower_bound(chr.to_string(), start)?;
        let mut batches = vec![lower_bound.prev()];
        while let Some((start_val, index)) = lower_bound.next() {
            if start_val > &end {
                break;
            } else {
                batches.push(Some((start_val, index)));
            }
        }
        Some(
            batches
                .iter()
                .flatten()
                .map(|(_key, index)| (*index).clone())
                .collect(),
        )
    }
}

/// Low-level reader for BSXplorer IPC file (made by [crate::io::bsx::write::BSXWriter]).
pub struct BSXReader {
    reader: IpcFileReader<File>,
    metadata: FileMetadata,
    pub index: BSXBTree,
}

impl BSXReader {
    pub fn new(mut handle: File, projection: Option<Vec<usize>>) -> Self {
        let temp_handle = handle.try_clone().unwrap();
        let mut temp_reader = IpcFileReader::new(temp_handle, projection.clone(), None);
        info!("Opened BSX file");
        let mut index = BSXBTree::new();

        // Create index
        for (num, df) in {
            (0..temp_reader.blocks_total())
                .map(|i| temp_reader.read_df_at(i))
                .enumerate()
        } {
            let df = df.expect("failed to read df");

            let start = df
                .column("position")
                .unwrap()
                .u32()
                .unwrap()
                .first()
                .unwrap();

            let chr_col = df.column("chr").unwrap().categorical().unwrap();
            let chr_mapping = chr_col.get_rev_map();

            let chr = chr_mapping
                .get(chr_col.physical().get(0).unwrap())
                .to_owned();

            index.insert(chr, start, num);
        }
        info!(
            "Indexed BSX file with {} batches",
            temp_reader.blocks_total()
        );

        // Reinitialize reader
        handle.seek(SeekFrom::Start(0)).unwrap();
        let reader = IpcFileReader::new(handle, projection, None);
        let metadata = reader.metadata().clone();

        Self {
            reader,
            metadata,
            index,
        }
    }
}

impl Iterator for BSXReader {
    type Item = UBatch;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.next() {
            Some(record_batch) => {
                let record_batch = record_batch.expect("Could not create record batch");
                let df = DataFrame::try_from((record_batch, self.metadata.schema.as_ref()))
                    .expect("Could not convert record batch into DataFrame");
                Some(df.into())
            }
            None => None,
        }
    }

    fn count(self) -> usize
    where
        Self: Sized,
    {
        self.reader.blocks_total()
    }

    fn last(mut self) -> Option<Self::Item>
    where
        Self: Sized,
    {
        self.nth(self.reader.blocks_total() - 1)
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        Some(self.reader.read_df_at(n).expect("Could not read df").into())
    }
}
