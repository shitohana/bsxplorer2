use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use std::io;
use std::io::{Read, Seek};
use std::sync::{Arc, Mutex, RwLock};

use anyhow::{Context, Result};
use itertools::Itertools;
use log::{debug, info, warn};
use serde::Serialize;

use crate::data_structs::batch::BsxBatchMethods;
use crate::data_structs::batch::EncodedBsxBatch;
use crate::io::bsx::read::BsxFileReader;

/// A reader for multiple BSX files that handles concurrent batch reading
/// across multiple files with synchronized batch indexing.
pub struct MultiBsxFileReader<M, R>
where
    M: Hash + Serialize + Eq + Send + Sync + Debug + 'static,
    R: Read + Seek + Send + Sync + 'static,
{
    /// Current batch index position in the readers
    current_batch_idx: usize,
    /// Map of metadata keys to their corresponding file readers
    readers: HashMap<Arc<M>, Arc<RwLock<BsxFileReader<R>>>>,
    /// Cached vector of keys for fast iteration
    _keys_bound: Vec<Arc<M>>,
}

impl<M, R> MultiBsxFileReader<M, R>
where
    M: Hash + Serialize + Eq + Send + Sync + Debug + 'static,
    R: Read + Seek + Send + Sync + 'static,
{
    /// Creates a new MultiBsxFileReader from a map of metadata to readers.
    ///
    /// # Arguments
    /// * `readers` - HashMap mapping metadata to BsxFileReader instances
    ///
    /// # Returns
    /// * `Result<Self>` - A new MultiBsxFileReader or an error
    ///
    /// # Errors
    /// * If the readers map is empty
    /// * If the readers have different numbers of batches
    pub fn try_new(readers: HashMap<M, BsxFileReader<R>>) -> Result<Self> {
        if readers.is_empty() {
            return Err(io::Error::from(io::ErrorKind::InvalidInput)).context(
                "Cannot create MultiBsxFileReader with empty readers map",
            );
        }

        let batch_counts = readers
            .values()
            .map(|r| r.blocks_total())
            .collect_vec();
        debug!("Reader batch counts: {:?}", batch_counts);

        if !batch_counts.iter().all_equal() {
            return Err(anyhow::anyhow!(
                "Files have inconsistent batch counts: {:?}",
                batch_counts
            ));
        }

        info!(
            "Creating MultiBsxFileReader with {} readers, {} batches each",
            readers.len(),
            batch_counts[0]
        );

        let readers: HashMap<Arc<M>, Arc<RwLock<BsxFileReader<R>>>> =
            HashMap::from_iter(
                readers
                    .into_iter()
                    .map(|(key, reader)| {
                        (Arc::new(key), Arc::new(RwLock::new(reader)))
                    }),
            );

        let _keys_bound = readers.keys().cloned().collect_vec();

        Ok(Self {
            _keys_bound,
            readers,
            current_batch_idx: 0,
        })
    }

    /// Sets the current batch index
    ///
    /// # Arguments
    /// * `batch_idx` - The batch index to set
    fn set_batch(
        &mut self,
        batch_idx: usize,
    ) {
        debug!("Setting batch index to {}", batch_idx);
        self.current_batch_idx = batch_idx;
    }

    /// Returns the current batch index
    ///
    /// # Returns
    /// * `usize` - Current batch index
    pub fn current_batch_idx(&self) -> usize {
        self.current_batch_idx
    }

    /// Retrieves a batch from all readers at the specified index
    ///
    /// # Arguments
    /// * `batch_idx` - The batch index to retrieve
    ///
    /// # Returns
    /// * `Result<Option<Vec<(Arc<M>, EncodedBsxBatch)>>>` - Vector of
    ///   (metadata, batch) pairs or None if no batches remain
    ///
    /// # Errors
    /// * If batches are missing from some readers
    /// * If batch retrieval fails for any reader
    pub fn get_batch(
        &mut self,
        batch_idx: usize,
    ) -> Result<Option<Vec<(Arc<M>, EncodedBsxBatch)>>> {
        self.set_batch(batch_idx);
        info!(
            "Getting batch {} from {} readers",
            batch_idx,
            self.readers.len()
        );

        let batch_results = Arc::new(Mutex::new(Vec::new())); // Thread-safe collection

        rayon::scope(|s| {
            for (thread_id, reader) in self.readers.iter() {
                let thread_id = Arc::clone(thread_id);
                let reader = Arc::clone(reader);
                let batch_results = Arc::clone(&batch_results);

                s.spawn(move |_| {
                    debug!("Spawned thread for reader {:?}", thread_id);
                    let batch = reader
                        .write()
                        .unwrap()
                        .get_batch(batch_idx);
                    let mut results = batch_results.lock().unwrap();
                    debug!("Thread for reader {:?} completed", thread_id);
                    results.push((thread_id, batch));
                });
            }
        });

        let batch_results = Arc::try_unwrap(batch_results)
            .expect("Arc::try_unwrap failed - this is a bug")
            .into_inner()
            .expect(
                "Mutex was poisoned - this indicates a panic in a reader \
                 thread",
            )
            .into_iter()
            .filter(|(_, batch)| batch.is_some())
            .map(|(key, batch)| (key, batch.unwrap()))
            .filter(|(_, batch)| {
                batch
                    .as_ref()
                    .map(|b| b.height() > 0)
                    .unwrap_or(false)
            })
            .collect_vec();

        if batch_results.is_empty() {
            debug!("No valid batches found at index {}", batch_idx);
            return Ok(None);
        }

        if batch_results.len() != self.readers.len() {
            warn!(
                "Batch {} missing from some readers. Got {} of {} expected \
                 batches",
                batch_idx,
                batch_results.len(),
                self.readers.len()
            );
            return Err(anyhow::anyhow!(
                "Batch {} missing from some readers ({}/{} available)",
                batch_idx,
                batch_results.len(),
                self.readers.len()
            ));
        }

        // Check for errors in any of the batches
        for (reader_id, batch_result) in batch_results.iter() {
            if let Err(e) = batch_result {
                return Err(anyhow::anyhow!(
                    "Error reading batch {} from reader {:?}: {}",
                    batch_idx,
                    reader_id,
                    e
                ));
            }
        }

        let result = batch_results
            .into_iter()
            .map(|(k, v)| (k, v.unwrap()))
            .collect();

        debug!("Successfully retrieved batch {}", batch_idx);
        Ok(Some(result))
    }

    /// Returns the total number of blocks/batches in the files
    ///
    /// # Returns
    /// * `usize` - Number of blocks in the files
    pub fn blocks_total(&self) -> usize {
        let total = self
            .readers
            .values()
            .next()
            .expect("Reader collection is empty - this is a bug")
            .read()
            .unwrap()
            .blocks_total();

        debug!("Total blocks: {}", total);
        total
    }
}

impl<M, R> Iterator for MultiBsxFileReader<M, R>
where
    M: Hash + Serialize + Eq + Send + Sync + Debug + 'static,
    R: Read + Seek + Send + Sync + 'static,
{
    type Item = Vec<(Arc<M>, EncodedBsxBatch)>;

    /// Advances to the next batch in all readers
    ///
    /// # Returns
    /// * `Option<Self::Item>` - The next batch or None if no batches remain
    fn next(&mut self) -> Option<Self::Item> {
        debug!("Iterator advancing to batch {}", self.current_batch_idx);
        let result = self
            .get_batch(self.current_batch_idx)
            .context(format!(
                "Failed to read batch at index {}",
                self.current_batch_idx
            ))
            .unwrap(); // Propagate errors in iterator
        self.current_batch_idx += 1;
        result
    }
}
