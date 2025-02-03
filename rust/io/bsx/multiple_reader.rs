use crate::bsx_batch::EncodedBsxBatch;
use crate::io::bsx::region_read::BsxFileReader;
use itertools::Itertools;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::error::Error;
use std::fmt::Debug;
use std::hash::Hash;
use std::io;
use std::io::{Read, Seek};
use std::sync::{Arc, Mutex, RwLock};

pub struct MultiBsxFileReader<M, R>
where
    M: Hash + Serialize + Eq + Send + Sync + Debug + 'static,
    R: Read + Seek + Send + Sync + 'static,
{
    current_batch_idx: usize,
    readers: HashMap<Arc<M>, Arc<RwLock<BsxFileReader<R>>>>,
    _keys_bound: Vec<Arc<M>>,
}

impl<M, R> MultiBsxFileReader<M, R>
where
    M: Hash + Serialize + Eq + Send + Sync + Debug + 'static,
    R: Read + Seek + Send + Sync + 'static,
{
    pub fn try_new(readers: HashMap<M, BsxFileReader<R>>) -> Result<Self, Box<dyn Error>> {
        if readers.is_empty() {
            return Err(io::Error::from(io::ErrorKind::InvalidInput).into());
        }
        if !(readers.values().map(|r| r.blocks_total()).all_equal()) {
            return Err("Different number of batches in files".into());
        }

        let readers: HashMap<Arc<M>, Arc<RwLock<BsxFileReader<R>>>> = HashMap::from_iter(
            readers
                .into_iter()
                .map(|(key, reader)| (Arc::new(key), Arc::new(RwLock::new(reader)))),
        );

        let _keys_bound = readers.keys().cloned().collect_vec();

        Ok(Self {
            _keys_bound,
            readers,
            current_batch_idx: 0,
        })
    }

    fn set_batch(&mut self, batch_idx: usize) {
        self.current_batch_idx = batch_idx;
    }

    pub fn get_batch(
        &mut self,
        batch_idx: usize,
    ) -> Result<Option<Vec<(Arc<M>, EncodedBsxBatch)>>, Box<dyn Error>> {
        self.set_batch(batch_idx);

        let batch_results = Arc::new(Mutex::new(Vec::new())); // Thread-safe collection

        rayon::scope(|s| {
            for (thread_id, reader) in self.readers.iter() {
                let thread_id = Arc::clone(thread_id);
                let reader = Arc::clone(reader);
                let batch_results = Arc::clone(&batch_results);

                s.spawn(move |_| {
                    let batch = reader.write().unwrap().get_batch(batch_idx);
                    let mut results = batch_results.lock().unwrap();
                    results.push((thread_id, batch));
                });
            }
        });

        let batch_results = Arc::try_unwrap(batch_results)
            .expect("Arc::try_unwrap failed")
            .into_inner()
            .expect("Mutex was poisoned")
            .into_iter()
            .filter(|(_, batch)| batch.is_some())
            .map(|(key, batch)| (key, batch.unwrap()))
            .collect_vec();

        if batch_results.is_empty() {
            return Ok(None);
        }
        if batch_results.len() != self.readers.len() {
            return Err(format!("Batches missing in some readers {}", batch_idx).into());
        }

        for (_, batch_result) in batch_results.iter() {
            if let Err(e) = batch_result {
                return Err(e.to_string().into());
            }
        }

        Ok(Some(
            batch_results
                .into_iter()
                .map(|(k, v)| (k, v.unwrap()))
                .collect(),
        ))
    }
}

impl<M, R> Iterator for MultiBsxFileReader<M, R>
where
    M: Hash + Serialize + Eq + Send + Sync + Debug + 'static,
    R: Read + Seek + Send + Sync + 'static,
{
    type Item = Vec<(Arc<M>, EncodedBsxBatch)>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.get_batch(self.current_batch_idx);
        self.current_batch_idx += 1;
        result.expect("Error while reading batches")
    }
}
