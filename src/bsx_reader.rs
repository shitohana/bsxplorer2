use std::collections::{BTreeMap, HashMap};
use std::collections::Bound::Excluded;
use std::io::{Read, Seek};
use polars::prelude::*;
use crate::utils::ipc_reader::IpcFileReader;
use crate::utils::types::{Context, Strand};

#[derive(Eq, Hash, PartialEq, Clone)]
pub struct BSXKey(String, Strand, Context);

pub type BSXIndex = HashMap<BSXKey, BTreeMap<u64, usize>>;

struct BSXCache(HashMap<u32, CacheItem>);
struct CacheItem {
    calls: u32,
    df: Arc<DataFrame>,
}

impl BSXCache {
    fn new(capacity: usize) -> Self {
        Self(HashMap::with_capacity(capacity))
    }
    fn insert(&mut self, key: u32, frame: Arc<DataFrame>) -> Result<(), PolarsError> {
        if self.0.len() < self.0.capacity() {
            self.0.insert(key, CacheItem{calls: 0, df: frame});
            Ok(())
        } else {
            let min_key = self.0.iter().min_by(|x, y| {x.1.calls.cmp(&y.1.calls)}).unwrap().0.clone();
            self.0.remove(&min_key);
            self.0.insert(key, CacheItem{calls: 0, df: frame});
            Ok(())
        }
    }
    
    fn get(&mut self, key: u32) -> Option<Arc<DataFrame>> {
        if let Some(item) = self.0.get_mut(&key) {
            (*item).calls += 1;
            Some(item.df.clone())
        } else { None }
    }
}

struct BSXReader<R: Read + Seek + Sync> {
    handle: IpcFileReader<R>,
    index: BSXIndex,
    cache: BSXCache,
}

// Constructors
impl<R: Read + Seek + Sync> BSXReader<R> {
    pub fn new(mut handle: IpcFileReader<R>, cache_capacity: usize) -> Result<Self, PolarsError> {
        let index = {
            let mut index: HashMap<BSXKey, BTreeMap<u64, usize>> = HashMap::new();
            for batch_num in 0..handle.metadata().blocks.len() {
                match handle.read_df_at(batch_num) {
                    Ok(df) => {
                        // Read first position of the batch
                        let pos = df
                            .column("position")?
                            .u64()?
                            .get(0)
                            .expect("could not get position");
                        // Create key instance for hashmap
                        let key: BSXKey = {
                            let chr_col = df.column("chr")?.categorical()?;
                            BSXKey(
                                chr_col
                                    .get_rev_map()
                                    .get(chr_col.physical().get(0).unwrap())
                                    .to_owned(),
                                Strand::from_bool(df.column("strand")?.bool()?.get(0)),
                                Context::from_bool(df.column("context")?.bool()?.get(0)),
                            )
                        };
                        // If key exists, append batch position
                        match index.get_mut(&key) {
                            Some(btree) => {
                                btree.insert(pos, batch_num);
                            }
                            None => {
                                let mut btree: BTreeMap<u64, usize> = BTreeMap::new();
                                btree.insert(pos, batch_num);
                                index.insert(key, btree);
                            }
                        };
                    }
                    Err(e) => panic!("error reading batch {}: {}", batch_num, e),
                };
            }
            index
        };
        let cache = BSXCache::new(cache_capacity);
        Ok(Self { handle, index, cache })
    }
}

// Getters
impl<R: Read + Seek + Sync> BSXReader<R> {
    pub fn get_batch_num(&mut self) -> usize {
        self.handle.metadata().blocks.len()
    }
    pub fn get_arrow_schema(&self) -> ArrowSchema {
        self.handle.schema().clone()
    }
    pub fn get_polars_schema(&self) -> Schema {
        let schema = self.get_arrow_schema();
        Schema::from_arrow_schema(&schema)
    }
}

// Reading
impl<R: Read + Seek + Sync> BSXReader<R> {
    pub fn get_batch(&mut self, n: u32) -> Arc<DataFrame> {
        if let Some(df) = self.cache.get(n) {
            df.clone()
        } else {
            let df = self.handle.read_df_at(n as usize).expect("could not get df");
            self.cache.insert(n, Arc::new(df)).expect("Failed to update cache");
            self.cache.get(n).expect("Failed to update cache").clone()
        }
    }
    
    pub fn get_index(&self, key: BSXKey, range: (Option<u64>, Option<u64>)) -> Option<Vec<u32>> {
        match self.index.get(&key) {
            Some(btree) => {
                // Init range
                let start = range.0
                    .unwrap_or(0);
                let end = range.1
                    .unwrap_or(*btree.iter().max().unwrap().0);
                // Get tree root
                let mut root = btree.lower_bound(Excluded(&start));
                // Create output vec
                let mut out: Vec<u32> = Vec::new();
                // If prev exists, append to out
                if let Some(prev) = root.peek_prev() {
                    out.push(*prev.1 as u32);
                }
                // While bounds satisfied, append to out
                while let Some(new_batch) = root.next() {
                    if *new_batch.0 <= end {
                        out.push(*new_batch.1 as u32);
                    }
                    else {
                        break;
                    };
                }
                Some(out)
            },
            None => None
        }
    }
    
    pub fn get_region(
        &mut self,
        indices: &[u32],
        start: u64,
        end: u64
    ) -> Result<DataFrame, PolarsError> {
        let dfs: Vec<LazyFrame> = indices
            .iter()
            .map(|i| self.get_batch(i.clone()))
            .map(|df| {
                df.as_ref().clone().lazy().filter(
                    col("position")
                        .gt_eq(lit(start))
                        .and(col("position").lt_eq(lit(end))),
                )
            })
            .collect();
        concat(dfs, UnionArgs::default())?.collect()
    }
}