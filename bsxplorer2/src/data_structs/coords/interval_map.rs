#![allow(unused)]
use hashbrown::HashMap;
use itertools::Itertools;
use rust_lapper::{
    Interval,
    Lapper,
};
use serde::de::DeserializeOwned;
use serde::{
    Deserialize,
    Serialize,
};

use super::Contig;
use crate::data_structs::typedef::PosType;
use crate::utils::THREAD_POOL;
use crate::BsxSmallStr;

#[derive(Clone, Serialize, Deserialize, Debug)]
#[serde(bound = "V: Serialize + DeserializeOwned")]
pub struct ContigIntervalMap<V>
where
    V: Sync + Send + Eq + Clone, {
    inner: HashMap<BsxSmallStr, Lapper<PosType, V>>,
}

impl<V> Default for ContigIntervalMap<V>
where
    V: Sync + Send + Eq + Clone,
{
    fn default() -> Self {
        Self {
            inner: HashMap::new(),
        }
    }
}

impl<V> FromIterator<(Contig, V)> for ContigIntervalMap<V>
where
    V: Sync + Send + Eq + Clone,
{
    fn from_iter<T: IntoIterator<Item = (Contig, V)>>(iter: T) -> Self {
        let multimap = iter
            .into_iter()
            .map(|(c, v)| {
                let key = c.seqname().clone();
                (key, (c, v))
            })
            .into_group_map();

        let mut inner = HashMap::with_capacity(multimap.len());
        for (chr, kv_pairs) in multimap.into_iter() {
            let imap = Lapper::new(
                kv_pairs
                    .into_iter()
                    .map(|(c, v)| {
                        Interval {
                            start: c.start(),
                            stop:  c.end(),
                            val:   v,
                        }
                    })
                    .collect_vec(),
            );
            inner.insert(chr, imap);
        }

        Self { inner }
    }
}

impl<V> ContigIntervalMap<V>
where
    V: Sync + Send + Eq + Clone,
{
    pub fn new() -> Self {
        Self::default()
    }

    pub fn into_inner(self) -> HashMap<BsxSmallStr, Lapper<PosType, V>> {
        self.inner
    }

    pub fn inner(&self) -> &HashMap<BsxSmallStr, Lapper<PosType, V>> {
        &self.inner
    }

    pub fn n_chr(&self) -> usize {
        self.inner().len()
    }

    pub fn chr_names(&self) -> Vec<BsxSmallStr> {
        self.inner().keys().cloned().collect()
    }

    pub fn insert(
        &mut self,
        key: Contig,
        value: V,
    ) {
        let imap = self
            .inner
            .entry(key.seqname().clone())
            .or_insert_with(|| Lapper::new(vec![]));
        imap.insert(Interval {
            start: key.start(),
            stop:  key.end(),
            val:   value,
        });
    }

    pub fn find(
        &self,
        key: &Contig,
    ) -> Option<Vec<&V>> {
        self.inner.get(key.seqname()).map(|imap| {
            imap.find(key.start(), key.end())
                .map(|e| (&e.val))
                .collect_vec()
        })
    }

    pub fn union(
        &mut self,
        other: &Self,
    ) {
        for (chr, imap) in other.inner.iter() {
            self.inner
                .entry(chr.clone())
                .and_modify(|old| {
                    old.union(imap);
                })
                .or_insert(imap.clone());
        }
    }

    pub fn get_breakpoints(&self) -> Vec<(BsxSmallStr, Vec<u32>)> {
        self.inner()
            .iter()
            .map(|(chr, lap)| {
                (
                    chr.clone(),
                    lap.iter().map(|i| vec![i.start, i.stop]).concat(),
                )
            })
            .collect_vec()
    }
}
