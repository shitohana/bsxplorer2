use hashbrown::HashMap;
use itertools::Itertools;
use multimap::MultiMap;

use crate::data_structs::annotation::{GffEntry, GffEntryAttributes};
use crate::data_structs::coords::Contig;
use crate::data_structs::typedef::{BsxSmallStr, SeqNameStr, SeqPosNum};
use crate::io::bsx::BatchIndex;

const UPSTREAM_TYPE_NAME: &str = "upstream";
const DOWNSTREAM_TYPE_NAME: &str = "downstream";

pub struct AnnotStore {
    id_map:       HashMap<BsxSmallStr, GffEntry>,
    parent_map:   MultiMap<BsxSmallStr, BsxSmallStr>,
    children_map: MultiMap<BsxSmallStr, BsxSmallStr>,
}

impl Default for AnnotStore {
    fn default() -> Self { Self::new() }
}

impl AnnotStore {
    pub fn new() -> Self {
        Self {
            id_map:       Default::default(),
            parent_map:   Default::default(),
            children_map: Default::default(),
        }
    }

    pub fn insert(
        &mut self,
        entry: GffEntry,
    ) -> Option<()> {
        if self.id_map.contains_key(&entry.id) {
            None
        }
        else {
            self.id_map
                .insert(entry.id.clone(), entry.clone());
            for parent in entry
                .attributes
                .parent
                .as_ref()
                .unwrap_or(&vec![])
            {
                self.parent_map
                    .insert(entry.id.clone(), parent.to_owned());
                self.children_map
                    .insert(parent.to_owned(), entry.id.clone());
            }
            Some(())
        }
    }

    pub fn remove(
        &mut self,
        id: &BsxSmallStr,
    ) -> Option<GffEntry> {
        let removed = self.id_map.remove(id);
        let _removed_parents = self.parent_map.remove(id);
        let _removed_children = self.children_map.remove(id);

        removed
    }

    pub fn iter(&self) -> AnnotIterator {
        AnnotIterator {
            store: self,
            iter:  Box::from(self.id_map.clone().into_iter()),
        }
    }

    pub fn iter_sorted<S: SeqNameStr, P: SeqPosNum>(
        &self,
        index: &BatchIndex<S, P>,
    ) -> AnnotIterator {
        let iterator: Box<dyn Iterator<Item = _>> =
            Box::from(
                self.iter()
                    .into_group_map_by(|(_id, entry)| {
                        index.get_chr_index(&S::from(
                            entry.contig.seqname().as_str(),
                        ))
                    })
                    .into_iter()
                    .sorted_by_cached_key(|(key, _value)| *key)
                    .map(|(_key, value)| {
                        value.into_iter().sorted_by_cached_key(
                            |(_id, entry)| entry.contig.start(),
                        )
                    })
                    .flatten()
                    .map(|(id, value)| (id.clone(), value.clone())),
            );

        AnnotIterator {
            store: self,
            iter:  iterator,
        }
    }

    pub fn get_parents(
        &self,
        id: &BsxSmallStr,
    ) -> Option<&Vec<BsxSmallStr>> {
        self.parent_map.get_vec(id)
    }

    pub fn get_children(
        &self,
        id: &BsxSmallStr,
    ) -> Option<&Vec<BsxSmallStr>> {
        self.children_map.get_vec(id)
    }

    pub fn add_upstream<F: Fn(&GffEntry) -> bool>(
        &mut self,
        selector: F,
        length: u32,
    ) {
        let mut entries = vec![];
        for (id, entry) in self.id_map.iter() {
            if selector(entry) {
                let upstream_end = entry.contig.start();
                let upstream_start = upstream_end.saturating_sub(length);
                let contig = Contig::new(
                    entry.contig.seqname(),
                    upstream_start,
                    upstream_end,
                    entry.contig.strand(),
                );

                let upstream_entry = GffEntry::new(
                    contig,
                    None,
                    Some(
                        if entry.feature_type.is_empty() {
                            UPSTREAM_TYPE_NAME.into()
                        }
                        else {
                            format!(
                                "{}_{}",
                                entry.feature_type, UPSTREAM_TYPE_NAME
                            )
                            .into()
                        },
                    ),
                    None,
                    None,
                    Some(
                        GffEntryAttributes::default()
                            .with_parent(Some(vec![id.clone()])),
                    ),
                );
                entries.push(upstream_entry);
            }
        }

        for entry in entries {
            self.insert(entry);
        }
    }

    pub fn add_downstream<F: Fn(&GffEntry) -> bool>(
        &mut self,
        selector: F,
        length: u32,
    ) {
        let mut entries = vec![];
        for (id, entry) in self.id_map.iter() {
            if selector(entry) {
                let downstream_start = entry.contig.end();
                let downstream_end = downstream_start + length;
                let contig = Contig::new(
                    entry.contig.seqname(),
                    downstream_start,
                    downstream_end,
                    entry.contig.strand(),
                );

                let downstream_entry = GffEntry::new(
                    contig,
                    None,
                    Some(
                        if entry.feature_type.is_empty() {
                            DOWNSTREAM_TYPE_NAME.into()
                        }
                        else {
                            format!(
                                "{}_{}",
                                entry.feature_type, DOWNSTREAM_TYPE_NAME
                            )
                            .into()
                        },
                    ),
                    None,
                    None,
                    Some(
                        GffEntryAttributes::default()
                            .with_parent(Some(vec![id.clone()])),
                    ),
                );
                entries.push(downstream_entry);
            }
        }

        for entry in entries {
            self.insert(entry);
        }
    }

    pub fn add_flanks<F: Fn(&GffEntry) -> bool + Clone>(
        &mut self,
        selector: F,
        length: u32,
    ) {
        self.add_upstream(selector.clone(), length);
        self.add_downstream(selector.clone(), length);
    }

    pub fn id_map(&self) -> &HashMap<BsxSmallStr, GffEntry> { &self.id_map }

    pub fn parent_map(&self) -> &MultiMap<BsxSmallStr, BsxSmallStr> {
        &self.parent_map
    }

    pub fn children_map(&self) -> &MultiMap<BsxSmallStr, BsxSmallStr> {
        &self.children_map
    }

    pub fn len(&self) -> usize { self.id_map.len() }
}

pub struct AnnotIterator<'a> {
    store: &'a AnnotStore,
    iter:  Box<dyn Iterator<Item = (BsxSmallStr, GffEntry)>>,
}

impl<'a> Iterator for AnnotIterator<'a> {
    type Item = (BsxSmallStr, GffEntry);

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.store.len(), Some(self.store.len()))
    }

    fn count(self) -> usize
    where
        Self: Sized, {
        self.store.len()
    }
}

impl FromIterator<GffEntry> for AnnotStore {
    fn from_iter<T: IntoIterator<Item = GffEntry>>(iter: T) -> Self {
        let mut new_self = Self::new();
        iter.into_iter().for_each(|entry| {
            new_self.insert(entry);
        });
        new_self
    }
}

impl FromIterator<(String, GffEntry)> for AnnotStore {
    fn from_iter<T: IntoIterator<Item = (String, GffEntry)>>(iter: T) -> Self {
        let mut new_self = Self::new();
        iter.into_iter()
            .for_each(|(_id, entry)| {
                new_self.insert(entry);
            });
        new_self
    }
}
