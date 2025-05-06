use hashbrown::HashMap;
use multimap::MultiMap;

use crate::data_structs::annotation::{GffEntry, GffEntryAttributes};
use crate::data_structs::coords::Contig;
use crate::data_structs::typedef::BsxSmallStr;

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
            iter:  self.id_map.iter(),
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
}

pub struct AnnotIterator<'a> {
    store: &'a AnnotStore,
    iter:  hashbrown::hash_map::Iter<'a, BsxSmallStr, GffEntry>,
}

impl<'a> Iterator for AnnotIterator<'a> {
    type Item = (&'a BsxSmallStr, &'a GffEntry);

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
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
