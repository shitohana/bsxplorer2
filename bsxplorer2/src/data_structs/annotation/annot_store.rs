use hashbrown::HashMap;
use itertools::Itertools;
use multimap::MultiMap;

use crate::data_structs::annotation::{GffEntry, GffEntryAttributes};
use crate::data_structs::coords::Contig;
use crate::data_structs::enums::Strand;
use crate::data_structs::typedef::{BsxSmallStr, SeqNameStr, SeqPosNum};
use crate::io::bsx::BatchIndex;

const UPSTREAM_TYPE_NAME: &str = "upstream";
const DOWNSTREAM_TYPE_NAME: &str = "downstream";

/// A store for genomic annotations.
///
/// Provides functionality to store, retrieve, and manipulate GFF entries.
///
/// # Examples
///
/// ```
/// use bsxplorer2::data_structs::annotation::annot_store::AnnotStore;
/// use bsxplorer2::data_structs::annotation::GffEntry;
/// # use bsxplorer2::data_structs::annotation::annot_store::create_test_entry;
///
/// // Create a new annotation store
/// let mut store = AnnotStore::new();
///
/// // Add entries to the store
/// let entry = create_test_entry("gene1", "chr1", 1000, 2000);
/// store.insert(entry);
///
/// // Check length
/// assert_eq!(store.len(), 1);
/// ```
pub struct AnnotStore {
    id_map:       HashMap<BsxSmallStr, GffEntry>,
    parent_map:   MultiMap<BsxSmallStr, BsxSmallStr>,
    children_map: MultiMap<BsxSmallStr, BsxSmallStr>,
}

impl Default for AnnotStore {
    fn default() -> Self {
        Self::new()
    }
}

impl AnnotStore {
    /// Creates a new empty `AnnotStore`.
    ///
    /// # Examples
    ///
    /// ```
    /// use bsxplorer2::data_structs::annotation::annot_store::AnnotStore;
    ///
    /// let store = AnnotStore::new();
    /// assert_eq!(store.len(), 0);
    /// ```
    pub fn new() -> Self {
        Self {
            id_map:       Default::default(),
            parent_map:   Default::default(),
            children_map: Default::default(),
        }
    }

    /// Inserts a GFF entry into the store.
    ///
    /// Returns `Some(())` if the entry was inserted successfully, or `None` if
    /// an entry with the same ID already exists.
    ///
    /// # Examples
    ///
    /// ```
    /// use bsxplorer2::data_structs::annotation::annot_store::AnnotStore;
    /// use bsxplorer2::data_structs::annotation::GffEntry;
    /// # use bsxplorer2::data_structs::annotation::annot_store::create_test_entry;
    ///
    /// let mut store = AnnotStore::new();
    ///
    /// // Insert an entry
    /// let entry1 = create_test_entry("gene1", "chr1", 1000, 2000);
    /// assert!(store.insert(entry1).is_some());
    ///
    /// // Try to insert an entry with the same ID
    /// let entry2 = create_test_entry("gene1", "chr1", 3000, 4000);
    /// assert!(store.insert(entry2).is_none());
    /// ```
    pub fn insert(
        &mut self,
        entry: GffEntry,
    ) -> Option<()> {
        if self.id_map.contains_key(&entry.id) {
            None
        }
        else {
            self.id_map.insert(entry.id.clone(), entry.clone());
            for parent in entry.attributes.parent.as_ref().unwrap_or(&vec![]) {
                self.parent_map.insert(entry.id.clone(), parent.to_owned());
                self.children_map
                    .insert(parent.to_owned(), entry.id.clone());
            }
            Some(())
        }
    }

    /// Removes a GFF entry from the store by ID.
    ///
    /// Returns the removed entry if found, or `None` if no entry with the given
    /// ID exists.
    ///
    /// # Examples
    ///
    /// ```
    /// use bsxplorer2::data_structs::annotation::annot_store::AnnotStore;
    /// # use bsxplorer2::data_structs::annotation::annot_store::create_test_entry;
    ///
    /// let mut store = AnnotStore::new();
    /// let entry = create_test_entry("gene1", "chr1", 1000, 2000);
    /// store.insert(entry);
    ///
    /// // Remove the entry
    /// let id: bsxplorer2::data_structs::typedef::BsxSmallStr = "gene1".into();
    /// let removed = store.remove(&id);
    /// assert!(removed.is_some());
    /// assert_eq!(store.len(), 0);
    /// ```
    pub fn remove(
        &mut self,
        id: &BsxSmallStr,
    ) -> Option<GffEntry> {
        let removed = self.id_map.remove(id);
        let _removed_parents = self.parent_map.remove(id);
        let _removed_children = self.children_map.remove(id);

        removed
    }

    /// Returns an iterator over all entries in the store.
    pub fn iter(&self) -> AnnotIterator {
        AnnotIterator {
            store: self,
            iter:  Box::from(self.id_map.clone().into_iter()),
        }
    }

    /// Returns an iterator over all entries, sorted by chromosome and position.
    pub fn iter_sorted<S: SeqNameStr, P: SeqPosNum>(
        &self,
        index: &BatchIndex<S, P>,
    ) -> AnnotIterator {
        let iterator: Box<dyn Iterator<Item = _>> = Box::from(
            self.iter()
                .into_group_map_by(|(_id, entry)| {
                    index.get_chr_index(&S::from(entry.contig.seqname().as_str()))
                })
                .into_iter()
                .sorted_by_cached_key(|(key, _value)| *key)
                .map(|(_key, value)| {
                    value
                        .into_iter()
                        .sorted_by_cached_key(|(_id, entry)| entry.contig.start())
                })
                .flatten()
                .map(|(id, value)| (id.clone(), value.clone())),
        );

        AnnotIterator {
            store: self,
            iter:  iterator,
        }
    }

    /// Gets the parents of an entry.
    ///
    /// # Examples
    ///
    /// ```
    /// use bsxplorer2::data_structs::annotation::annot_store::AnnotStore;
    /// # use bsxplorer2::data_structs::annotation::annot_store::create_test_entry_with_parent;
    ///
    /// let mut store = AnnotStore::new();
    ///
    /// // Create parent entry
    /// let parent = create_test_entry_with_parent("gene1", "chr1", 1000, 5000, None);
    /// store.insert(parent);
    ///
    /// // Create child entry with parent reference
    /// let child = create_test_entry_with_parent("exon1", "chr1", 1000, 2000, Some(vec!["gene1"]));
    /// store.insert(child);
    ///
    /// // Get parents of the child
    /// let id: bsxplorer2::data_structs::typedef::BsxSmallStr = "exon1".into();
    /// let parents = store.get_parents(&id);
    /// assert!(parents.is_some());
    /// assert_eq!(parents.unwrap()[0], "gene1");
    /// ```
    pub fn get_parents(
        &self,
        id: &BsxSmallStr,
    ) -> Option<&Vec<BsxSmallStr>> {
        self.parent_map.get_vec(id)
    }

    /// Gets the children of an entry.
    ///
    /// # Examples
    ///
    /// ```
    /// use bsxplorer2::data_structs::annotation::annot_store::AnnotStore;
    /// # use bsxplorer2::data_structs::annotation::annot_store::create_test_entry_with_parent;
    ///
    /// let mut store = AnnotStore::new();
    ///
    /// // Create parent entry
    /// let parent = create_test_entry_with_parent("gene1", "chr1", 1000, 5000, None);
    /// store.insert(parent);
    ///
    /// // Create child entry with parent reference
    /// let child = create_test_entry_with_parent("exon1", "chr1", 1000, 2000, Some(vec!["gene1"]));
    /// store.insert(child);
    ///
    /// // Get children of the parent
    /// let id: bsxplorer2::data_structs::typedef::BsxSmallStr = "gene1".into();
    /// let children = store.get_children(&id);
    /// assert!(children.is_some());
    /// assert_eq!(children.unwrap()[0], "exon1");
    /// ```
    pub fn get_children(
        &self,
        id: &BsxSmallStr,
    ) -> Option<&Vec<BsxSmallStr>> {
        self.children_map.get_vec(id)
    }

    /// Adds upstream regions to selected entries.
    ///
    /// # Examples
    ///
    /// ```
    /// use bsxplorer2::data_structs::annotation::annot_store::AnnotStore;
    /// # use bsxplorer2::data_structs::annotation::annot_store::create_test_entry;
    ///
    /// let mut store = AnnotStore::new();
    /// let entry = create_test_entry("gene1", "chr1", 1000, 2000);
    /// store.insert(entry);
    ///
    /// // Add upstream regions for all entries
    /// store.add_upstream(|_| true, 500);
    ///
    /// // Should now have 2 entries: the original and the upstream region
    /// assert_eq!(store.len(), 2);
    ///
    /// // Verify that at least one entry has "upstream" in its feature type
    /// let has_upstream = store.iter().any(|(_, entry)|
    ///     entry.feature_type.contains("upstream"));
    /// assert!(has_upstream);
    /// ```
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
                    entry.contig.seqname().to_owned(),
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
                            format!("{}_{}", entry.feature_type, UPSTREAM_TYPE_NAME)
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

    /// Adds downstream regions to selected entries.
    ///
    /// # Examples
    ///
    /// ```
    /// use bsxplorer2::data_structs::annotation::annot_store::AnnotStore;
    /// # use bsxplorer2::data_structs::annotation::annot_store::create_test_entry;
    ///
    /// let mut store = AnnotStore::new();
    /// let entry = create_test_entry("gene1", "chr1", 1000, 2000);
    /// store.insert(entry);
    ///
    /// // Add downstream regions for all entries
    /// store.add_downstream(|_| true, 500);
    ///
    /// // Should now have 2 entries: the original and the downstream region
    /// assert_eq!(store.len(), 2);
    ///
    /// // Verify that at least one entry has "downstream" in its feature type
    /// let has_downstream = store.iter().any(|(_, entry)|
    ///     entry.feature_type.contains("downstream"));
    /// assert!(has_downstream);
    /// ```
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
                    entry.contig.seqname().to_owned(),
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
                            format!("{}_{}", entry.feature_type, DOWNSTREAM_TYPE_NAME)
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

    /// Adds both upstream and downstream regions to selected entries.
    ///
    /// # Examples
    ///
    /// ```
    /// use bsxplorer2::data_structs::annotation::annot_store::AnnotStore;
    /// # use bsxplorer2::data_structs::annotation::annot_store::create_test_entry;
    /// # use arcstr::ArcStr;
    ///
    /// let mut store = AnnotStore::new();
    /// let entry = create_test_entry("gene1", "chr1", 1000, 2000);
    /// store.insert(entry);
    ///
    /// // Add both upstream and downstream regions for all entries
    /// store.add_flanks(|entry| entry.feature_type == ArcStr::from("gene"), 500);
    ///
    /// // Should now have 3 entries: the original, upstream, and downstream regions
    /// assert_eq!(store.len(), 3);
    /// ```
    pub fn add_flanks<F: Fn(&GffEntry) -> bool + Clone>(
        &mut self,
        selector: F,
        length: u32,
    ) {
        self.add_upstream(selector.clone(), length);
        self.add_downstream(selector.clone(), length);
    }

    /// Returns a reference to the internal ID map.
    pub fn id_map(&self) -> &HashMap<BsxSmallStr, GffEntry> {
        &self.id_map
    }

    /// Returns a reference to the internal parent map.
    pub fn parent_map(&self) -> &MultiMap<BsxSmallStr, BsxSmallStr> {
        &self.parent_map
    }

    /// Returns a reference to the internal children map.
    pub fn children_map(&self) -> &MultiMap<BsxSmallStr, BsxSmallStr> {
        &self.children_map
    }

    /// Returns the number of entries in the store.
    ///
    /// # Examples
    ///
    /// ```
    /// use bsxplorer2::data_structs::annotation::annot_store::AnnotStore;
    /// # use bsxplorer2::data_structs::annotation::annot_store::create_test_entry;
    ///
    /// let mut store = AnnotStore::new();
    /// assert_eq!(store.len(), 0);
    ///
    /// store.insert(create_test_entry("gene1", "chr1", 1000, 2000));
    /// assert_eq!(store.len(), 1);
    /// ```
    pub fn len(&self) -> usize {
        self.id_map.len()
    }
}

/// An iterator over the entries in an `AnnotStore`.
pub struct AnnotIterator<'a> {
    store: &'a AnnotStore,
    iter:  Box<dyn Iterator<Item = (BsxSmallStr, GffEntry)>>,
}

impl<'a> Iterator for AnnotIterator<'a> {
    type Item = (BsxSmallStr, GffEntry);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }

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
        iter.into_iter().for_each(|(_id, entry)| {
            new_self.insert(entry);
        });
        new_self
    }
}

/// Helper function to create a test GffEntry.
///
/// This function is not meant to be used directly but is used by the doctests.
#[doc(hidden)]
pub fn create_test_entry(
    id: &str,
    seqname: &str,
    start: u32,
    end: u32,
) -> GffEntry {
    let contig =
        Contig::new(arcstr::ArcStr::from(seqname), start, end, Strand::Forward);

    let attrs = GffEntryAttributes::default().with_id::<BsxSmallStr>(Some(id.into()));

    GffEntry::new(
        contig,
        Some(arcstr::ArcStr::from("test")),
        Some(arcstr::ArcStr::from("gene")),
        None,
        None,
        Some(attrs),
    )
}

/// Helper function to create a test GffEntry with parent information.
///
/// This function is not meant to be used directly but is used by the doctests.
#[doc(hidden)]
pub fn create_test_entry_with_parent(
    id: &str,
    seqname: &str,
    start: u32,
    end: u32,
    parent: Option<Vec<&str>>,
) -> GffEntry {
    let contig =
        Contig::new(arcstr::ArcStr::from(seqname), start, end, Strand::Forward);

    let mut attrs =
        GffEntryAttributes::default().with_id::<BsxSmallStr>(Some(id.into()));

    if let Some(parents) = parent {
        attrs = attrs.with_parent(Some(
            parents.into_iter().map(|p| BsxSmallStr::from(p)).collect(),
        ));
    }

    GffEntry::new(
        contig,
        Some(arcstr::ArcStr::from("test")),
        Some(arcstr::ArcStr::from("gene")),
        None,
        None,
        Some(attrs),
    )
}

#[cfg(test)]
mod tests {
    use arcstr::ArcStr;

    use super::*;

    /// Create a mock BatchIndex for testing
    fn create_mock_index() -> BatchIndex<ArcStr, u32> {
        let mut index = BatchIndex::new();

        // Insert contigs for chr1 and chr2 with batch indices
        index.insert(
            Contig::new(ArcStr::from("chr1"), 0, 10000, Strand::Forward),
            0,
        );
        index.insert(
            Contig::new(ArcStr::from("chr2"), 0, 20000, Strand::Forward),
            1,
        );

        index
    }

    #[test]
    fn test_insert_and_remove() {
        let mut store = AnnotStore::new();

        // Insert an entry
        let entry = create_test_entry("gene1", "chr1", 1000, 2000);
        assert!(store.insert(entry.clone()).is_some());
        assert_eq!(store.len(), 1);

        // Try to insert the same entry again
        assert!(store.insert(entry).is_none());
        assert_eq!(store.len(), 1);

        // Remove the entry
        let id: BsxSmallStr = "gene1".into();
        let removed = store.remove(&id);
        assert!(removed.is_some());
        assert_eq!(store.len(), 0);
    }

    #[test]
    fn test_parent_child_relationships() {
        let mut store = AnnotStore::new();

        // Create parent entry
        let parent = create_test_entry_with_parent("gene1", "chr1", 1000, 5000, None);
        store.insert(parent);

        // Create child entries with parent references
        let child1 = create_test_entry_with_parent(
            "exon1",
            "chr1",
            1000,
            2000,
            Some(vec!["gene1"]),
        );
        let child2 = create_test_entry_with_parent(
            "exon2",
            "chr1",
            3000,
            4000,
            Some(vec!["gene1"]),
        );
        store.insert(child1);
        store.insert(child2);

        // Check parent relationships
        let parent_id: BsxSmallStr = "gene1".into();
        let children = store.get_children(&parent_id).unwrap();
        assert_eq!(children.len(), 2);
        assert!(children.contains(&"exon1".into()));
        assert!(children.contains(&"exon2".into()));

        // Check child relationships
        let child_id: BsxSmallStr = "exon1".into();
        let parents = store.get_parents(&child_id).unwrap();
        assert_eq!(parents.len(), 1);
        assert_eq!(parents[0], "gene1");
    }

    #[test]
    fn test_add_upstream_downstream() {
        let mut store = AnnotStore::new();
        let entry = create_test_entry("gene1", "chr1", 1000, 2000);
        store.insert(entry);

        // Add upstream regions
        store.add_upstream(|entry| entry.feature_type == ArcStr::from("gene"), 500);
        assert_eq!(store.len(), 2);

        // Add downstream regions
        store.add_downstream(|entry| entry.feature_type == ArcStr::from("gene"), 500);
        assert_eq!(store.len(), 3);

        // Check types
        let has_upstream = store
            .iter()
            .any(|(_, entry)| entry.feature_type.contains("upstream"));
        let has_downstream = store
            .iter()
            .any(|(_, entry)| entry.feature_type.contains("downstream"));
        assert!(has_upstream);
        assert!(has_downstream);
    }

    #[test]
    fn test_iter_sorted() {
        let mut store = AnnotStore::new();

        // Add entries in different chromosomes
        store.insert(create_test_entry("gene1", "chr1", 3000, 4000));
        store.insert(create_test_entry("gene2", "chr1", 1000, 2000));
        store.insert(create_test_entry("gene3", "chr2", 1000, 2000));

        // Create mock index
        let index = create_mock_index();

        // Get sorted iterator
        let sorted_entries: Vec<_> = store.iter_sorted(&index).collect();

        // Check sorting - first by chromosome, then by start position
        assert_eq!(sorted_entries[0].0, "gene2"); // chr1, pos 1000
        assert_eq!(sorted_entries[1].0, "gene1"); // chr1, pos 3000
        assert_eq!(sorted_entries[2].0, "gene3"); // chr2, pos 1000
    }
}
