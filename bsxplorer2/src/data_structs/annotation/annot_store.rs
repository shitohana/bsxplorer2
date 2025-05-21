use std::error::Error;
use std::io::Read;
use std::ops::Range;

use arcstr::ArcStr;
use bio::data_structures::interval_tree::IntervalTree;
use hashbrown::HashMap;
use id_tree::{InsertBehavior, Node, NodeId, Tree, TreeBuilder};
use itertools::Itertools;
use log::warn;
use multimap::MultiMap;

use super::RawGffEntry;
use crate::data_structs::annotation::{GffEntry, GffEntryAttributes};
use crate::data_structs::coords::{Contig, GenomicPosition};
use crate::data_structs::typedef::BsxSmallStr;
use crate::io::bsx::BatchIndex;

const UPSTREAM_TYPE_NAME: &str = "upstream";
const DOWNSTREAM_TYPE_NAME: &str = "downstream";

#[derive(Clone, Debug)]
pub struct HcAnnotStore {
    ids:          Vec<ArcStr>,
    entries:      HashMap<ArcStr, GffEntry>,
    tree:         Tree<ArcStr>,
    tree_ids:     HashMap<ArcStr, NodeId>,
    tree_root_id: NodeId,
    interval_map: HashMap<BsxSmallStr, IntervalTree<GenomicPosition, ArcStr>>,
}

impl<A> From<A> for HcAnnotStore
where
    A: AsRef<[GffEntry]>,
{
    fn from(value: A) -> Self {
        let mut new = Self::with_capacity(value.as_ref().len());
        for entry in value.as_ref() {
            new.insert(entry.clone()).expect("Failed to insert node");
        }
        new
    }
}

impl FromIterator<GffEntry> for HcAnnotStore {
    fn from_iter<T: IntoIterator<Item = GffEntry>>(iter: T) -> Self {
        let mut store = Self::new();
        for entry in iter {
            store.insert(entry).expect("Failed to insert node");
        }
        store
    }
}

const TREE_ROOT_ID: ArcStr = arcstr::literal!("BSX_ROOT_NODE");

impl HcAnnotStore {
    fn new() -> Self {
        let mut tree = Tree::new();
        let root_id = tree
            .insert(Node::new(TREE_ROOT_ID.clone()), InsertBehavior::AsRoot)
            .unwrap();

        Self {
            ids: Vec::new(),
            entries: HashMap::new(),
            tree,
            tree_ids: HashMap::new(),
            tree_root_id: root_id,
            interval_map: HashMap::new(),
        }
    }

    fn from_gff<R: Read>(handle: R) -> Result<Self, Box<dyn Error>> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .flexible(true)
            .ascii()
            .from_reader(handle);

        let mut entries = Vec::new();
        for result in reader.deserialize() {
            let entry: RawGffEntry = result?;
            entries.push(GffEntry::try_from(entry)?);
        }

        Ok(Self::from_iter(entries))
    }

    pub fn from_bed<R: Read>(handle: R) -> Result<Self, Box<dyn Error>> {
        use bio::io::bed;

        let mut reader = bed::Reader::new(handle);
        let mut annot_store = Self::new();

        for record in reader.records() {
            let record = record?;

            let entry = GffEntry::from(record);
            annot_store.insert(entry)?;
        }

        Ok(annot_store)
    }

    fn with_capacity(capacity: usize) -> Self {
        let mut tree = TreeBuilder::new().with_node_capacity(capacity).build();
        let root_id = tree
            .insert(Node::new(TREE_ROOT_ID.clone()), InsertBehavior::AsRoot)
            .unwrap();

        Self {
            ids: Vec::with_capacity(capacity),
            entries: HashMap::with_capacity(capacity),
            tree,
            tree_ids: HashMap::with_capacity(capacity),
            tree_root_id: root_id,
            interval_map: HashMap::new(),
        }
    }

    /// # Note
    /// Name "BSX_ROOT_NODE" is reserved
    fn insert(
        &mut self,
        entry: GffEntry,
    ) -> Result<(), Box<dyn Error>> {
        let id = ArcStr::from(entry.id().as_str());
        self.ids.push(id.clone());

        self.interval_map
            .entry(entry.contig().seqname().to_owned())
            .and_modify(|tree| {
                tree.insert(Range::<_>::from(entry.contig()), id.clone());
            })
            .or_insert_with(|| {
                let mut tree = IntervalTree::new();
                tree.insert(Range::<_>::from(entry.contig()), id.clone());
                tree
            });

        if let Some(parent_ids) = entry.attributes().parent() {
            if !parent_ids.is_empty() {
                let parent = parent_ids.first().unwrap();
                if parent == "BSX_ROOT_NODE" {
                    return Err("Name \"BSX_ROOT_NODE\" is reserved".into());
                }

                let parent_node_id =
                    self.tree_ids.get(&ArcStr::from(parent.as_str())).unwrap();
                let new_node_id = self.tree.insert(
                    Node::new(id.clone()),
                    InsertBehavior::UnderNode(parent_node_id),
                )?;
                self.tree_ids.insert(id.clone(), new_node_id);

                if parent_ids.len() > 1 {
                    warn!(
                        "Multiple parent IDs are currently not supported. EntryID: {}",
                        id
                    );
                }
            }
            else {
                let new_node_id = self.tree.insert(
                    Node::new(id.clone()),
                    InsertBehavior::UnderNode(&self.tree_root_id),
                )?;
                self.tree_ids.insert(id.clone(), new_node_id);
            }
        }

        self.entries.insert(id.clone(), entry);
        Ok(())
    }

    fn genomic_query(
        &self,
        contig: &Contig,
    ) -> Option<Vec<ArcStr>> {
        self.interval_map.get(contig.seqname()).map(|tree| {
            tree.find(Range::<_>::from(contig))
                .map(|e| e.data().clone())
                .collect_vec()
        })
    }

    fn get_feature_types(&self) -> Vec<&str> {
        self.entries
            .values()
            .map(|e| e.feature_type().as_str())
            .unique()
            .collect()
    }

    fn get_entry(
        &self,
        id: &ArcStr,
    ) -> Option<&GffEntry> {
        self.entries.get(id)
    }

    fn get_nodeid(
        &self,
        id: &ArcStr,
    ) -> Option<&NodeId> {
        self.tree_ids.get(id)
    }

    fn get_node(
        &self,
        id: &ArcStr,
    ) -> Option<&Node<ArcStr>> {
        self.tree_ids
            .get(id)
            .and_then(|node_id| self.tree.get(node_id).ok())
    }

    fn get_children(
        &self,
        id: &ArcStr,
    ) -> Option<Vec<ArcStr>> {
        self.get_node(id).map(|node| {
            node.children()
                .into_iter()
                .map(|child_id| self.tree.get(child_id).unwrap().data().clone())
                .collect()
        })
    }

    fn get_parent(
        &self,
        id: &ArcStr,
    ) -> Option<ArcStr> {
        self.get_node(id)
            .map(|node| {
                node.parent()
                    .map(|parent| self.tree.get(parent).unwrap().data().clone())
            })
            .flatten()
    }

    fn sort_self(&mut self) -> Result<(), Box<dyn Error>> {
        let ids_list = self.tree_ids.values().cloned().collect_vec();
        for node_id in ids_list.iter() {
            self.tree.sort_children_by_key(node_id, |node| {
                self.entries
                    .get(node.data())
                    .map(|entry| entry.contig().start())
                    .unwrap_or(0)
            })?;
        }
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }
}

/// A store for genomic annotations.
///
/// Provides functionality to store, retrieve, and manipulate GFF entries.
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
    pub fn iter_sorted(
        &self,
        index: &BatchIndex,
    ) -> AnnotIterator {
        let iterator: Box<dyn Iterator<Item = _>> = Box::from(
            self.iter()
                .into_group_map_by(|(_id, entry)| {
                    index.get_chr_index(entry.contig.seqname())
                })
                .into_iter()
                .sorted_by_cached_key(|(key, _value)| *key)
                .flat_map(|(_key, value)| {
                    value
                        .into_iter()
                        .sorted_by_cached_key(|(_id, entry)| entry.contig.start())
                })
                .map(|(id, value)| (id.clone(), value.clone())),
        );

        AnnotIterator {
            store: self,
            iter:  iterator,
        }
    }

    /// Gets the parents of an entry.
    pub fn get_parents(
        &self,
        id: &BsxSmallStr,
    ) -> Option<&Vec<BsxSmallStr>> {
        self.parent_map.get_vec(id)
    }

    /// Gets the children of an entry.
    pub fn get_children(
        &self,
        id: &BsxSmallStr,
    ) -> Option<&Vec<BsxSmallStr>> {
        self.children_map.get_vec(id)
    }

    /// Adds upstream regions to selected entries.
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
                    BsxSmallStr::from(entry.contig.seqname().as_str()),
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
                    BsxSmallStr::from(entry.contig.seqname().as_str()),
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
    pub fn len(&self) -> usize {
        self.id_map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.id_map.is_empty()
    }
}

/// An iterator over the entries in an `AnnotStore`.
pub struct AnnotIterator<'a> {
    store: &'a AnnotStore,
    iter:  Box<dyn Iterator<Item = (BsxSmallStr, GffEntry)>>,
}

impl Iterator for AnnotIterator<'_> {
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


#[cfg(test)]
mod tests {
    use arcstr::ArcStr;

    use super::*;
    use crate::data_structs::enums::Strand;


    /// Helper function to create a test GffEntry.
    ///
    /// This function is not meant to be used directly but is used by the
    /// doctests.
    #[doc(hidden)]
    pub fn create_test_entry(
        id: &str,
        seqname: &str,
        start: u32,
        end: u32,
    ) -> GffEntry {
        let contig =
            Contig::new(BsxSmallStr::from(seqname), start, end, Strand::Forward);

        let attrs =
            GffEntryAttributes::default().with_id::<BsxSmallStr>(Some(id.into()));

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
    /// This function is not meant to be used directly but is used by the
    /// doctests.
    #[doc(hidden)]
    pub fn create_test_entry_with_parent(
        id: &str,
        seqname: &str,
        start: u32,
        end: u32,
        parent: Option<Vec<&str>>,
    ) -> GffEntry {
        let contig =
            Contig::new(BsxSmallStr::from(seqname), start, end, Strand::Forward);

        let mut attrs =
            GffEntryAttributes::default().with_id::<BsxSmallStr>(Some(id.into()));

        if let Some(parents) = parent {
            attrs = attrs.with_parent(Some(
                parents.into_iter().map(BsxSmallStr::from).collect(),
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

    /// Create a mock BatchIndex for testing
    fn create_mock_index() -> BatchIndex {
        let mut index = BatchIndex::new();

        // Insert contigs for chr1 and chr2 with batch indices
        index.insert(
            Contig::new(BsxSmallStr::from("chr1"), 0, 10000, Strand::Forward),
            0,
        );
        index.insert(
            Contig::new(BsxSmallStr::from("chr2"), 0, 20000, Strand::Forward),
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

    use id_tree::InsertBehavior::*;
    #[test]
    fn tree_test() {
        let mut tree: Tree<i32> = TreeBuilder::new().with_node_capacity(5).build();

        let root_id: NodeId = tree.insert(Node::new(0), AsRoot).unwrap();
        let child_id: NodeId = tree.insert(Node::new(1), UnderNode(&root_id)).unwrap();
        tree.insert(Node::new(2), UnderNode(&root_id)).unwrap();
        tree.insert(Node::new(3), UnderNode(&child_id)).unwrap();
        let child2_id = tree.insert(Node::new(4), UnderNode(&child_id)).unwrap();
        tree.insert(Node::new(3), UnderNode(&child2_id)).unwrap();

        assert_eq!(tree.children(&child_id).unwrap().count(), 2);
        assert_eq!(tree.children(&child2_id).unwrap().count(), 1);
        println!("{:?}", tree.ancestors(&child2_id).unwrap().collect_vec())
    }
}
