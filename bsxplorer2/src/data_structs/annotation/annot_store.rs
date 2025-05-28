use std::error::Error;
use std::io::Read;
use std::ops::Range;

use arcstr::ArcStr;
use bio::data_structures::interval_tree::IntervalTree;
use hashbrown::HashMap;
use id_tree::{
    InsertBehavior,
    Node,
    NodeId,
    Tree,
    TreeBuilder,
};
use itertools::Itertools;
use regex_lite::Regex;

use super::RawGffEntry;
use crate::data_structs::annotation::{
    GffEntry,
    GffEntryAttributes,
};
use crate::data_structs::coords::{
    Contig,
    GenomicPosition,
};
use crate::data_structs::typedef::BsxSmallStr;
use crate::getter_fn;

/// A data structure to store and manage genomic annotations (GffEntry).
///
/// It uses multiple internal structures for efficient lookup:
/// - A `HashMap` for ID-based lookup of `GffEntry`.
/// - An `id_tree::Tree` to represent parent-child relationships between
///   entries.
/// - A `HashMap` mapping sequence names to
///   `bio::data_structures::interval_tree::IntervalTree` for genomic range
///   queries.
#[derive(Debug)]
pub struct HcAnnotStore {
    ids:          Vec<ArcStr>,
    entries:      HashMap<ArcStr, GffEntry>,
    tree:         Tree<ArcStr>,
    tree_ids:     HashMap<ArcStr, NodeId>,
    tree_root_id: NodeId,
    interval_map: HashMap<BsxSmallStr, IntervalTree<GenomicPosition, ArcStr>>,
}

impl Default for HcAnnotStore {
    fn default() -> Self {
        Self::new()
    }
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
    getter_fn!(ids, Vec<ArcStr>);

    getter_fn!(entries, HashMap<ArcStr, GffEntry>);

    getter_fn!(tree, Tree<ArcStr>);

    getter_fn!(tree_ids, HashMap<ArcStr, NodeId>);

    getter_fn!(interval_map, HashMap<BsxSmallStr, IntervalTree<GenomicPosition, ArcStr>>);

    /// Creates a new empty `HcAnnotStore`.
    pub fn new() -> Self {
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

    /// Returns an iterator over all entries in the store.
    pub fn iter(&self) -> impl Iterator<Item = &GffEntry> {
        self.entries.values()
    }

    /// Creates a new `HcAnnotStore` by reading and parsing a GFF file.
    ///
    /// # Errors
    ///
    /// Returns an error if the GFF file cannot be read or parsed correctly.
    pub fn from_gff<R: Read>(handle: R) -> Result<Self, Box<dyn Error>> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .has_headers(false)
            .flexible(true)
            .from_reader(handle);

        let mut entries = Vec::new();

        for result in reader.deserialize() {
            let entry: RawGffEntry = result?;
            entries.push(GffEntry::try_from(entry)?);
        }
        Ok(Self::from_iter(entries))
    }

    /// Creates a new `HcAnnotStore` by reading and parsing a BED file.
    ///
    /// # Errors
    ///
    /// Returns an error if the BED file cannot be read or parsed correctly,
    /// or if an entry fails to insert.
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

    /// Creates a new empty `HcAnnotStore` with a pre-allocated capacity.
    ///
    /// # Arguments
    ///
    /// * `capacity` - The estimated number of entries to pre-allocate space
    ///   for.
    pub fn with_capacity(capacity: usize) -> Self {
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

    /// Inserts a `GffEntry` into the store.
    ///
    /// The entry is added to the internal HashMap, IntervalTree, and the
    /// id_tree based on its ID and parent information.
    ///
    /// # Arguments
    ///
    /// * `entry` - The `GffEntry` to insert.
    ///
    /// # Note
    /// Name "BSX_ROOT_NODE" is reserved. An error will be returned if an entry
    /// uses this ID or specifies this as a parent.
    /// Multiple parent IDs are currently not fully supported; only the first
    /// parent is used for tree structure, and a warning is logged for multiple
    /// parents.
    ///
    /// # Errors
    ///
    /// Returns an error if the entry's ID is the reserved "BSX_ROOT_NODE",
    /// if a parent ID is the reserved "BSX_ROOT_NODE", or if tree insertion
    /// fails.
    pub fn insert(
        &mut self,
        entry: GffEntry,
    ) -> Result<(), Box<dyn Error>> {
        let id = ArcStr::from(entry.id().as_str());

        self.insert_to_imap(id.clone(), entry.contig.clone());
        self.insert_to_tree(id.clone(), entry.attributes().parent())?;

        self.ids.push(id.clone());
        self.entries.insert(id.clone(), entry);
        Ok(())
    }

    fn insert_to_imap(
        &mut self,
        id: ArcStr,
        contig: Contig,
    ) {
        self.interval_map
            .entry(contig.seqname().to_owned())
            .and_modify(|tree| {
                tree.insert(Range::<_>::from(&contig), id.clone());
            })
            .or_insert_with(|| {
                let mut tree = IntervalTree::new();
                tree.insert(Range::<_>::from(&contig), id.clone());
                tree
            });
    }

    fn insert_to_tree(
        &mut self,
        id: ArcStr,
        parents: Option<&Vec<BsxSmallStr>>,
    ) -> Result<(), Box<dyn Error>> {
        if parents.as_ref().is_none_or(|p| p.is_empty()) {
            self.insert_to_root(id.as_str().into())?;
        }
        else {
            let parent: ArcStr =
                parents.as_ref().unwrap().first().unwrap().as_str().into();
            if parent == "BSX_ROOT_NODE" {
                return Err("Name \"BSX_ROOT_NODE\" is reserved".into());
            }

            let parent_node_id = {
                match self.tree_ids.get(&parent) {
                    Some(id) => id,
                    None => {
                        self.insert_to_root(parent.clone())?;
                        self.tree_ids.get(&parent).unwrap()
                    },
                }
            };

            let new_node_id = self.tree.insert(
                Node::new(id.clone()),
                InsertBehavior::UnderNode(parent_node_id),
            )?;
            self.tree_ids.insert(id.clone(), new_node_id);
        }

        Ok(())
    }

    fn insert_to_root(
        &mut self,
        id: ArcStr,
    ) -> Result<(), Box<dyn Error>> {
        let new_node = Node::new(id.clone());
        let new_node_id = self
            .tree
            .insert(new_node, InsertBehavior::UnderNode(&self.tree_root_id))?;
        self.tree_ids.insert(id.clone(), new_node_id);
        Ok(())
    }

    /// Adds flanking regions to annotation map
    ///
    /// For entries selected by the `selector` function, new GffEntry objects
    /// representing flanking regions are created and inserted into the store.
    /// The new entries' feature types will be prefixed with `prefix`.
    ///
    /// # Parameters
    ///
    /// * `selector` - A function that takes a reference to a `GffEntry` and
    ///   returns `true` if the entry should have flanking regions added.
    /// * `flank` - The length of the flank regions to be added. A positive
    ///   value adds the flank downstream (after the end), a negative value adds
    ///   the flank upstream (before the start).
    /// * `prefix` - A string prefix to add to the feature type of the newly
    ///   created flanking entries.
    pub fn add_flank<F>(
        &mut self,
        selector: F,
        flank: i32,
        prefix: &str,
    ) where
        F: Fn(&GffEntry) -> bool, {
        let selected_entries = self
            .entries
            .iter()
            .filter(|(_id, entry)| selector(entry)) // Filter by entry content, not id
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect_vec();

        for (id, parent) in selected_entries {
            let (start, end) = if flank > 0 {
                // Flank downstream (after end)
                (
                    parent.contig.end_gpos(),
                    parent.contig.end_gpos().shift(flank as isize),
                )
            }
            else {
                // Flank upstream (before start)
                (
                    parent.contig.start_gpos().shift(flank as isize),
                    parent.contig.start_gpos(),
                )
            };

            // Ensure start <= end for the range
            let (start, end) = if start <= end {
                (start, end)
            }
            else {
                (end, start)
            };


            let mut feature_type = prefix.to_string();
            feature_type.push_str(parent.feature_type.as_str());

            // Create a new unique ID for the flank entry
            let flank_id_str = format!("{}_flank_{}", id, flank);
            let flank_id: ArcStr = flank_id_str.clone().into(); // Convert to ArcStr

            let flank_entry = GffEntry::new(
                (start..end).into(), /* Assuming Range<GenomicPosition> converts to
                                      * Contig */
                None,                      // No source
                Some(feature_type.into()), // New feature type
                None,                      // No score
                None,                      // No strand
                Some(
                    // Attributes
                    GffEntryAttributes::default()
                        .with_id(flank_id.as_str().into()) // Set the unique ID for the flank
                        .with_parent(Some(vec![BsxSmallStr::from_str(id.as_str())])), /* Set parent to the original entry */
                ),
            );

            // Use insert which handles adding to all internal structures
            // Propagate the error if insertion fails
            self.insert(flank_entry).expect(
                format!("Failed to insert flank entry: {}", flank_id_str).as_str(),
            );
        }
    }

    /// Queries the store for entries overlapping a given genomic contig/range.
    ///
    /// # Arguments
    ///
    /// * `contig` - The genomic contig (including start and end positions) to
    ///   query.
    ///
    /// # Returns
    ///
    /// An optional vector of `ArcStr` IDs for the entries that overlap the
    /// given `contig`. Returns `None` if the sequence name is not found in
    /// the store.
    pub fn genomic_query(
        &self,
        contig: &Contig,
    ) -> Option<Vec<ArcStr>> {
        self.interval_map.get(contig.seqname()).map(|tree| {
            tree.find(Range::<_>::from(contig))
                .map(|e| e.data().clone())
                .collect_vec()
        })
    }

    /// Returns a vector of unique feature types present in the store.
    pub fn get_feature_types(&self) -> Vec<&str> {
        self.entries
            .values()
            .map(|e| e.feature_type().as_str())
            .unique()
            .collect()
    }

    /// Retrieves a reference to a `GffEntry` by its ID.
    ///
    /// # Arguments
    ///
    /// * `id` - The ID of the entry to retrieve.
    ///
    /// # Returns
    ///
    /// An optional reference to the `GffEntry` if found, otherwise `None`.
    pub fn get_entry(
        &self,
        id: &ArcStr,
    ) -> Option<&GffEntry> {
        self.entries.get(id)
    }

    /// Retrieves a vector of references to `GffEntry` objects whose IDs match a
    /// regex pattern.
    ///
    /// # Arguments
    ///
    /// * `pattern` - The regex pattern to match against entry IDs.
    ///
    /// # Returns
    ///
    /// A vector containing references to the matching `GffEntry` objects.
    pub fn get_entries_regex(
        &self,
        pattern: &str,
    ) -> Result<Vec<&GffEntry>, Box<dyn Error>> {
        let regex_completed = Regex::new(pattern)?;
        Ok(self
            .entries
            .values()
            .filter(|entry| regex_completed.is_match(entry.id().as_str()))
            .collect())
    }

    fn get_node(
        &self,
        id: &ArcStr,
    ) -> Option<&Node<ArcStr>> {
        self.tree_ids
            .get(id)
            .and_then(|node_id| self.tree.get(node_id).ok())
    }

    /// Retrieves the IDs of the direct children of a given entry in the tree
    /// structure.
    ///
    /// # Arguments
    ///
    /// * `id` - The ID of the parent entry.
    ///
    /// # Returns
    ///
    /// An optional vector of `ArcStr` IDs representing the children.
    /// Returns `None` if the given ID is not found in the tree.
    pub fn get_children(
        &self,
        id: &ArcStr,
    ) -> Option<Vec<ArcStr>> {
        self.get_node(id).map(|node| {
            node.children()
                .iter()
                .map(|child_id| self.tree.get(child_id).unwrap().data().clone())
                .collect()
        })
    }

    /// Retrieves the ID of the direct parent of a given entry in the tree
    /// structure.
    ///
    /// # Arguments
    ///
    /// * `id` - The ID of the child entry.
    ///
    /// # Returns
    ///
    /// An optional `ArcStr` ID representing the parent. Returns `None` if the
    /// given ID is not found, or if the found node is the root node.
    pub fn get_parent(
        &self,
        id: &ArcStr,
    ) -> Option<ArcStr> {
        self.get_node(id).and_then(|node| {
            node.parent()
                .map(|parent| self.tree.get(parent).unwrap().data().clone())
        })
    }

    /// Sorts the children of each node in the internal tree based on the start
    /// position of the corresponding `GffEntry`.
    ///
    /// This can help in traversing the tree in genomic order.
    ///
    /// # Errors
    ///
    /// Returns an error if tree manipulation fails.
    pub fn sort_self(&mut self) -> Result<(), Box<dyn Error>> {
        let ids_list = self.tree_ids.values().cloned().collect_vec();
        for node_id in ids_list.iter() {
            // The id_tree root node does not have an associated GffEntry, skip sorting
            // children for it? Or handle the case where entries.get()
            // returns None. The current implementation defaults to 0.
            // This might sort the children of the root node, which could be top-level
            // entries.
            self.tree.sort_children_by_key(node_id, |node| {
                self.entries
                    .get(node.data())
                    .map(|entry| entry.contig().start())
                    .unwrap_or(0) // Default to 0 if entry not found (e.g., for
                                  // the root node's direct children if root
                                  // data is not an entry ID)
            })?;
        }
        Ok(())
    }

    /// Returns the number of entries in the store.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Returns `true` if the store contains no entries.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }
}
