#![allow(unused)]

use std::collections::VecDeque;
use std::error::Error;
use std::io::Read;
use std::ops::Range;
use std::sync::Arc;

use anyhow::{
    anyhow,
    bail,
};
use arcstr::ArcStr;
use hashbrown::HashMap;
use id_tree::{
    InsertBehavior,
    Node,
    NodeId,
    NodeIdError,
    RemoveBehavior,
    Tree,
    TreeBuilder,
};
use itertools::{
    Either,
    Itertools,
};
use regex_lite::Regex;
use slotmap::{
    new_key_type,
    KeyData,
    SlotMap,
};

use super::RawGffEntry;
use crate::data_structs::annotation::{
    GffEntry,
    GffEntryAttributes,
};
use crate::data_structs::coords::{
    Contig,
    ContigIntervalMap,
    GenomicPosition,
};
use crate::data_structs::typedef::{
    BsxSmallStr,
    PosType,
};
use crate::getter_fn;

// TODO: Somehow change the parent-child id find, so we dont have to
// generate ids as strings

new_key_type! {
    pub struct EntryId;
}

struct EntryTree {
    tree:          Tree<EntryId>,
    tree_node_ids: HashMap<EntryId, Arc<NodeId>>,
    tree_root_id:  Arc<NodeId>,
}

impl FromIterator<(EntryId, Option<EntryId>)> for EntryTree {
    fn from_iter<T: IntoIterator<Item = (EntryId, Option<EntryId>)>>(iter: T) -> Self {
        let mut new_self = Self::new();
        new_self.append(iter).unwrap();
        new_self
    }
}

impl Default for EntryTree {
    fn default() -> Self {
        Self::new()
    }
}

impl EntryTree {
    fn new() -> Self {
        let mut tree = Tree::new();
        let tree_root = EntryId::from(KeyData::from_ffi(u64::MAX));
        let tree_root_node = Node::new(tree_root);
        let tree_root_id =
            Arc::new(tree.insert(tree_root_node, InsertBehavior::AsRoot).unwrap());
        let tree_node_ids = HashMap::from_iter([(tree_root, tree_root_id.clone())]);

        Self {
            tree,
            tree_root_id,
            tree_node_ids,
        }
    }

    pub fn append<I: IntoIterator<Item = (EntryId, Option<EntryId>)>>(
        &mut self,
        iter: I,
    ) -> anyhow::Result<()> {
        let mut queue = VecDeque::from_iter(iter);
        let mut last_len = queue.len();
        let mut count = 0;

        while let Some((child, parent)) = queue.pop_front() {
            if let Some(parent_id) = parent {
                if !self.tree_node_ids.contains_key(&parent_id) {
                    queue.push_back((child, Some(parent_id)));
                }
                else {
                    self.insert_under(child, parent_id).unwrap()
                }
            }
            else {
                self.insert_to_root(child).unwrap()
            }

            count += 1;

            if count == last_len && queue.len() >= last_len {
                bail!("Some childs have unexistent parents")
            }
            else {
                count = 0;
                last_len = queue.len();
            }
        }

        Ok(())
    }

    fn insert_to_root(
        &mut self,
        id: EntryId,
    ) -> Result<(), NodeIdError> {
        if !self.tree_node_ids.contains_key(&id) {
            let node = Node::new(id);
            let node_id =
                Arc::new(self.tree.insert(
                    node,
                    InsertBehavior::UnderNode(self.tree_root_id.as_ref()),
                )?);
            self.tree_node_ids.insert(id, node_id);
        }
        Ok(())
    }

    fn insert_under(
        &mut self,
        child_id: EntryId,
        parent_id: EntryId,
    ) -> anyhow::Result<()> {
        if !self.tree_node_ids.contains_key(&parent_id) {
            bail!("Parent id {:?} does not exist in the tree", parent_id)
        }
        match self.get_parent(child_id) {
            Some(existing_parent_id) if existing_parent_id != parent_id => {
                bail!(
                    "This child id {:?} has already got parent {:?}",
                    child_id,
                    existing_parent_id
                )
            },
            Some(existing_parent_id) if existing_parent_id == parent_id => {
                return Ok(())
            },
            _ => {},
        }

        let node = Node::new(child_id);
        let parent = unsafe { self.tree_node_ids.get(&parent_id).unwrap_unchecked() };
        let node_id = self.tree.insert(node, InsertBehavior::UnderNode(&parent))?;

        self.tree_node_ids.insert(child_id, node_id.into());
        Ok(())
    }

    fn remove(
        &mut self,
        id: EntryId,
    ) -> anyhow::Result<()> {
        self.tree_node_ids
            .remove(&id)
            .ok_or(anyhow!("Such id did not exist in a tree"))
            .and_then(|node_id| {
                Arc::try_unwrap(node_id)
                    .map_err(|_| anyhow!("Other references to NodeId exist"))
            })
            .and_then(|node_id| {
                self.tree
                    .remove_node(node_id, RemoveBehavior::DropChildren)
                    .map_err(|e| anyhow!(e))
            })
            .map(|_| ())
    }

    fn get_node(
        &self,
        id: EntryId,
    ) -> Option<&Node<EntryId>> {
        self.tree_node_ids
            .get(&id)
            .map(|node_id| self.tree.get(&node_id).expect("Should not fail"))
    }

    pub fn get_parent(
        &self,
        child_id: EntryId,
    ) -> Option<EntryId> {
        self.get_node(child_id)
            .and_then(Node::parent)
            .map(|parent_node_id| {
                self.tree
                    .get(parent_node_id)
                    .expect("Should not fail")
                    .data()
                    .to_owned()
            })
    }

    pub fn get_children(
        &self,
        parent_id: EntryId,
    ) -> Option<Vec<EntryId>> {
        self.get_node(parent_id).map(|node| {
            node.children()
                .iter()
                .map(|child_id| {
                    self.tree
                        .get(child_id)
                        .expect("Should not fail")
                        .data()
                        .to_owned()
                })
                .collect_vec()
        })
    }
}

pub struct HcAnnotStore {
    entries:      SlotMap<EntryId, GffEntry>,
    tree:         Either<EntryTree, Option<EntryTree>>,
    interval_map:
        Either<ContigIntervalMap<EntryId>, Option<ContigIntervalMap<EntryId>>>,
}

impl Default for HcAnnotStore {
    fn default() -> Self {
        Self {
            entries:      Default::default(),
            tree:         Either::Right(None),
            interval_map: Either::Right(None),
        }
    }
}

impl FromIterator<GffEntry> for HcAnnotStore {
    fn from_iter<T: IntoIterator<Item = GffEntry>>(iter: T) -> Self {
        let mut new = Self::new();
        for e in iter {
            new.insert(e);
        }
        new
    }
}

impl HcAnnotStore {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn keys(&self) -> slotmap::basic::Keys<'_, EntryId, GffEntry> {
        self.entries.keys()
    }

    pub fn values(&self) -> slotmap::basic::Values<'_, EntryId, GffEntry> {
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
            annot_store.insert(entry);
        }

        Ok(annot_store)
    }

    fn invalidate_maps(&mut self) {
        if let Either::Left(val) =
            std::mem::replace(&mut self.tree, Either::Right(None))
        {
            self.tree = Either::Right(Some(val))
        }
        if let Either::Left(val) =
            std::mem::replace(&mut self.interval_map, Either::Right(None))
        {
            self.interval_map = Either::Right(Some(val))
        }
    }

    /// Inserts a `GffEntry` into the store.
    ///
    /// The entry is added to the internal HashMap, IntervalTree, and the
    /// id_tree based on its ID and parent information.
    ///
    /// # Arguments
    /// * `entry` - The `GffEntry` to insert.
    pub fn insert(
        &mut self,
        entry: GffEntry,
    ) -> EntryId {
        self.invalidate_maps();
        self.entries.insert(entry)
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
            let flank_id_str = format!("{}_flank_{}", parent.id(), flank);
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
                        .with_id(flank_id) // Set the unique ID for the flank
                        .with_parent(vec![parent.id().clone()]), /* Set parent to
                                                                  * the original
                                                                  * entry */
                ),
            );

            // Use insert which handles adding to all internal structures
            // Propagate the error if insertion fails
            self.insert(flank_entry);
        }
    }

    pub fn init_tree(&mut self) -> anyhow::Result<()> {
        if self.tree.is_left() {
            return Ok(());
        }

        let mut tree = std::mem::replace(&mut self.tree, Either::Right(None))
            .right()
            .unwrap()
            .unwrap_or_default();

        let gffid2entryid: HashMap<&smallstr::SmallString<[u8; 20]>, EntryId> =
            HashMap::from_iter(
                self.entries
                    .iter()
                    .flat_map(|(k, v)| v.attributes().id().zip(Some(k))),
            );
        let iter = self
            .entries
            .iter()
            .map(|(id, entry)| {
                match entry.attributes().parent() {
                    Some(parents) if parents.is_empty() => Ok((id, None)),
                    Some(parents) => {
                        let first_parent = parents.first().unwrap();
                        gffid2entryid
                            .get(first_parent)
                            .ok_or(anyhow!("No such parent id {}", first_parent))
                            .map(|(parent)| (id, Some(parent.clone())))
                    },
                    None => Ok((id, None)),
                }
            })
            .collect::<anyhow::Result<Vec<_>>>()?;

        tree.append(iter)?;
        self.tree = Either::Left(tree);

        Ok(())
    }

    pub fn init_imap(&mut self) {
        if self.tree.is_left() {
            return;
        }

        let mut imap = std::mem::replace(&mut self.interval_map, Either::Right(None))
            .right()
            .unwrap()
            .unwrap_or_default();

        let iter = self
            .entries
            .iter()
            .map(|(id, entry)| (entry.contig().clone(), id.clone()));
        let new_imap = ContigIntervalMap::from_iter(iter);
        imap.union(&new_imap);

        self.interval_map = Either::Left(imap);
    }

    pub fn iter(&self) -> slotmap::basic::Iter<EntryId, GffEntry> {
        self.entries.iter()
    }

    pub fn get(
        &self,
        id: EntryId,
    ) -> Option<&GffEntry> {
        self.entries.get(id)
    }

    pub fn get_feature_types(&self) -> std::collections::HashMap<&str, Vec<EntryId>> {
        self.iter()
            .map(|(id, e)| (e.feature_type().as_str(), id))
            .into_group_map()
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
        let regex_compliled = Regex::new(pattern)?;
        Ok(self
            .entries
            .values()
            .filter(|entry| regex_compliled.is_match(entry.id().as_str()))
            .collect())
    }

    /// Retrieves the IDs of the direct children of a given entry in the tree
    /// structure.
    ///
    /// # Arguments
    ///
    /// * `id` - The ID of the parent entry.
    pub fn get_children(
        &self,
        id: EntryId,
    ) -> anyhow::Result<Option<Vec<EntryId>>> {
        if let Some(tree) = self.tree.as_ref().left() {
            Ok(tree.get_children(id))
        }
        else {
            bail!("Tree is not initialized. Call .init_tree() first")
        }
    }

    /// Retrieves the ID of the direct parent of a given entry in the tree
    /// structure.
    ///
    /// # Arguments
    ///
    /// * `id` - The ID of the child entry.
    pub fn get_parent(
        &self,
        id: &EntryId,
    ) -> anyhow::Result<Option<EntryId>> {
        if let Some(tree) = self.tree.as_ref().left() {
            Ok(tree.get_parent(*id))
        }
        else {
            bail!("Tree is not initialized. Call .init_tree() first")
        }
    }

    pub fn genomic_query(
        &self,
        contig: &Contig,
    ) -> anyhow::Result<Vec<EntryId>> {
        self.interval_map
            .as_ref()
            .left()
            .map(|imap| {
                imap.find(contig)
                    .unwrap_or_default()
                    .into_iter()
                    .cloned()
                    .collect()
            })
            .ok_or(anyhow!(
                "Interval map is not initialized. Call .init_imap() first"
            ))
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
