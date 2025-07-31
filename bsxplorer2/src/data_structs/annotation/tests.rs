use std::str::FromStr;

use arcstr::ArcStr;
use rstest::{fixture, rstest};

use super::*;
use super::annot_store::EntryTree;
use crate::data_structs::annotation::gff_entry::RawGffEntry;
use crate::data_structs::typedef::BsxSmallStr;

#[test]
fn test_gff_entry_attributes_serialization() {
    let mut attributes = GffEntryAttributes::default();
    attributes.id = Some("gene123".into());
    attributes.name = Some(vec!["my_gene".into()]);
    attributes
        .other
        .insert("custom".to_string(), "value".to_string());

    let serialized = attributes.to_string();
    assert_eq!(serialized, "ID=gene123;Name=my_gene;custom=value");
}

#[test]
fn test_gff_entry_attributes_deserialization() {
    let gff_string =
        "ID=gene123;Name=my_gene;Alias=another_name,other_name;custom=value;";
    let deserialized = GffEntryAttributes::from_str(gff_string).unwrap();

    let mut expected_attributes = GffEntryAttributes::default();
    expected_attributes.id = Some("gene123".into());
    expected_attributes.name = Some(vec!["my_gene".into()]);
    expected_attributes.alias = Some(vec!["another_name".into(), "other_name".into()]);
    expected_attributes
        .other
        .insert("custom".to_string(), "value".to_string());

    assert_eq!(deserialized, expected_attributes);
}

#[test]
fn test_raw_gff_entry_conversion() {
    let raw_gff_entry = RawGffEntry {
        seqid:        "chr1".into(),
        source:       "test".into(),
        feature_type: "gene".into(),
        start:        100,
        end:          200,
        score:        Some(0.9),
        strand:       '+',
        phase:        Some(0),
        attributes:   "ID=gene123".to_string(),
    };

    let gff_entry = GffEntry::try_from(raw_gff_entry).unwrap();

    assert_eq!(
        gff_entry.contig.seqname().to_owned(),
        BsxSmallStr::from_str("chr1").unwrap()
    );
    assert_eq!(gff_entry.contig.start(), 100);
    assert_eq!(gff_entry.contig.end(), 200);
    assert_eq!(gff_entry.source, ArcStr::from("test"));
    assert_eq!(gff_entry.feature_type, ArcStr::from("gene"));
    assert_eq!(gff_entry.score, Some(0.9));
    assert_eq!(gff_entry.phase, Some(0));
    assert_eq!(gff_entry.attributes.id, Some("gene123".into()));
}

#[fixture]
fn entry_tree_entries() -> Vec<(EntryId, Option<EntryId>)> {
    vec![
        (1.into(), None),
        (2.into(), Some(1.into())),
        (3.into(), Some(1.into())),
        (4.into(), Some(2.into())),
        (5.into(), Some(2.into())),
        (6.into(), Some(3.into())),
    ]
}

#[fixture]
fn simple_entry_tree(entry_tree_entries: Vec<(EntryId, Option<EntryId>)>) -> EntryTree {
    let mut tree = EntryTree::new();
    tree.append(entry_tree_entries).unwrap();
    tree
}

#[rstest]
fn test_entry_tree_from_iter(entry_tree_entries: Vec<(EntryId, Option<EntryId>)>) {
    let tree = EntryTree::from_iter(entry_tree_entries);

    // Verify parents
    assert_eq!(tree.get_parent(1), Some(u64::MAX.into())); // Root node is child of internal root
    assert_eq!(tree.get_parent(2), Some(1.into()));
    assert_eq!(tree.get_parent(3), Some(1.into()));
    assert_eq!(tree.get_parent(4), Some(2.into()));
    assert_eq!(tree.get_parent(5), Some(2.into()));
    assert_eq!(tree.get_parent(6), Some(3.into()));
    assert_eq!(tree.get_parent(99), None); // Non-existent entry

    // Verify children
    let mut children_of_1 = tree.get_children(1).unwrap();
    children_of_1.sort();
    assert_eq!(children_of_1, vec![2.into(), 3.into()]);

    let mut children_of_2 = tree.get_children(2).unwrap();
    children_of_2.sort();
    assert_eq!(children_of_2, vec![4.into(), 5.into()]);

    let mut children_of_3 = tree.get_children(3).unwrap();
    children_of_3.sort();
    assert_eq!(children_of_3, vec![6.into()]);

    assert!(tree.get_children(4).unwrap().is_empty()); // Leaf node
    assert!(tree.get_children(99).is_none());          // Non-existent entry

    // Verify children of the internal root
    let mut children_of_internal_root = tree.get_children(u64::MAX).unwrap();
    children_of_internal_root.sort();
    assert_eq!(children_of_internal_root, vec![1.into()]);
}

#[rstest]
fn test_entry_tree_append_existing_parents(mut simple_entry_tree: EntryTree) {
    // Append new nodes with existing parents
    let new_entries = vec![
        (7.into(), Some(1.into())), // Child of 1
        (8.into(), Some(4.into())), // Child of 4 (leaf node)
    ];
    simple_entry_tree.append(new_entries).unwrap();

    // Verify new parent relationships
    assert_eq!(simple_entry_tree.get_parent(7), Some(1.into()));
    assert_eq!(simple_entry_tree.get_parent(8), Some(4.into()));

    // Verify children of parent 1
    let mut children_of_1 = simple_entry_tree.get_children(1).unwrap();
    children_of_1.sort();
    assert_eq!(children_of_1, vec![2.into(), 3.into(), 7.into()]);

    // Verify children of parent 4
    let mut children_of_4 = simple_entry_tree.get_children(4).unwrap();
    children_of_4.sort();
    assert_eq!(children_of_4, vec![8.into()]);

    // Append a new root node (which becomes a child of EntryId::MAX)
    simple_entry_tree.append(vec![(9.into(), None)]).unwrap();
    assert_eq!(simple_entry_tree.get_parent(9), Some(u64::MAX.into()));

    let mut children_of_internal_root = simple_entry_tree.get_children(u64::MAX).unwrap();
    children_of_internal_root.sort();
    assert!(children_of_internal_root.contains(&1.into()));
    assert!(children_of_internal_root.contains(&9.into()));
}

#[rstest]
fn test_entry_tree_append_unexistent_parents(mut simple_entry_tree: EntryTree) {
    // Attempt to append with an unexistent parent
    let new_entries = vec![
        (10.into(), Some(99.into())), // 99 does not exist
    ];
    let result = simple_entry_tree.append(new_entries);
    assert!(result.is_err());
    assert_eq!(
        result.unwrap_err().to_string(),
        "Some children have unexistent parents"
    );

    // Ensure the tree state hasn't changed from the original (no 10 added)
    assert_eq!(simple_entry_tree.get_parent(10), None);
    assert_eq!(simple_entry_tree.get_children(10), None);
}

#[rstest]
fn test_entry_tree_insert_under() {
    let mut tree = EntryTree::new();
    tree.insert_to_root(1).unwrap(); // Add a root node (child of EntryId::MAX)

    // Insert 2 under 1
    tree.insert_under(2, 1).unwrap();
    assert_eq!(tree.get_parent(2), Some(1.into()));
    let mut children_of_1 = tree.get_children(1).unwrap();
    children_of_1.sort();
    assert_eq!(children_of_1, vec![2.into()]);

    // Try to insert under non-existent parent
    let _ = tree.insert_under(3, 99);

    // Insert 4 under 2
    tree.insert_under(4, 2).unwrap();
    assert_eq!(tree.get_parent(4), Some(2.into()));
    let mut children_of_2 = tree.get_children(2).unwrap();
    children_of_2.sort();
    assert_eq!(children_of_2, vec![4.into()]);

    // Try to insert a child that already has a different parent
    // First, insert 5 as a child of 1
    tree.insert_under(5, 1).unwrap();
    assert_eq!(tree.get_parent(5), Some(1.into()));
    // Now try to insert 5 under 2, should fail
    let result = tree.insert_under(5, 2);
    assert!(result.is_err());

    // Try to insert a child that already has the same parent
    let result = tree.insert_under(2, 1);
    assert!(result.is_ok()); // Should be Ok(())
    assert_eq!(tree.get_parent(2), Some(1.into())); // No change in parent
}

#[rstest]
fn test_entry_tree_remove(mut simple_entry_tree: EntryTree) {
    // Before removal:
    // EntryId::MAX -> 1
    // 1 -> 2, 3
    // 2 -> 4, 5
    // 3 -> 6

    // Remove a leaf node (4)
    simple_entry_tree.remove(4).unwrap();
    assert_eq!(simple_entry_tree.get_parent(4), None); // Should be gone
    let children_of_2 = simple_entry_tree.get_children(2).unwrap();
    assert_eq!(children_of_2, vec![5.into()]); // Parent 2 should no longer have 4 as child

    // Remove a node with children (2)
    simple_entry_tree.remove(2).unwrap();
    assert_eq!(simple_entry_tree.get_parent(2), None); // Node 2 should be gone
    assert_eq!(simple_entry_tree.get_parent(5), None); // Child 5 should also be gone due to DropChildren
    let mut children_of_1 = simple_entry_tree.get_children(1).unwrap();
    children_of_1.sort();
    assert_eq!(children_of_1, vec![3.into()]); // Parent 1 should no longer have 2 as child

    // Try to remove non-existent id
    let result = simple_entry_tree.remove(99);
    assert!(result.is_err());
    assert_eq!(
        result.unwrap_err().to_string(),
        "Such id did not exist in a tree"
    );

    // Remove the user-defined root (1)
    simple_entry_tree.remove(1).unwrap();
    assert_eq!(simple_entry_tree.get_parent(1), None); // Node 1 should be gone
    assert_eq!(simple_entry_tree.get_parent(3), None); // Child 3 should also be gone
    assert_eq!(simple_entry_tree.get_parent(6), None); // Child 6 should also be gone

    // After removing user-defined root 1, the children of the internal `EntryId::MAX` root
    // should be empty as 1 was its only child.
    let children_of_internal_root = simple_entry_tree.get_children(u64::MAX).unwrap();
    assert!(children_of_internal_root.is_empty());
}

#[rstest]
fn test_entry_tree_get_parent(simple_entry_tree: EntryTree) {
    // Existing parents
    assert_eq!(simple_entry_tree.get_parent(2), Some(1.into()));
    assert_eq!(simple_entry_tree.get_parent(4), Some(2.into()));
    assert_eq!(simple_entry_tree.get_parent(6), Some(3.into()));

    // User-defined root node (parent is the internal EntryId::MAX root)
    assert_eq!(simple_entry_tree.get_parent(1), Some(u64::MAX.into()));

    // Non-existent node
    assert_eq!(simple_entry_tree.get_parent(99), None);

    // Leaf node
    assert_eq!(simple_entry_tree.get_parent(4), Some(2.into()));
}

#[rstest]
fn test_entry_tree_get_children(simple_entry_tree: EntryTree) {
    // Node with children
    let mut children_of_1 = simple_entry_tree.get_children(1).unwrap();
    children_of_1.sort();
    assert_eq!(children_of_1, vec![2.into(), 3.into()]);

    let mut children_of_2 = simple_entry_tree.get_children(2).unwrap();
    children_of_2.sort();
    assert_eq!(children_of_2, vec![4.into(), 5.into()]);

    // Leaf node (no children)
    assert!(simple_entry_tree.get_children(4).unwrap().is_empty());
    assert!(simple_entry_tree.get_children(5).unwrap().is_empty());
    assert!(simple_entry_tree.get_children(6).unwrap().is_empty());

    // Non-existent node
    assert_eq!(simple_entry_tree.get_children(99), None);

    // Get children of the implicit root EntryId::MAX
    let mut children_of_max_root = simple_entry_tree.get_children(u64::MAX).unwrap();
    children_of_max_root.sort();
    assert_eq!(children_of_max_root, vec![1.into()]);
}
