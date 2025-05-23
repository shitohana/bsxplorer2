use std::str::FromStr;

use arcstr::ArcStr;

use super::*;
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
    let gff_string = "ID=gene123;Name=my_gene;Alias=another_name;custom=value";
    let deserialized = GffEntryAttributes::from_str(gff_string).unwrap();

    let mut expected_attributes = GffEntryAttributes::default();
    expected_attributes.id = Some("gene123".into());
    expected_attributes.name = Some(vec!["my_gene".into()]);
    expected_attributes.alias = Some(vec!["another_name".into()]);
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
        BsxSmallStr::from_str("chr1")
    );
    assert_eq!(gff_entry.contig.start(), 100);
    assert_eq!(gff_entry.contig.end(), 200);
    assert_eq!(gff_entry.source, ArcStr::from("test"));
    assert_eq!(gff_entry.feature_type, ArcStr::from("gene"));
    assert_eq!(gff_entry.score, Some(0.9));
    assert_eq!(gff_entry.phase, Some(0));
    assert_eq!(gff_entry.attributes.id, Some("gene123".into()));
}

/// Helper function to create a test GffEntry.
///
/// This function is not meant to be used directly but is used by the
/// doctests.
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

use id_tree::InsertBehavior::*;
use id_tree::{Node, NodeId, Tree, TreeBuilder};
use itertools::Itertools;
use crate::data_structs::coords::Contig;
use crate::data_structs::enums::Strand;
use crate::io::bsx::BatchIndex;

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