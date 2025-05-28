use std::fs::File;
use std::path::PathBuf;
use std::str::FromStr;

use arcstr::ArcStr;

use super::*;
use crate::prelude::Strand;
use crate::{data_structs::annotation::gff_entry::RawGffEntry, prelude::Contig};
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
    let gff_string = "ID=gene123;Name=my_gene;Alias=another_name,other_name;custom=value;";
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

#[test]
fn test_basic_construction() {
    HcAnnotStore::new();
    HcAnnotStore::default();

    let entries = vec![
        GffEntry::new(Contig::new("123".into(), 123, 456, Strand::None), None, None, None, None, None),
        GffEntry::new(Contig::new("123".into(), 345, 567, Strand::None), None, None, None, None, None)
    ];
    HcAnnotStore::from_iter(entries);
}

#[test]
fn test_with_capacity() {
    let capacity = 1000;
    let annot_store = HcAnnotStore::with_capacity(capacity);
    assert!(annot_store.entries().capacity() >= capacity);
    assert!(annot_store.tree().capacity() >= capacity);
    assert!(annot_store.tree_ids().capacity() >= capacity);
}

use rstest::{fixture, rstest};

#[fixture]
fn bed_store() -> HcAnnotStore {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/genes.bed");
    HcAnnotStore::from_bed(File::open(path).unwrap()).unwrap()
}

#[fixture]
#[once]
fn gff_store_static() -> HcAnnotStore {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/annot.gff");
    HcAnnotStore::from_gff(File::open(path).unwrap()).unwrap()
}

#[fixture]
#[once]
fn bed_store_static() -> HcAnnotStore {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/genes.bed");
    HcAnnotStore::from_bed(File::open(path).unwrap()).unwrap()
}

#[rstest]
fn test_read_bed(bed_store_static: &HcAnnotStore) {
    assert!(!bed_store_static.is_empty())
}

#[rstest]
fn test_read_gff(gff_store_static: &HcAnnotStore) {
    assert!(!gff_store_static.is_empty())
}

#[rstest]
fn test_search_regex(gff_store_static: &HcAnnotStore) {
    let res = gff_store_static.get_entries_regex(r"gene-AT1G010\d0").unwrap();
    assert!(!res.is_empty());
    assert_eq!(res.len(), 9);
}

#[rstest]
fn test_get_feature_types(bed_store_static: &HcAnnotStore) {
    assert!(!bed_store_static.get_feature_types().is_empty())
}

#[rstest]
fn test_add_flank(bed_store: HcAnnotStore) {
    let mut store = bed_store;
    let initial_count = store.entries().len();
    assert!(initial_count > 0, "Test BED file should not be empty");

    let flank_len = 50;
    let prefix = "flank_";

    // Add flanks (positive and negative)
    store.add_flank(|_| true, flank_len, prefix);
    let count_after_flank = store.entries().len();
    // Should add two flanks (positive and negative) for every original entry
    assert_eq!(count_after_flank, initial_count * 2, "Should add two flanks for every entry");

    // Find an original entry (one whose feature type doesn't start with the prefix)
    let original_entry = store.entries().values().find(|e| !e.feature_type.starts_with(prefix)).expect("Should find an original entry");
    let original_id = original_entry.id.clone(); // Clone id for easier use

    // --- Test get_children ---
    let children_ids = store.get_children(&original_id.as_str().into()).expect("Original entry should have children");
    assert_eq!(children_ids.len(), 1, "Original entry should have two flank children");

    // Identify the positive and negative flank IDs among the children
    let expected_pos_flank_id_str = format!("{}_flank_{}", original_id, flank_len);

    let pos_flank_id = children_ids.iter()
        .find(|id| id.as_str() == expected_pos_flank_id_str)
        .expect("Positive flank ID should be in children list");

    // Ensure both found IDs are distinct and match the expected formats
    assert_eq!(children_ids.contains(pos_flank_id), true);


    // --- Test get_parent for children ---
    let pos_parent_id = store.get_parent(pos_flank_id).expect("Positive flank should have a parent");
    assert_eq!(pos_parent_id, ArcStr::from(original_id.as_str()), "Positive flank's parent should be the original entry");

    // Test get_parent for a root node (original entry)
    let original_parent = store.get_parent(&original_id.as_str().into());
    assert_eq!(original_parent, Some("BSX_ROOT_NODE".into()), "Original entry should be a child of root node");

    // Test get_children for a leaf node (a flank entry)
    let pos_flank_children = store.get_children(pos_flank_id).expect("Positive flank entry should be found");
    assert!(pos_flank_children.is_empty(), "Flank entry should be a leaf node (no children)");

    // --- Verify properties of flank entry (checking positive flank as in original test) ---
    let pos_flank_entry = store.get_entry(pos_flank_id).expect("Positive flank entry should exist");

    // Verify the flank entry properties
    assert_eq!(pos_flank_entry.feature_type, format!("{}{}", prefix, original_entry.feature_type));
    // Check that the parent attribute list contains the original ID
    assert!(pos_flank_entry.attributes.parent.as_ref().map_or(false, |parents| parents.iter().any(|p| p.as_str() == original_id.as_str())));
    assert_eq!(pos_flank_entry.contig.seqname(), original_entry.contig.seqname());

    // Verify position for positive flank based on current code behavior: Start = original end, End = original end + flank_len
    assert_eq!(pos_flank_entry.contig.start(), original_entry.contig.end(), "Positive flank should start at original end");
    assert_eq!(pos_flank_entry.contig.end(), original_entry.contig.end() + flank_len as u32, "Positive flank should end at original end + flank_len");
}

#[rstest]
fn test_genomic_query(bed_store_static: &HcAnnotStore) {
    let res = bed_store_static.genomic_query(&Contig::new("NC_003070.9".into(), 0, 10000, Strand::None));
    assert!(res.is_some());
    assert_eq!(res.as_ref().unwrap().len(), 2)
}
