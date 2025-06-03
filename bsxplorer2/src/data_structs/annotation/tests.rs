use std::fs::File;
use std::path::PathBuf;
use std::str::FromStr;
use std::time::Instant;

use arcstr::ArcStr;

use super::*;
use crate::data_structs::annotation::gff_entry::RawGffEntry;
use crate::data_structs::typedef::BsxSmallStr;
use crate::prelude::{
    Contig,
    Strand,
};

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
        GffEntry::new(
            Contig::new("123".into(), 123, 456, Strand::None),
            None,
            None,
            None,
            None,
            None,
        ),
        GffEntry::new(
            Contig::new("123".into(), 345, 567, Strand::None),
            None,
            None,
            None,
            None,
            None,
        ),
    ];
    HcAnnotStore::from_iter(entries);
}

use rstest::{
    fixture,
    rstest,
};

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
fn test_add_flanks(mut bed_store: HcAnnotStore) {
    let start = Instant::now();
    bed_store.add_flank(|_| true, 2000, "upstream_");
    println!("Elapsed: {}", start.elapsed().as_millis())
}

#[rstest]
fn test_read_gff(gff_store_static: &HcAnnotStore) {
    assert!(!gff_store_static.is_empty())
}

#[rstest]
fn test_search_regex(gff_store_static: &HcAnnotStore) {
    let res = gff_store_static
        .get_entries_regex(r"gene-AT1G010\d0")
        .unwrap();
    assert!(!res.is_empty());
    assert_eq!(res.len(), 9);
}

#[rstest]
fn test_get_feature_types(bed_store_static: &HcAnnotStore) {
    assert!(!bed_store_static.get_feature_types().is_empty())
}

#[rstest]
fn test_genomic_query(mut bed_store: HcAnnotStore) {
    bed_store.init_imap();
    let res = bed_store
        .genomic_query(&Contig::new("NC_003070.9".into(), 0, 10000, Strand::None))
        .unwrap();
    assert!(!res.is_empty());
    assert_eq!(res.len(), 2)
}

#[rstest]
fn test_init_tree(mut bed_store: HcAnnotStore) {
    bed_store.init_tree().unwrap()
}

#[rstest]
fn test_init_imap(mut bed_store: HcAnnotStore) {
    bed_store.init_imap()
}
