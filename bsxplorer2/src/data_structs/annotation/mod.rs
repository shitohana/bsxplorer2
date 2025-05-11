pub mod annot_store;
pub mod gff_entry;

pub use annot_store::AnnotStore;
pub use gff_entry::{GffEntry, GffEntryAttributes};

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use arcstr::ArcStr;
    use hashbrown::HashMap;

    use super::*;
    use crate::data_structs::annotation::gff_entry::RawGffEntry;
    use crate::data_structs::coords::Contig;
    use crate::data_structs::enums::Strand;

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
            "ID=gene123;Name=my_gene;Alias=another_name;custom=value";
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

        assert_eq!(gff_entry.contig.seqname().to_owned(), ArcStr::from("chr1"));
        assert_eq!(gff_entry.contig.start(), 100);
        assert_eq!(gff_entry.contig.end(), 200);
        assert_eq!(gff_entry.source, ArcStr::from("test"));
        assert_eq!(gff_entry.feature_type, ArcStr::from("gene"));
        assert_eq!(gff_entry.score, Some(0.9));
        assert_eq!(gff_entry.phase, Some(0));
        assert_eq!(gff_entry.attributes.id, Some("gene123".into()));
    }

    #[test]
    fn test_annot_store_insert_and_retrieve() {
        let mut store = AnnotStore::new();

        let contig = Contig::new(ArcStr::from("chr1"), 1, 100, Strand::Forward);
        let attributes = GffEntryAttributes {
            id:            Some("gene1".into()),
            name:          Some(vec!["my_gene".into()]),
            alias:         None,
            parent:        None,
            target:        None,
            gap:           None,
            derives_from:  None,
            note:          None,
            dbxref:        None,
            ontology_term: None,
            other:         HashMap::new(),
        };
        let entry = GffEntry::new(
            contig,
            Some(ArcStr::from("test_source")),
            Some(ArcStr::from("gene")),
            Some(0.5),
            Some(1),
            Some(attributes),
        );

        store.insert(entry.clone());
        let retrieved_entry = store.id_map().get("gene1").unwrap();
        assert_eq!(retrieved_entry, &entry);
    }

    #[test]
    fn test_annot_store_parent_child_relationships() {
        let mut store = AnnotStore::new();

        let contig1 =
            Contig::new(ArcStr::from("chr1"), 1, 100, Strand::Forward);
        let attributes1 = GffEntryAttributes {
            id:            Some("gene1".into()),
            name:          None,
            alias:         None,
            parent:        None,
            target:        None,
            gap:           None,
            derives_from:  None,
            note:          None,
            dbxref:        None,
            ontology_term: None,
            other:         HashMap::new(),
        };
        let entry1 = GffEntry::new(
            contig1,
            Some(ArcStr::from("test_source")),
            Some(ArcStr::from("gene")),
            Some(0.5),
            Some(1),
            Some(attributes1),
        );

        let contig2 =
            Contig::new(ArcStr::from("chr1"), 10, 50, Strand::Forward);
        let attributes2 = GffEntryAttributes {
            id:            Some("exon1".into()),
            name:          None,
            alias:         None,
            parent:        Some(vec!["gene1".into()]),
            target:        None,
            gap:           None,
            derives_from:  None,
            note:          None,
            dbxref:        None,
            ontology_term: None,
            other:         HashMap::new(),
        };
        let entry2 = GffEntry::new(
            contig2,
            Some(ArcStr::from("test_source")),
            Some(ArcStr::from("exon")),
            Some(0.7),
            Some(2),
            Some(attributes2),
        );

        store.insert(entry1);
        store.insert(entry2);

        let parents = store
            .get_parents(&"exon1".into())
            .unwrap();
        assert_eq!(parents, &vec!["gene1".to_string()]);

        let children = store
            .get_children(&"gene1".into())
            .unwrap();
        assert_eq!(children, &vec!["exon1".to_string()]);
    }
}
