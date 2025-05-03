use std::collections::HashMap;
use std::fmt::Debug;
use std::str::FromStr;
use std::sync::Arc;
use std::{f64, fmt};

use arcstr::ArcStr;
use itertools::Itertools;
use multimap::MultiMap;
use serde::de::{self, Deserializer, Visitor};
use serde::{Deserialize, Serialize};

use super::coords::Contig;
use super::enums::Strand;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct GffEntryAttributes {
    id:            Option<String>,
    name:          Option<Vec<String>>,
    alias:         Option<Vec<String>>,
    parent:        Option<Vec<String>>,
    target:        Option<Vec<Contig<Arc<str>, u32>>>,
    gap:           Option<Vec<String>>,
    derives_from:  Option<Vec<String>>,
    note:          Option<Vec<String>>,
    dbxref:        Option<Vec<String>>,
    ontology_term: Option<Vec<String>>,
    other:         HashMap<String, String>,
}

impl fmt::Display for GffEntryAttributes {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        let mut serialized = serde_json::to_string(self).unwrap();
        serialized.pop();
        serialized.remove(0);
        write!(f, "{}", serialized)
    }
}

impl FromStr for GffEntryAttributes {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let formatted = format!("\"{}\"", s);
        let res: Self = serde_json::from_str(&formatted)?;
        Ok(res)
    }
}


impl Serialize for GffEntryAttributes {
    fn serialize<S>(
        &self,
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer, {
        let name_string = self.name.as_ref().map(|v| v.join(","));
        let alias_string = self.alias.as_ref().map(|v| v.join(","));
        let parent_string = self
            .parent
            .as_ref()
            .map(|v| v.join(","));
        let target_string = self.target.as_ref().map(|v| {
            v.iter()
                .map(|c| format!("{} {} {}", c.seqname(), c.start(), c.end()))
                .join(",")
        });
        let gap_string = self.gap.as_ref().map(|v| v.join(","));
        let derives_from_string = self
            .derives_from
            .as_ref()
            .map(|v| v.join(","));
        let note_string = self.note.as_ref().map(|v| v.join(","));
        let dbxref_string = self
            .dbxref
            .as_ref()
            .map(|v| v.join(","));
        let ontology_term_string = self
            .ontology_term
            .as_ref()
            .map(|v| v.join(","));
        let other_string = self
            .other
            .iter()
            .map(|(k, v)| format!("{}={}", k, v))
            .collect::<Vec<_>>()
            .join(";");

        let mut serialized = String::new();
        if let Some(id) = self.id.as_ref() {
            serialized.push_str(&format!("ID={}", id));
        }
        if let Some(name) = name_string {
            serialized.push_str(&format!(";Name={}", name));
        }
        if let Some(alias) = alias_string {
            serialized.push_str(&format!(";Alias={}", alias));
        }
        if let Some(parent) = parent_string {
            serialized.push_str(&format!(";Parent={}", parent));
        }
        if let Some(target) = target_string {
            serialized.push_str(&format!(";Target={}", target));
        }
        if let Some(gap) = gap_string {
            serialized.push_str(&format!(";Gap={}", gap));
        }
        if let Some(derives_from) = derives_from_string {
            serialized.push_str(&format!(";Derives_from={}", derives_from));
        }
        if let Some(note) = note_string {
            serialized.push_str(&format!(";Note={}", note));
        }
        if let Some(dbxref) = dbxref_string {
            serialized.push_str(&format!(";Dbxref={}", dbxref));
        }
        if let Some(ontology_term) = ontology_term_string {
            serialized.push_str(&format!(";Ontology_term={}", ontology_term));
        }
        if !other_string.is_empty() {
            serialized.push_str(&format!(";{}", other_string));
        }

        serializer.serialize_str(&serialized)
    }
}

impl<'de> serde::Deserialize<'de> for GffEntryAttributes {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>, {
        struct GffEntryAttributesVisitor;

        impl Visitor<'_> for GffEntryAttributesVisitor {
            type Value = GffEntryAttributes;

            fn expecting(
                &self,
                formatter: &mut fmt::Formatter,
            ) -> fmt::Result {
                formatter.write_str("a string of GFF attributes")
            }

            fn visit_str<E>(
                self,
                value: &str,
            ) -> Result<GffEntryAttributes, E>
            where
                E: de::Error, {
                let mut attributes = GffEntryAttributes::default();
                for pair in value.split(';') {
                    if pair.is_empty() {
                        continue;
                    }

                    let mut parts = pair.splitn(2, '=');
                    let key = parts
                        .next()
                        .ok_or_else(|| E::custom("Missing key"))?;
                    let value = parts.next();

                    match key {
                        "ID" => {
                            attributes.id = value.map(|s| s.to_string());
                        },
                        "Name" => {
                            attributes.name = value.map(|s| {
                                s.split(',')
                                    .map(|s| s.to_string())
                                    .collect()
                            });
                        },
                        "Alias" => {
                            attributes.alias = value.map(|s| {
                                s.split(',')
                                    .map(|s| s.to_string())
                                    .collect()
                            });
                        },
                        "Parent" => {
                            attributes.parent = value.map(|s| {
                                s.split(',')
                                    .map(|s| s.to_string())
                                    .collect()
                            });
                        },
                        "Target" => {
                            // TODO
                            attributes.target = None; // Setting to None for now
                                                      // because contig parsing
                                                      // is unimplemented
                        },
                        "Gap" => {
                            attributes.gap = value.map(|s| {
                                s.split(',')
                                    .map(|s| s.to_string())
                                    .collect()
                            });
                        },
                        "Derives_from" => {
                            attributes.derives_from = value.map(|s| {
                                s.split(',')
                                    .map(|s| s.to_string())
                                    .collect()
                            });
                        },
                        "Note" => {
                            attributes.note = value.map(|s| {
                                s.split(',')
                                    .map(|s| s.to_string())
                                    .collect()
                            });
                        },
                        "Dbxref" => {
                            attributes.dbxref = value.map(|s| {
                                s.split(',')
                                    .map(|s| s.to_string())
                                    .collect()
                            });
                        },
                        "Ontology_term" => {
                            attributes.ontology_term = value.map(|s| {
                                s.split(',')
                                    .map(|s| s.to_string())
                                    .collect()
                            });
                        },
                        _ => {
                            if let Some(val) = value {
                                attributes
                                    .other
                                    .insert(key.to_string(), val.to_string());
                            }
                        },
                    }
                }

                Ok(attributes)
            }
        }

        deserializer.deserialize_str(GffEntryAttributesVisitor)
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RawGffEntry {
    seqid:        String,
    source:       String,
    feature_type: String,
    start:        u32,
    end:          u32,
    score:        Option<f64>,
    strand:       char,
    phase:        Option<u8>,
    attributes:   String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GffEntry {
    contig:       Contig<ArcStr, u32>,
    source:       ArcStr,
    feature_type: ArcStr,
    score:        Option<f64>,
    phase:        Option<u8>,
    attributes:   GffEntryAttributes,
    id:           String,
}

impl From<bio::io::bed::Record> for GffEntry {
    fn from(value: bio::io::bed::Record) -> Self {
        let contig = Contig::from(value.clone());
        let contig = Contig::new(
            ArcStr::from(contig.seqname()),
            contig.start(),
            contig.end(),
            contig.strand(),
        );
        Self::new(
            contig,
            None,
            None,
            value
                .score()
                .map(|s| s.parse::<f64>().unwrap_or(f64::NAN)),
            None,
            None,
        )
    }
}

impl GffEntry {
    pub fn new(
        contig: Contig<ArcStr, u32>,
        source: Option<ArcStr>,
        feature_type: Option<ArcStr>,
        score: Option<f64>,
        phase: Option<u8>,
        attributes: Option<GffEntryAttributes>,
    ) -> Self {
        let attributes = attributes.unwrap_or_default();
        let id = attributes
            .id
            .as_ref()
            .cloned()
            .unwrap_or(uuid::Uuid::default().to_string());
        Self {
            contig,
            score,
            phase,
            id,
            attributes,
            source: source.unwrap_or_default(),
            feature_type: feature_type.unwrap_or_default(),
        }
    }
}

impl TryFrom<RawGffEntry> for GffEntry {
    type Error = anyhow::Error;

    fn try_from(value: RawGffEntry) -> Result<Self, Self::Error> {
        let seqid = ArcStr::from(value.seqid);
        let source = ArcStr::from(value.source);
        let strand = Strand::from(value.strand.to_string());
        let feature_type = ArcStr::from(value.feature_type);
        let attributes =
            GffEntryAttributes::from_str(value.attributes.as_str())?;

        Ok(GffEntry::new(
            Contig::new(seqid, value.start, value.end, strand),
            Some(source),
            Some(feature_type),
            value.score,
            value.phase,
            Some(attributes),
        ))
    }
}

impl TryFrom<GffEntry> for RawGffEntry {
    type Error = serde_json::Error;

    fn try_from(value: GffEntry) -> Result<Self, Self::Error> {
        let attributes = value.attributes.to_string();
        Ok(RawGffEntry {
            seqid: value.contig.seqname().to_string(),
            source: value.source.to_string(),
            feature_type: value.feature_type.to_string(),
            start: value.contig.start(),
            end: value.contig.end(),
            score: value.score,
            strand: value.contig.strand().into(),
            phase: value.phase,
            attributes,
        })
    }
}

pub struct AnnotStore {
    id_map:       HashMap<String, GffEntry>,
    parent_map:   MultiMap<String, String>,
    children_map: MultiMap<String, String>,
}

impl Default for AnnotStore {
    fn default() -> Self { Self::new() }
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
        id: &String,
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
        id: &String,
    ) -> Option<&Vec<String>> {
        self.parent_map.get_vec(id)
    }

    pub fn get_children(
        &self,
        id: &String,
    ) -> Option<&Vec<String>> {
        self.children_map.get_vec(id)
    }
}

pub struct AnnotIterator<'a> {
    store: &'a AnnotStore,
    iter:  std::collections::hash_map::Iter<'a, String, GffEntry>,
}

impl<'a> Iterator for AnnotIterator<'a> {
    type Item = (&'a String, &'a GffEntry);

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
}

#[cfg(test)]
mod tests {
    use arcstr::ArcStr;

    use super::*;

    #[test]
    fn test_gff_entry_attributes_serialization() {
        let mut attributes = GffEntryAttributes::default();
        attributes.id = Some("gene123".to_string());
        attributes.name = Some(vec!["my_gene".to_string()]);
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
        expected_attributes.id = Some("gene123".to_string());
        expected_attributes.name = Some(vec!["my_gene".to_string()]);
        expected_attributes.alias = Some(vec!["another_name".to_string()]);
        expected_attributes
            .other
            .insert("custom".to_string(), "value".to_string());

        assert_eq!(deserialized, expected_attributes);
    }

    #[test]
    fn test_raw_gff_entry_conversion() {
        let raw_gff_entry = RawGffEntry {
            seqid:        "chr1".to_string(),
            source:       "test".to_string(),
            feature_type: "gene".to_string(),
            start:        100,
            end:          200,
            score:        Some(0.9),
            strand:       '+',
            phase:        Some(0),
            attributes:   "ID=gene123".to_string(),
        };

        let gff_entry = GffEntry::try_from(raw_gff_entry).unwrap();

        assert_eq!(gff_entry.contig.seqname(), ArcStr::from("chr1"));
        assert_eq!(gff_entry.contig.start(), 100);
        assert_eq!(gff_entry.contig.end(), 200);
        assert_eq!(gff_entry.source, ArcStr::from("test"));
        assert_eq!(gff_entry.feature_type, ArcStr::from("gene"));
        assert_eq!(gff_entry.score, Some(0.9));
        assert_eq!(gff_entry.phase, Some(0));
        assert_eq!(gff_entry.attributes.id, Some("gene123".to_string()));
    }

    #[test]
    fn test_annot_store_insert_and_retrieve() {
        let mut store = AnnotStore::new();

        let contig = Contig::new(ArcStr::from("chr1"), 1, 100, Strand::Forward);
        let attributes = GffEntryAttributes {
            id:            Some("gene1".to_string()),
            name:          Some(vec!["my_gene".to_string()]),
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
        let retrieved_entry = store.id_map.get("gene1").unwrap();
        assert_eq!(retrieved_entry, &entry);
    }

    #[test]
    fn test_annot_store_parent_child_relationships() {
        let mut store = AnnotStore::new();

        let contig1 =
            Contig::new(ArcStr::from("chr1"), 1, 100, Strand::Forward);
        let attributes1 = GffEntryAttributes {
            id:            Some("gene1".to_string()),
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
            id:            Some("exon1".to_string()),
            name:          None,
            alias:         None,
            parent:        Some(vec!["gene1".to_string()]),
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
            .get_parents(&"exon1".to_string())
            .unwrap();
        assert_eq!(parents, &vec!["gene1".to_string()]);

        let children = store
            .get_children(&"gene1".to_string())
            .unwrap();
        assert_eq!(children, &vec!["exon1".to_string()]);
    }
}
