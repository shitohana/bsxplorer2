use std::fmt::{
    Debug,
    Write,
};
use std::hash::Hash;
use std::str::FromStr;
use std::{
    f64,
    fmt,
};

use anyhow::anyhow;
use arcstr::ArcStr;
use hashbrown::HashMap;
use nanoid::nanoid;
use paste::paste;
use serde::de::{
    self,
    Deserializer,
    Visitor,
};
use serde::{
    Deserialize,
    Serialize,
};

use crate::data_structs::coords::Contig;
use crate::data_structs::enums::Strand;
use crate::data_structs::typedef::BsxSmallStr;
use crate::{
    getter_fn,
    with_field_fn,
};

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct GffEntryAttributes {
    pub id:            Option<BsxSmallStr>,
    pub name:          Option<Vec<BsxSmallStr>>,
    pub alias:         Option<Vec<BsxSmallStr>>,
    pub parent:        Option<Vec<BsxSmallStr>>,
    pub target:        Option<Vec<Contig>>,
    pub gap:           Option<Vec<String>>,
    pub derives_from:  Option<Vec<String>>,
    pub note:          Option<Vec<String>>,
    pub dbxref:        Option<Vec<BsxSmallStr>>,
    pub ontology_term: Option<Vec<String>>,
    pub other:         HashMap<String, String>,
}

impl Hash for GffEntryAttributes {
    fn hash<H: std::hash::Hasher>(
        &self,
        state: &mut H,
    ) {
        self.id.hash(state);
        // 'other' field is intentionally excluded from hashing
    }
}

macro_rules! with_vec_field_fn {
    ($name: ident, $item: ty, $op: expr) => {
        paste! {
            #[cfg_attr(coverage_nightly, coverage(off))]
            pub fn [<with_$name>]<S: $item, I: IntoIterator<Item = S>>(
                mut self,
                values: I,
            ) -> Self {
                self.$name = Some(values.into_iter().map($op).collect());
                self
            }
        }
    };

    ($name: ident, $type: ty) => {
        paste! {
            #[cfg_attr(coverage_nightly, coverage(off))]
            pub fn [<with_$name>]<I: IntoIterator<Item = $type>>(
                mut self,
                values: I,
            ) -> Self {
                self.$name = Some(values.into_iter().map(|s| s.to_owned()).collect());
                self
            }
        }
    };
}

impl GffEntryAttributes {
    with_vec_field_fn!(target, Contig);

    with_vec_field_fn!(gap, String);

    with_vec_field_fn!(derives_from, String);

    with_vec_field_fn!(note, String);

    with_vec_field_fn!(ontology_term, String);

    with_vec_field_fn!(dbxref, AsRef<str>, |v| v.as_ref().into());

    with_vec_field_fn!(name, AsRef<str>, |v| v.as_ref().into());

    with_vec_field_fn!(alias, AsRef<str>, |v| v.as_ref().into());

    with_vec_field_fn!(parent, AsRef<str>, |v| v.as_ref().into());

    with_field_fn!(other, HashMap<String, String>);

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn id(&self) -> Option<&BsxSmallStr> {
        self.id.as_ref()
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn name(&self) -> Option<&Vec<BsxSmallStr>> {
        self.name.as_ref()
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn alias(&self) -> Option<&Vec<BsxSmallStr>> {
        self.alias.as_ref()
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn parent(&self) -> Option<&Vec<BsxSmallStr>> {
        self.parent.as_ref()
    }

    /// Sets the ID attribute.
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn with_id<S: AsRef<str>>(
        mut self,
        id: S,
    ) -> Self {
        self.id = Some(id.as_ref().into());
        self
    }
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
        let mut attributes = GffEntryAttributes::default();
        for pair in s.split(';') {
            if pair.is_empty() {
                continue;
            }

            let mut parts = pair.splitn(2, '=');
            let key = parts.next().ok_or(anyhow!("Missing key"))?;
            let value = parts.next();

            match key {
                "ID" => {
                    attributes.id = value.map(|s| s.into());
                },
                "Name" => {
                    attributes.name =
                        value.map(|s| s.split(',').map(|s| s.into()).collect());
                },
                "Alias" => {
                    attributes.alias =
                        value.map(|s| s.split(',').map(|s| s.into()).collect());
                },
                "Parent" => {
                    attributes.parent =
                        value.map(|s| s.split(',').map(|s| s.into()).collect());
                },
                "Target" => {
                    // TODO
                    eprintln!("Target attribute parsing is not currently implemented");
                    attributes.target = None; // Setting to None for now
                                              // because contig parsing
                                              // is unimplemented
                },
                "Gap" => {
                    attributes.gap =
                        value.map(|s| s.split(',').map(|s| s.to_string()).collect());
                },
                "Derives_from" => {
                    attributes.derives_from =
                        value.map(|s| s.split(',').map(|s| s.to_string()).collect());
                },
                "Note" => {
                    attributes.note =
                        value.map(|s| s.split(',').map(|s| s.to_string()).collect());
                },
                "Dbxref" => {
                    attributes.dbxref =
                        value.map(|s| s.split(',').map(|s| s.into()).collect());
                },
                "Ontology_term" => {
                    attributes.ontology_term =
                        value.map(|s| s.split(',').map(|s| s.to_string()).collect());
                },
                _ => {
                    if let Some(val) = value {
                        attributes.other.insert(key.to_string(), val.to_string());
                    }
                },
            }
        }

        Ok(attributes)
    }
}

impl Serialize for GffEntryAttributes {
    fn serialize<S>(
        &self,
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer, {
        let mut serialized = String::with_capacity(128);
        let mut first = true;

        macro_rules! write_attr {
            // For fields that are Vec<T: Display>
            ($field:expr, $key:literal) => {
                if let Some(val) = $field.as_ref() {
                    if !first {
                        serialized.push(';');
                    }
                    else {
                        first = false;
                    }
                    write!(serialized, "{}=", $key)
                        .map_err(serde::ser::Error::custom)?;
                    let mut first_val = true;
                    for item in val {
                        if !first_val {
                            serialized.push(',');
                        }
                        else {
                            first_val = false;
                        }
                        write!(serialized, "{}", item)
                            .map_err(serde::ser::Error::custom)?;
                    }
                }
            };
            // For fields needing custom formatting (like Target)
            ($field:expr, $key:literal, $formatter:expr) => {
                if let Some(val) = $field.as_ref() {
                    if !first {
                        serialized.push(';');
                    }
                    else {
                        first = false;
                    }
                    write!(serialized, "{}=", $key)
                        .map_err(serde::ser::Error::custom)?;
                    let mut first_val = true;
                    for item in val {
                        if !first_val {
                            serialized.push(',');
                        }
                        else {
                            first_val = false;
                        }
                        // Apply the provided formatter closure
                        write!(serialized, "{}", $formatter(item))
                            .map_err(serde::ser::Error::custom)?;
                    }
                }
            };
        }
        if let Some(id) = self.id.as_ref() {
            // No need to check 'first' here, ID is always first if present
            write!(serialized, "ID={}", id).map_err(serde::ser::Error::custom)?;
            first = false; // Mark that we've written the first attribute
        }

        write_attr!(self.name, "Name");
        write_attr!(self.alias, "Alias");
        write_attr!(self.parent, "Parent");
        write_attr!(self.target, "Target", |c: &Contig| {
            format!("{} {} {}", c.seqname(), c.start(), c.end())
        });
        write_attr!(self.gap, "Gap");
        write_attr!(self.derives_from, "Derives_from");
        write_attr!(self.note, "Note");
        write_attr!(self.dbxref, "Dbxref");
        write_attr!(self.ontology_term, "Ontology_term");

        let mut sorted_other: Vec<_> = self.other.iter().collect();
        sorted_other.sort_unstable_by_key(|(k, _)| *k);

        for (k, v) in sorted_other {
            if !first {
                serialized.push(';');
            }
            else {
                first = false;
            }
            // Use write! for potentially better performance than format! +
            // push_str
            write!(serialized, "{}={}", k, v).map_err(serde::ser::Error::custom)?;
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
                GffEntryAttributes::from_str(value).map_err(serde::de::Error::custom)
            }
        }

        deserializer.deserialize_str(GffEntryAttributesVisitor)
    }
}
fn deserialize_optional_f64<'de, D>(deserializer: D) -> Result<Option<f64>, D::Error>
where
    D: serde::Deserializer<'de>, {
    let s = String::deserialize(deserializer)?;
    if s == "." {
        Ok(None)
    }
    else {
        s.parse::<f64>().map(Some).map_err(|e| {
            serde::de::Error::custom(format!("Failed to parse f64: {}", e))
        })
    }
}

fn deserialize_optional_u8<'de, D>(deserializer: D) -> Result<Option<u8>, D::Error>
where
    D: serde::Deserializer<'de>, {
    let s = String::deserialize(deserializer)?;
    if s == "." {
        Ok(None)
    }
    else {
        s.parse::<u8>()
            .map(Some)
            .map_err(|e| serde::de::Error::custom(format!("Failed to parse u8: {}", e)))
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RawGffEntry {
    pub seqid:        BsxSmallStr,
    pub source:       BsxSmallStr,
    pub feature_type: BsxSmallStr,
    pub start:        u32,
    pub end:          u32,
    #[serde(deserialize_with = "deserialize_optional_f64")]
    pub score:        Option<f64>,
    pub strand:       char,
    #[serde(deserialize_with = "deserialize_optional_u8")]
    pub phase:        Option<u8>,
    pub attributes:   String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GffEntry {
    pub contig:       Contig,
    pub source:       ArcStr,
    pub feature_type: ArcStr,
    pub score:        Option<f64>,
    pub phase:        Option<u8>,
    pub attributes:   GffEntryAttributes,
    pub id:           BsxSmallStr,
}

impl From<bio::io::bed::Record> for GffEntry {
    fn from(value: bio::io::bed::Record) -> Self {
        let contig = Contig::from(value.clone());
        Self::new(
            contig,
            None,
            None,
            value.score().map(|s| s.parse::<f64>().unwrap_or(f64::NAN)),
            None,
            None,
        )
    }
}

impl GffEntry {
    getter_fn!(contig, Contig);

    getter_fn!(source, ArcStr);

    getter_fn!(feature_type, ArcStr);

    getter_fn!(*score, Option<f64>);

    getter_fn!(*phase, Option<u8>);

    getter_fn!(attributes, GffEntryAttributes);

    getter_fn!(id, BsxSmallStr);

    pub fn new(
        contig: Contig,
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
            .unwrap_or(nanoid!(16).into());
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
        let seqid = value.seqid;
        let source = ArcStr::from(value.source.as_str());
        let strand = Strand::from_str(value.strand.to_string().as_str()).unwrap();
        let feature_type = ArcStr::from(value.feature_type.as_str());
        let attributes = GffEntryAttributes::from_str(value.attributes.as_str())?;

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
            seqid: value.contig.seqname().as_str().into(),
            source: value.source.as_str().into(),
            feature_type: value.feature_type.as_str().into(),
            start: value.contig.start(),
            end: value.contig.end(),
            score: value.score,
            strand: value.contig.strand().into(),
            phase: value.phase,
            attributes,
        })
    }
}
