use crate::data_structs::batch::{BatchType, BsxTypeTag};
use crate::utils::types::{Context, IPCEncodedEnum, Strand};
use itertools::Itertools;
use polars::df;
use polars::frame::DataFrame;

#[derive(Debug, Clone)]
pub struct Entry(u32, Strand, Context);

impl Ord for Entry {
    fn cmp(
        &self,
        other: &Self,
    ) -> std::cmp::Ordering {
        self.0
            .cmp(&other.0)
            .then_with(|| self.1.cmp(&other.1))
            .then_with(|| self.2.cmp(&other.2))
    }
}

impl PartialOrd for Entry {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Entry {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.0 == other.0 && self.1 == other.1 && self.2 == other.2
    }
}

impl Eq for Entry {}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct ContextData {
    entries: std::collections::BTreeSet<Entry>,
}

impl ContextData {
    pub fn new(
        positions: Vec<u32>,
        strands: Vec<Strand>,
        contexts: Vec<Context>,
    ) -> Self {
        assert_eq!(positions.len(), strands.len());
        assert_eq!(positions.len(), contexts.len());

        let mut entries = std::collections::BTreeSet::new();
        for i in 0..positions.len() {
            entries.insert(Entry(positions[i], strands[i], contexts[i]));
        }

        Self { entries }
    }

    pub fn from_sequence(seq: Vec<u8>) -> Self {
        let mut new_self = Self {
            entries: std::collections::BTreeSet::new(),
        };
        new_self.read_sequence(&seq, 1);
        new_self
    }

    pub fn empty() -> Self {
        Self {
            entries: std::collections::BTreeSet::new(),
        }
    }

    pub fn read_sequence(
        &mut self,
        seq: &[u8],
        start: u32,
    ) {
        let uppercase = seq.to_ascii_uppercase();

        for (shift, trinuc) in uppercase.windows(3).enumerate() {
            match trinuc {
                &[b'C', b'G', _] => {
                    self.entries.insert(Entry(
                        start + shift as u32,
                        Strand::Forward,
                        Context::CG,
                    ));
                    self.entries.insert(Entry(
                        start + shift as u32 + 1,
                        Strand::Reverse,
                        Context::CG,
                    ));
                },
                &[b'C', _, b'G'] => {
                    self.entries.insert(Entry(
                        start + shift as u32,
                        Strand::Forward,
                        Context::CHG,
                    ));
                    self.entries.insert(Entry(
                        start + shift as u32 + 2,
                        Strand::Reverse,
                        Context::CHG,
                    ));
                },
                &[b'C', _, _] => {
                    self.entries.insert(Entry(
                        start + shift as u32,
                        Strand::Forward,
                        Context::CHH,
                    ));
                },
                &[_, _, b'G'] => {
                    self.entries.insert(Entry(
                        start + shift as u32 + 2,
                        Strand::Reverse,
                        Context::CHH,
                    ));
                },
                _ => continue,
            }
        }
        if let Some(&[b'C', b'G']) =
            uppercase.get(uppercase.len().saturating_sub(2)..)
        {
            self.entries.insert(Entry(
                start + uppercase.len() as u32 - 1,
                Strand::Reverse,
                Context::CG,
            ));
        }
    }

    pub fn drain_until(
        &mut self,
        pos: u32,
    ) -> Self {
        let split_point = self
            .entries
            .iter()
            .find(|entry| entry.0 > pos)
            .cloned();

        let mut drained = std::collections::BTreeSet::new();
        if let Some(split_entry) = split_point {
            let remaining = self.entries.split_off(&split_entry);
            std::mem::swap(&mut drained, &mut self.entries);
            self.entries = remaining;
        } else {
            std::mem::swap(&mut drained, &mut self.entries);
        }

        Self { entries: drained }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn take(self) -> (Vec<u32>, Vec<Strand>, Vec<Context>) {
        let mut positions = Vec::with_capacity(self.entries.len());
        let mut strands = Vec::with_capacity(self.entries.len());
        let mut contexts = Vec::with_capacity(self.entries.len());

        for entry in self.entries {
            positions.push(entry.0);
            strands.push(entry.1);
            contexts.push(entry.2);
        }

        (positions, strands, contexts)
    }

    pub fn to_df<B: BsxTypeTag>(self) -> DataFrame {
        let (positions, strands, contexts) = self.take();

        match B::type_enum() {
            BatchType::Decoded => {
                df!(
                    "position" => positions.into_iter().map(|x| x as u64).collect_vec(),
                    "strand" => strands.into_iter().map(|x| x.to_string()).collect_vec(),
                    "context" => contexts.into_iter().map(|x| x.to_string()).collect_vec(),
                ).unwrap()
            }
            BatchType::Encoded => {
                df!(
                    "position" => positions,
                    "strand" => strands.into_iter().map(|x| x.to_bool()).collect_vec(),
                    "context" => contexts.into_iter().map(|x| x.to_bool()).collect_vec(),
                ).unwrap()
            }
        }
    }
}
