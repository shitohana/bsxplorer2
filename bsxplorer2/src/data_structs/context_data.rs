use itertools::Itertools;
use polars::df;
use polars::frame::DataFrame;
#[allow(unused_imports)]
use polars::prelude::*;

use super::typedef::PosType;
use crate::data_structs::enums::{
    Context,
    Strand,
};

#[derive(Debug, Clone)]
pub struct Entry(PosType, Strand, Context);

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
    /// Creates a new `ContextData` from vectors of positions, strands, and
    /// contexts.
    pub fn new(
        positions: Vec<PosType>,
        strands: Vec<Strand>,
        contexts: Vec<Context>,
    ) -> Self {
        assert_eq!(positions.len(), strands.len());
        assert_eq!(positions.len(), contexts.len());

        let mut entries = std::collections::BTreeSet::new();
        #[allow(unsafe_code)]
        unsafe {
            for i in 0..positions.len() {
                entries.insert(Entry(
                    *positions.get_unchecked(i),
                    *strands.get_unchecked(i),
                    *contexts.get_unchecked(i),
                ));
            }
        }


        Self { entries }
    }

    /// Creates a `ContextData` from a sequence of bytes.
    pub fn from_sequence(seq: &[u8]) -> Self {
        let mut new_self = Self {
            entries: std::collections::BTreeSet::new(),
        };
        new_self.read_sequence(seq, 1);
        new_self
    }

    /// Creates an empty `ContextData`.
    pub fn empty() -> Self {
        Self {
            entries: std::collections::BTreeSet::new(),
        }
    }

    /// Reads a sequence and populates the `ContextData` with entries.
    pub fn read_sequence(
        &mut self,
        seq: &[u8],
        start: PosType,
    ) {
        for (shift, trinuc) in seq.windows(3).enumerate() {
            match trinuc {
                [c1, c2, _]
                    if (c1.eq_ignore_ascii_case(&b'C')
                        && c2.eq_ignore_ascii_case(&b'G')) =>
                {
                    self.entries.insert(Entry(
                        start + shift as PosType,
                        Strand::Forward,
                        Context::CG,
                    ));
                    self.entries.insert(Entry(
                        start + shift as PosType + 1,
                        Strand::Reverse,
                        Context::CG,
                    ));
                },
                [c1, c2, c3]
                    if (c1.eq_ignore_ascii_case(&b'C')
                        && c3.eq_ignore_ascii_case(&b'G')
                        && !c2.eq_ignore_ascii_case(&b'C')) =>
                {
                    self.entries.insert(Entry(
                        start + shift as PosType + 2,
                        Strand::Reverse,
                        Context::CHG,
                    ));
                },
                [c1, c2, c3]
                    if (c1.eq_ignore_ascii_case(&b'C')
                        && c3.eq_ignore_ascii_case(&b'G')
                        && !c2.eq_ignore_ascii_case(&b'G')) =>
                {
                    self.entries.insert(Entry(
                        start + shift as PosType,
                        Strand::Forward,
                        Context::CHG,
                    ));
                },
                [c1, c2, c3]
                    if (c1.eq_ignore_ascii_case(&b'C')
                        && !c2.eq_ignore_ascii_case(&b'G')
                        && !c3.eq_ignore_ascii_case(&b'G')) =>
                {
                    self.entries.insert(Entry(
                        start + shift as PosType,
                        Strand::Forward,
                        Context::CHH,
                    ));
                },
                [c1, c2, c3]
                    if (c3.eq_ignore_ascii_case(&b'G')
                        && !c1.eq_ignore_ascii_case(&b'C')
                        && !c2.eq_ignore_ascii_case(&b'C')) =>
                {
                    self.entries.insert(Entry(
                        start + shift as PosType + 2,
                        Strand::Reverse,
                        Context::CHH,
                    ));
                },
                _ => continue,
            }
        }
        if seq.len() >= 2
            && seq
                .get(seq.len().saturating_sub(2))
                .map(|v| v.eq_ignore_ascii_case(&b'C'))
                .unwrap_or(false)
            && seq
                .get(seq.len().saturating_sub(1))
                .map(|v| v.eq_ignore_ascii_case(&b'G'))
                .unwrap_or(false)
        {
            self.entries.insert(Entry(
                start + seq.len() as PosType - 2,
                Strand::Forward,
                Context::CG,
            ));
            self.entries.insert(Entry(
                start + seq.len() as PosType - 1,
                Strand::Reverse,
                Context::CG,
            ));
        }
    }

    /// Drains entries until a given position, returning the drained entries.
    pub fn drain_until(
        &mut self,
        pos: u32,
    ) -> Self {
        let split_point = self.entries.iter().find(|entry| entry.0 > pos).cloned();

        let mut drained = std::collections::BTreeSet::new();
        if let Some(split_entry) = split_point {
            let remaining = self.entries.split_off(&split_entry);
            std::mem::swap(&mut drained, &mut self.entries);
            self.entries = remaining;
        }
        else {
            std::mem::swap(&mut drained, &mut self.entries);
        }

        Self { entries: drained }
    }

    /// Returns the number of entries in the `ContextData`.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Checks if the `ContextData` is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Takes the data, returning vectors of positions, strands, and contexts.
    pub fn take(self) -> (Vec<PosType>, Vec<Strand>, Vec<Context>) {
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

    /// Converts the `ContextData` to a Polars `DataFrame`.
    pub fn to_df(self) -> DataFrame {
        let (positions, strands, contexts) = self.take();

        df!(
            "position" => positions,
            "strand" => strands.into_iter().map(Option::<bool>::from).collect_vec(),
            "context" => contexts.into_iter().map(Option::<bool>::from).collect_vec(),
        )
        .unwrap()
    }
}
