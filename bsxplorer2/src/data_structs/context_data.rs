use itertools::Itertools;
use polars::df;
use polars::frame::DataFrame;
#[allow(unused_imports)]
use polars::prelude::*;

use crate::data_structs::batch::{BatchType, BsxTypeTag};
use crate::data_structs::enums::{Context, IPCEncodedEnum, Strand}; // Added for testing DataFrame creation

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
    /// Creates a new `ContextData` from vectors of positions, strands, and
    /// contexts.
    pub fn new(
        positions: Vec<u32>,
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
        start: u32,
    ) {
        for (shift, trinuc) in seq.windows(3).enumerate() {
            match trinuc {
                [c1, c2, _]
                    if (c1.eq_ignore_ascii_case(&b'C')
                        && c2.eq_ignore_ascii_case(&b'G')) =>
                {
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
                [c1, c2, c3]
                    if (c1.eq_ignore_ascii_case(&b'C')
                        && c3.eq_ignore_ascii_case(&b'G')
                        && !c2.eq_ignore_ascii_case(&b'C')) =>
                {
                    self.entries.insert(Entry(
                        start + shift as u32 + 2,
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
                        start + shift as u32,
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
                        start + shift as u32,
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
                        start + shift as u32 + 2,
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
                start + seq.len() as u32 - 2,
                Strand::Forward,
                Context::CG,
            ));
            self.entries.insert(Entry(
                start + seq.len() as u32 - 1,
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

    /// Converts the `ContextData` to a Polars `DataFrame`.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_structs::batch::{BatchType, BsxTypeTag};
    use crate::data_structs::enums::{Context, Strand};

    // Mock structs implementing BsxTypeTag for testing to_df
    struct MockDecodedBsxBatch;
    impl BsxTypeTag for MockDecodedBsxBatch {
        fn type_name() -> &'static str {
            "mock_decoded"
        }

        fn type_enum() -> BatchType {
            BatchType::Decoded
        }
    }

    struct MockEncodedBsxBatch;
    impl BsxTypeTag for MockEncodedBsxBatch {
        fn type_name() -> &'static str {
            "mock_encoded"
        }

        fn type_enum() -> BatchType {
            BatchType::Encoded
        }
    }

    #[test]
    fn test_new() {
        let positions = vec![10, 20, 30];
        let strands = vec![Strand::Forward, Strand::Reverse, Strand::None];
        let contexts = vec![Context::CG, Context::CHG, Context::CHH];
        let data = ContextData::new(
            positions.clone(),
            strands.clone(),
            contexts.clone(),
        );

        assert_eq!(data.len(), 3);
        let (p, s, c) = data.take(); // Order should be preserved by BTreeSet
        assert_eq!(p, vec![10, 20, 30]);
        assert_eq!(s, vec![Strand::Forward, Strand::Reverse, Strand::None]);
        assert_eq!(c, vec![Context::CG, Context::CHG, Context::CHH]);
    }

    #[test]
    fn test_from_sequence_basic() {
        // Sequence: "CCGTAG", start=1
        // Tracing based on the code's specific logic:
        // Shift 0 (pos 1): CCG
        //   [C,_,G] [1]!=G (C!=G): Yes. Entry(1, F, CHG)
        // Shift 1 (pos 2): CGT
        //   [C,G,_]: Yes. Entry(2, F, CG), Entry(3, R, CG)
        //   [C,_,_] [1]!=G && [2]!=G (G!=G && T!=G): No.
        // Shift 2 (pos 3): GTA - No C.. or ..G match.
        // Shift 3 (pos 4): TAG
        //   [_,_,G] [0]!=C && [1]!=C (T!=C && A!=C): Yes. Entry(1+3+2, R, CHH)
        // -> (6, R, CHH). Suffix "AG": No CG.
        // Unique entries: (1, F, CHG), (2, F, CG), (3, R, CG), (6, R, CHH)
        // Sorted: (1, F, CHG), (2, F, CG), (3, R, CG), (6, R, CHH)
        // Pos: 1, 2, 3, 6
        // Strands: F, F, R, R
        // Contexts: CHG, CG, CG, CHH

        let seq = b"CCGTAG";
        let data = ContextData::from_sequence(seq);
        let (p, s, c) = data.take();

        assert_eq!(p, vec![1, 2, 3, 6]);
        assert_eq!(s, vec![
            Strand::Forward,
            Strand::Forward,
            Strand::Reverse,
            Strand::Reverse
        ]);
        assert_eq!(c, vec![
            Context::CHG,
            Context::CG,
            Context::CG,
            Context::CHH
        ]);
    }

    #[test]
    fn test_read_sequence_start_position() {
        let seq = b"CCGTAG"; // Expect (1, F, CHG), (2, F, CG), (3, R, CG), (6, R, CHH) when start=1
        let mut data = ContextData::empty();
        data.read_sequence(seq, 10); // Expect positions 10 + (original_pos - 1)

        // Original positions: 1, 2, 3, 6
        // New positions (start=10): 10+(1-1)=10, 10+(2-1)=11, 10+(3-1)=12,
        // 10+(6-1)=15 Expected entries: (10, F, CHG), (11, F, CG), (12,
        // R, CG), (15, R, CHH) Sorted: (10, F, CHG), (11, F, CG), (12,
        // R, CG), (15, R, CHH)

        let (p, s, c) = data.take();
        assert_eq!(p, vec![10, 11, 12, 15]);
        assert_eq!(s, vec![
            Strand::Forward,
            Strand::Forward,
            Strand::Reverse,
            Strand::Reverse
        ]);
        assert_eq!(c, vec![
            Context::CHG,
            Context::CG,
            Context::CG,
            Context::CHH
        ]);
    }

    #[test]
    fn test_read_sequence_suffix_cg() {
        // Test the suffix rule handling CG at the end
        // Sequence: "ATCG", start = 1
        // Shift 0 (pos 1): ATC
        //   [C,_,_] [1]!=G && [2]!=G (A!=G && T!=G): Yes. Entry(1, F, CHH).
        // Shift 1 (pos 2): TCG
        //   [C,G,_]: Yes. Entry(1+1, F, CG) -> (2, F, CG). Entry(1+1+1, R, CG)
        // -> (3, R, CG). Suffix "CG" at end (len=4, index 3): Yes. Pos
        // 1 + 4 - 1 = 4. Entry(4, R, CG). Unique entries: (1, F, CHH),
        // (2, F, CG), (3, R, CG), (4, R, CG) Sorted: (1, F, CHH), (2,
        // F, CG), (3, R, CG), (4, R, CG) Pos: 1, 2, 3, 4
        // Strands: F, F, R, R
        // Contexts: CHH, CG, CG, CG

        let mut data = ContextData::empty();
        data.read_sequence(b"ATCG", 1);
        let (p, s, c) = data.take();

        assert_eq!(p, vec![3, 4]);
        assert_eq!(s, vec![Strand::Forward, Strand::Reverse]);
        assert_eq!(c, vec![Context::CG, Context::CG]);
    }

    #[test]
    fn test_read_sequence_empty_or_short() {
        let mut data = ContextData::empty();
        data.read_sequence(b"", 1);
        assert!(data.is_empty());

        let mut data = ContextData::empty();
        data.read_sequence(b"A", 1); // Length < 3
        assert!(data.is_empty());

        let mut data = ContextData::empty();
        data.read_sequence(b"AT", 1); // Length < 3
        assert!(data.is_empty());

        let mut data = ContextData::empty();
        data.read_sequence(b"ATC", 1); // No trinuc matches, no CG suffix
        assert!(data.is_empty());
    }

    #[test]
    fn test_drain_until() {
        let mut data = ContextData::new(
            vec![10, 20, 30, 40, 50],
            vec![
                Strand::Forward,
                Strand::Reverse,
                Strand::Forward,
                Strand::Reverse,
                Strand::Forward,
            ],
            vec![
                Context::CG,
                Context::CHG,
                Context::CHH,
                Context::CG,
                Context::CHG,
            ],
        );

        // Drain up to pos 30 (exclusive)
        let drained = data.drain_until(20);
        assert_eq!(drained.len(), 2);
        let (p, s, c) = drained.take();
        assert_eq!(p, vec![10, 20]); // Entries with pos < 30
        assert_eq!(s, vec![Strand::Forward, Strand::Reverse]);
        assert_eq!(c, vec![Context::CG, Context::CHG]);

        // Check remaining data
        assert_eq!(data.len(), 3);
        let (p, s, c) = data.take();
        assert_eq!(p, vec![30, 40, 50]); // Entries with pos >= 30
        assert_eq!(s, vec![Strand::Forward, Strand::Reverse, Strand::Forward]);
        assert_eq!(c, vec![Context::CHH, Context::CG, Context::CHG]);

        // Drain again, draining all remaining
        let mut data = ContextData::new(
            vec![30, 40, 50],
            vec![Strand::Forward, Strand::Reverse, Strand::Forward],
            vec![Context::CHH, Context::CG, Context::CHG],
        );
        let drained_all = data.drain_until(100);
        assert_eq!(drained_all.len(), 3);
        let (p, s, c) = drained_all.take();
        assert_eq!(p, vec![30, 40, 50]);
        assert_eq!(s, vec![Strand::Forward, Strand::Reverse, Strand::Forward]);
        assert_eq!(c, vec![Context::CHH, Context::CG, Context::CHG]);

        // Check remaining data (should be empty)
        assert_eq!(data.len(), 0);
        assert!(data.is_empty());

        // Drain from empty data
        let mut empty_data = ContextData::empty();
        let drained_empty = empty_data.drain_until(100);
        assert_eq!(drained_empty.len(), 0);
        assert!(drained_empty.is_empty());
        assert_eq!(empty_data.len(), 0);
        assert!(empty_data.is_empty());

        // Drain with pos smaller than any entry
        let mut data = ContextData::new(
            vec![10, 20],
            vec![Strand::Forward, Strand::Reverse],
            vec![Context::CG, Context::CHG],
        );
        let drained_none = data.drain_until(5);
        assert_eq!(drained_none.len(), 0);
        assert!(drained_none.is_empty());
        assert_eq!(data.len(), 2); // Original data remains
    }

    #[test]
    fn test_len_and_is_empty() {
        let empty_data = ContextData::empty();
        assert_eq!(empty_data.len(), 0);
        assert!(empty_data.is_empty());

        let data = ContextData::new(
            vec![10, 20, 30],
            vec![Strand::Forward, Strand::Reverse, Strand::None],
            vec![Context::CG, Context::CHG, Context::CHH],
        );
        assert_eq!(data.len(), 3);
        assert!(!data.is_empty());
    }

    #[test]
    fn test_take() {
        let positions = vec![10, 20, 30];
        let strands = vec![Strand::Forward, Strand::Reverse, Strand::None];
        let contexts = vec![Context::CG, Context::CHG, Context::CHH];
        let data = ContextData::new(
            positions.clone(),
            strands.clone(),
            contexts.clone(),
        );

        let (p, s, c) = data.take();
        // BTreeSet orders by Entry Ord: pos -> strand -> context
        // (10, F, CG), (20, R, CHG), (30, None, CHH)
        assert_eq!(p, vec![10, 20, 30]);
        assert_eq!(s, vec![Strand::Forward, Strand::Reverse, Strand::None]);
        assert_eq!(c, vec![Context::CG, Context::CHG, Context::CHH]);

        // Ensure the original struct is effectively empty after take (data is
        // moved out) Cloning before take is needed if we want to check
        // the original struct afterwards
        let data_to_take = ContextData::new(positions, strands, contexts);
        let _ = data_to_take.take(); // Move data out
    }

    #[test]
    fn test_to_df_decoded() {
        let positions = vec![10, 20];
        let strands = vec![Strand::Forward, Strand::Reverse];
        let contexts = vec![Context::CG, Context::CHG];
        let data = ContextData::new(positions, strands, contexts);

        let df = data.to_df::<MockDecodedBsxBatch>();

        let expected_df = df!(
            "position" => &[10u64, 20], // Decoded uses UInt64 for position
            "strand" => &["+", "-"],
            "context" => &["CG", "CHG"],
        )
        .unwrap();

        assert_eq!(df, expected_df);
    }

    #[test]
    fn test_to_df_encoded() {
        let positions = vec![10, 20];
        let strands = vec![Strand::Forward, Strand::Reverse];
        let contexts = vec![Context::CG, Context::CHG];
        let data = ContextData::new(positions, strands, contexts);

        let df = data.to_df::<MockEncodedBsxBatch>();

        // Encoded uses bool for strand/context and UInt32 for position
        let expected_df = df!(
            "position" => &[10u32, 20],
            "strand" => &[Some(true), Some(false)], // Forward -> Some(true), Reverse -> Some(false)
            "context" => &[Some(true), Some(false)], // CG -> Some(true), CHG -> Some(false)
        )
        .unwrap();

        assert_eq!(df, expected_df);
    }

    #[test]
    fn test_to_df_empty() {
        let data = ContextData::empty();

        let df_decoded = data
            .clone()
            .to_df::<MockDecodedBsxBatch>();
        let expected_df_decoded = df!(
            "position" => Series::new_empty("position".into(), &DataType::UInt64),
            "strand" => Series::new_empty("strand".into(), &DataType::String),
            "context" => Series::new_empty("context".into(), &DataType::String),
        ).unwrap();
        assert_eq!(df_decoded.schema(), expected_df_decoded.schema());
        assert_eq!(df_decoded.height(), 0);

        let df_encoded = data.to_df::<MockEncodedBsxBatch>();
        let expected_df_encoded = df!(
            "position" => Series::new_empty("position".into(), &DataType::UInt32),
            "strand" => Series::new_empty("strand".into(), &DataType::Boolean),
            "context" => Series::new_empty("context".into(), &DataType::Boolean),
        ).unwrap();
        assert_eq!(df_encoded.schema(), expected_df_encoded.schema());
        assert_eq!(df_encoded.height(), 0);
    }
}
