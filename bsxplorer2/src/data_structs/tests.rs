mod enums_tests {
    use std::str::FromStr;

    use crate::data_structs::enums::*;

    // --- Context Tests ---

    #[test]
    fn test_context_from_str() {
        assert_eq!(Context::from_str("CG").unwrap(), Context::CG);
        assert_eq!(Context::from_str("cg").unwrap(), Context::CG);
        assert_eq!(Context::from_str("CHG").unwrap(), Context::CHG);
        assert_eq!(Context::from_str("chg").unwrap(), Context::CHG);
        assert_eq!(Context::from_str("CHH").unwrap(), Context::CHH);
        assert_eq!(Context::from_str("chh").unwrap(), Context::CHH);
    }

    #[test]
    fn test_context_ipc_encoded() {
        assert_eq!(Context::from(Some(true)), Context::CG);
        assert_eq!(Context::from(Some(false)), Context::CHG);
        assert_eq!(Context::from(None), Context::CHH);

        assert_eq!(Option::<bool>::from(Context::CG), Some(true));
        assert_eq!(Option::<bool>::from(Context::CHG), Some(false));
        assert_eq!(Option::<bool>::from(Context::CHH), None);
    }

    // --- Strand Tests ---

    #[test]
    fn test_strand_from_str() {
        assert_eq!(Strand::from_str("+").unwrap(), Strand::Forward);
        assert_eq!(Strand::from_str("-").unwrap(), Strand::Reverse);
        assert_eq!(Strand::from_str(".").unwrap(), Strand::None);
        assert_eq!(Strand::from_str("forward").unwrap(), Strand::None); // Defaults to None
        assert_eq!(Strand::from_str("").unwrap(), Strand::None);
        assert_eq!(Strand::from_str("AnythingElse").unwrap(), Strand::None);
    }

    #[test]
    fn test_strand_ipc_encoded() {
        assert_eq!(Strand::from(Some(true)), Strand::Forward);
        assert_eq!(Strand::from(Some(false)), Strand::Reverse);
        assert_eq!(Strand::from(None), Strand::None);

        assert_eq!(Option::<bool>::from(Strand::Forward), Some(true));
        assert_eq!(Option::<bool>::from(Strand::Reverse), Some(false));
        assert_eq!(Option::<bool>::from(Strand::None), None);
    }
}

mod context_data_tests {
    use polars::prelude::*;
    use rstest::rstest;

    use crate::data_structs::context_data::ContextData;
    use crate::data_structs::enums::{
        Context,
        Strand,
    };

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

    #[rstest]
    #[case(b"", 1)]
    #[case(b"A", 1)] // Length < 3
    #[case(b"AT", 1)] // Length < 3
    #[case(b"ATC", 1)] // No trinuc matches, no CG suffix
    fn test_read_sequence_empty_or_short(
        #[case] seq: &[u8],
        #[case] start: u32,
    ) {
        let mut data = ContextData::empty();
        data.read_sequence(seq, start);
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
    fn test_to_df() {
        let positions = vec![10, 20];
        let strands = vec![Strand::Forward, Strand::Reverse];
        let contexts = vec![Context::CG, Context::CHG];
        let data = ContextData::new(positions, strands, contexts);

        let df = data.to_df();

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

        let df_decoded = data.clone().to_df();
        let expected_df_decoded = df!(
            "position" => Series::new_empty("position".into(), &DataType::UInt32),
            "strand" => Series::new_empty("strand".into(), &DataType::Boolean),
            "context" => Series::new_empty("context".into(), &DataType::Boolean),
        )
        .unwrap();
        assert_eq!(df_decoded.schema(), expected_df_decoded.schema());
        assert_eq!(df_decoded.height(), 0);
    }
}
