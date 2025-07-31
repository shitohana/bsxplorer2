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
        // Test fallback to CHH for any other string
        assert_eq!(Context::from_str("unknown").unwrap(), Context::CHH);
        assert_eq!(Context::from_str("").unwrap(), Context::CHH);
        assert_eq!(Context::from_str("xyz").unwrap(), Context::CHH);
    }

    #[test]
    fn test_context_display() {
        assert_eq!(Context::CG.to_string(), "CG");
        assert_eq!(Context::CHG.to_string(), "CHG");
        assert_eq!(Context::CHH.to_string(), "CHH");
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

    #[test]
    fn test_context_serialization() {
        // Test JSON serialization/deserialization
        let cg = Context::CG;
        let serialized = serde_json::to_string(&cg).unwrap();
        assert_eq!(serialized, "\"CG\"");
        let deserialized: Context = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, Context::CG);

        let chg = Context::CHG;
        let serialized = serde_json::to_string(&chg).unwrap();
        assert_eq!(serialized, "\"CHG\"");
        let deserialized: Context = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, Context::CHG);

        let chh = Context::CHH;
        let serialized = serde_json::to_string(&chh).unwrap();
        assert_eq!(serialized, "\"CHH\"");
        let deserialized: Context = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, Context::CHH);
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
    fn test_strand_display() {
        assert_eq!(Strand::Forward.to_string(), "+");
        assert_eq!(Strand::Reverse.to_string(), "-");
        assert_eq!(Strand::None.to_string(), ".");
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

    #[test]
    fn test_strand_from_bool() {
        assert_eq!(Strand::from(true), Strand::Forward);
        assert_eq!(Strand::from(false), Strand::Reverse);

        assert_eq!(bool::from(Strand::Forward), true);
        assert_eq!(bool::from(Strand::Reverse), false);
        // Note: bool::from(Strand::None) would panic with unimplemented!()
    }

    #[test]
    fn test_strand_from_char() {
        assert_eq!(char::from(Strand::Forward), '+');
        assert_eq!(char::from(Strand::Reverse), '-');
        assert_eq!(char::from(Strand::None), '.');
    }

    #[test]
    #[should_panic]
    fn test_strand_none_to_bool_panics() {
        let _: bool = Strand::None.into();
    }

    #[test]
    fn test_strand_serialization() {
        // Test JSON serialization/deserialization
        let forward = Strand::Forward;
        let serialized = serde_json::to_string(&forward).unwrap();
        assert_eq!(serialized, "\"+\"");
        let deserialized: Strand = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, Strand::Forward);

        let reverse = Strand::Reverse;
        let serialized = serde_json::to_string(&reverse).unwrap();
        assert_eq!(serialized, "\"-\"");
        let deserialized: Strand = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, Strand::Reverse);

        let none = Strand::None;
        let serialized = serde_json::to_string(&none).unwrap();
        assert_eq!(serialized, "\".\"");
        let deserialized: Strand = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, Strand::None);
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


#[cfg(test)]
mod methstats_tests {
    use hashbrown::HashMap;
    use assert_approx_eq::assert_approx_eq;
    use rstest::rstest;

    use crate::data_structs::enums::{
        Context,
        Strand,
    };
    use crate::data_structs::methstats::{
        MethAgg,
        RegionMethAgg,
    };
    use crate::data_structs::typedef::DensityType;

    // --- MethAgg Tests ---

    #[test]
    fn test_meth_agg_new() {
        let agg = MethAgg::new();
        assert_approx_eq!(agg.sum(), 0.0);
        assert_approx_eq!(agg.count(), 0.0);
    }

    #[test]
    fn test_meth_agg_default() {
        let agg = MethAgg::default();
        assert_approx_eq!(agg.sum(), 0.0);
        assert_approx_eq!(agg.count(), 0.0);
    }

    #[test]
    fn test_meth_agg_from_tuple() {
        let agg = MethAgg::from((5.0, 10.0));
        assert_approx_eq!(agg.sum(), 5.0);
        assert_approx_eq!(agg.count(), 10.0);
    }

    #[test]
    fn test_meth_agg_add_density() {
        let mut agg = MethAgg::new();
        agg.add_density(0.8);
        assert_approx_eq!(agg.sum(), 0.8);
        assert_approx_eq!(agg.count(), 1.0);

        agg.add_density(0.6);
        assert_approx_eq!(agg.sum(), 1.4);
        assert_approx_eq!(agg.count(), 2.0);
    }

    #[test]
    fn test_meth_agg_finalize() {
        let mut agg = MethAgg::new();
        agg.add_density(0.8);
        agg.add_density(0.6);
        assert_approx_eq!(agg.finalize(), 0.7); // (0.8 + 0.6) / 2
    }

    #[test]
    fn test_meth_agg_add_assign() {
        let mut agg1 = MethAgg::from((2.0, 3.0));
        let agg2 = MethAgg::from((4.0, 5.0));
        agg1 += agg2;
        assert_approx_eq!(agg1.sum(), 6.0);
        assert_approx_eq!(agg1.count(), 8.0);
    }

    #[test]
    fn test_meth_agg_add() {
        let agg1 = MethAgg::from((2.0, 3.0));
        let agg2 = MethAgg::from((4.0, 5.0));
        let result = agg1 + agg2;
        assert_approx_eq!(result.sum(), 6.0);
        assert_approx_eq!(result.count(), 8.0);
    }

    // --- RegionMethAgg Tests ---

    #[test]
    fn test_region_meth_agg_new() {
        let agg = RegionMethAgg::new();
        assert!(agg.is_empty());
        assert_eq!(agg.context().len(), 0);
        assert_eq!(agg.strand().len(), 0);
        assert_eq!(agg.coverage().len(), 0);
    }

    #[test]
    fn test_region_meth_agg_default() {
        let agg = RegionMethAgg::default();
        assert!(agg.is_empty());
    }

    #[test]
    fn test_region_meth_agg_full() {
        let mut coverage = HashMap::new();
        coverage.insert(10, 5);
        let mut context = HashMap::new();
        context.insert(Context::CG, MethAgg::from((8.0, 10.0)));
        let mut strand = HashMap::new();
        strand.insert(Strand::Forward, MethAgg::from((6.0, 8.0)));

        let agg = RegionMethAgg::full(coverage.clone(), context.clone(), strand.clone());
        assert!(!agg.is_empty());
        assert_eq!(agg.coverage(), &coverage);
        assert_eq!(agg.context(), &context);
        assert_eq!(agg.strand(), &strand);
    }

    #[test]
    fn test_region_meth_agg_add_cytosine() {
        let mut agg = RegionMethAgg::new();

        // Add cytosine with coverage
        agg.add_cytosine(8.0, Some(10), Context::CG, Strand::Forward);

        assert!(!agg.is_empty());
        assert_eq!(agg.coverage().get(&10), Some(&1));

        let context_agg = agg.context().get(&Context::CG).unwrap();
        assert_approx_eq!(context_agg.sum(), 0.8); // 8.0 / 10
        assert_approx_eq!(context_agg.count(), 1.0);

        let strand_agg = agg.strand().get(&Strand::Forward).unwrap();
        assert_approx_eq!(strand_agg.sum(), 0.8);
        assert_approx_eq!(strand_agg.count(), 1.0);
    }

    #[test]
    fn test_region_meth_agg_add_cytosine_no_coverage() {
        let mut agg = RegionMethAgg::new();

        // Add cytosine without coverage
        agg.add_cytosine(5.0, None, Context::CHG, Strand::Reverse);

        assert!(!agg.is_empty());
        assert_eq!(agg.coverage().len(), 0); // No coverage recorded

        let context_agg = agg.context().get(&Context::CHG).unwrap();
        assert_approx_eq!(context_agg.sum(), 5.0); // 5.0 / 1 (default when None)
        assert_approx_eq!(context_agg.count(), 1.0);

        let strand_agg = agg.strand().get(&Strand::Reverse).unwrap();
        assert_approx_eq!(strand_agg.sum(), 5.0);
        assert_approx_eq!(strand_agg.count(), 1.0);
    }

    #[test]
    fn test_region_meth_agg_multiple_cytosines() {
        let mut agg = RegionMethAgg::new();

        agg.add_cytosine(8.0, Some(10), Context::CG, Strand::Forward);
        agg.add_cytosine(6.0, Some(10), Context::CG, Strand::Forward);
        agg.add_cytosine(4.0, Some(5), Context::CHG, Strand::Reverse);

        // Check coverage counts
        assert_eq!(agg.coverage().get(&10), Some(&2));
        assert_eq!(agg.coverage().get(&5), Some(&1));

        // Check context aggregation
        let cg_agg = agg.context().get(&Context::CG).unwrap();
        assert_approx_eq!(cg_agg.sum(), 1.4); // (8.0/10) + (6.0/10) = 0.8 + 0.6
        assert_approx_eq!(cg_agg.count(), 2.0);

        let chg_agg = agg.context().get(&Context::CHG).unwrap();
        assert_approx_eq!(chg_agg.sum(), 0.8); // 4.0/5
        assert_approx_eq!(chg_agg.count(), 1.0);

        // Check strand aggregation
        let forward_agg = agg.strand().get(&Strand::Forward).unwrap();
        assert_approx_eq!(forward_agg.sum(), 1.4); // 0.8 + 0.6
        assert_approx_eq!(forward_agg.count(), 2.0);

        let reverse_agg = agg.strand().get(&Strand::Reverse).unwrap();
        assert_approx_eq!(reverse_agg.sum(), 0.8);
        assert_approx_eq!(reverse_agg.count(), 1.0);
    }

    #[test]
    fn test_region_meth_agg_finalize_context() {
        let mut agg = RegionMethAgg::new();
        agg.add_cytosine(8.0, Some(10), Context::CG, Strand::Forward);
        agg.add_cytosine(6.0, Some(10), Context::CG, Strand::Forward);
        agg.add_cytosine(4.0, Some(5), Context::CHG, Strand::Reverse);

        let finalized = agg.finalize_context();
        assert_approx_eq!(finalized.get(&Context::CG).unwrap(), &0.7); // 1.4 / 2.0
        assert_approx_eq!(finalized.get(&Context::CHG).unwrap(), &0.8); // 0.8 / 1.0
    }

    #[test]
    fn test_region_meth_agg_finalize_strand() {
        let mut agg = RegionMethAgg::new();
        agg.add_cytosine(8.0, Some(10), Context::CG, Strand::Forward);
        agg.add_cytosine(6.0, Some(10), Context::CG, Strand::Forward);
        agg.add_cytosine(4.0, Some(5), Context::CHG, Strand::Reverse);

        let finalized = agg.finalize_strand();
        assert_approx_eq!(finalized.get(&Strand::Forward).unwrap(), &0.7); // 1.4 / 2.0
        assert_approx_eq!(finalized.get(&Strand::Reverse).unwrap(), &0.8); // 0.8 / 1.0
    }

    #[test]
    fn test_region_meth_agg_add_assign() {
        let mut agg1 = RegionMethAgg::new();
        agg1.add_cytosine(8.0, Some(10), Context::CG, Strand::Forward);

        let mut agg2 = RegionMethAgg::new();
        agg2.add_cytosine(6.0, Some(10), Context::CG, Strand::Forward);
        agg2.add_cytosine(4.0, Some(5), Context::CHG, Strand::Reverse);

        agg1 += agg2;

        // Check combined coverage
        assert_eq!(agg1.coverage().get(&10), Some(&2));
        assert_eq!(agg1.coverage().get(&5), Some(&1));

        // Check combined context
        let cg_agg = agg1.context().get(&Context::CG).unwrap();
        assert_approx_eq!(cg_agg.sum(), 1.4); // 0.8 + 0.6
        assert_approx_eq!(cg_agg.count(), 2.0);

        let chg_agg = agg1.context().get(&Context::CHG).unwrap();
        assert_approx_eq!(chg_agg.sum(), 0.8);
        assert_approx_eq!(chg_agg.count(), 1.0);
    }

    #[test]
    fn test_region_meth_agg_add() {
        let mut agg1 = RegionMethAgg::new();
        agg1.add_cytosine(8.0, Some(10), Context::CG, Strand::Forward);

        let mut agg2 = RegionMethAgg::new();
        agg2.add_cytosine(6.0, Some(10), Context::CG, Strand::Forward);

        let result = agg1 + agg2;

        let cg_agg = result.context().get(&Context::CG).unwrap();
        assert_approx_eq!(cg_agg.sum(), 1.4); // 0.8 + 0.6
        assert_approx_eq!(cg_agg.count(), 2.0);
    }

    #[rstest]
    #[case(0.0, 0.0)]
    #[case(5.0, 1.0)]
    #[case(10.5, 3.5)]
    fn test_meth_agg_with_various_values(
        #[case] sum: DensityType,
        #[case] count: DensityType,
    ) {
        let agg = MethAgg::from((sum, count));
        assert_approx_eq!(agg.sum(), sum);
        assert_approx_eq!(agg.count(), count);
        if count > 0.0 {
            assert_eq!(agg.finalize(), sum / count);
        }
    }

    #[rstest]
    #[case(Context::CG, Strand::Forward)]
    #[case(Context::CHG, Strand::Reverse)]
    #[case(Context::CHH, Strand::None)]
    fn test_region_meth_agg_different_contexts_strands(
        #[case] context: Context,
        #[case] strand: Strand,
    ) {
        let mut agg = RegionMethAgg::new();
        agg.add_cytosine(5.0, Some(10), context, strand);

        assert!(agg.context().contains_key(&context));
        assert!(agg.strand().contains_key(&strand));

        let context_agg = agg.context().get(&context).unwrap();
        assert_approx_eq!(context_agg.finalize(), 0.5); // 5.0 / 10

        let strand_agg = agg.strand().get(&strand).unwrap();
        assert_approx_eq!(strand_agg.finalize(), 0.5);
    }
}
