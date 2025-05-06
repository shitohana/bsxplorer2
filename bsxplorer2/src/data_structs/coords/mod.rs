mod contig;
mod gpos;

pub use contig::Contig;
pub use gpos::GenomicPosition;

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;
    use std::ops::{Add, Sub};

    use bio::bio_types::annot::contig::Contig as BioContig;
    use bio::bio_types::annot::pos::SeqPosUnstranded;
    use bio::bio_types::strand::ReqStrand;
    use bio::io::bed::Record as BedRecord;
    use bio_types::annot::loc::Loc;
    use bio_types::strand::NoStrand;

    use super::*;
    use crate::data_structs::enums::Strand;

    // --- GenomicPosition Tests ---

    #[test]
    fn test_genomic_position_from_seq_pos_unstranded_u32() {
        let bio_pos =
            SeqPosUnstranded::new("chr1".to_string(), 100, NoStrand::Unknown);
        let gp: GenomicPosition<String, u32> = bio_pos.into();
        assert_eq!(gp.seqname(), "chr1");
        assert_eq!(gp.position(), 100);
    }

    #[test]
    fn test_genomic_position_from_seq_pos_unstranded_u64() {
        let bio_pos = SeqPosUnstranded::new(
            "chrX".to_string(),
            1_000_000,
            NoStrand::Unknown,
        );
        let gp: GenomicPosition<String, u64> = bio_pos.into();
        assert_eq!(gp.seqname(), "chrX");
        assert_eq!(gp.position(), 1_000_000);
    }

    #[test]
    fn test_genomic_position_into_seq_pos_unstranded_u32() {
        let gp = GenomicPosition::new("chr1".to_string(), 100u32);
        let bio_pos: SeqPosUnstranded = gp.into();
        assert_eq!(bio_pos.refid(), "chr1");
        assert_eq!(bio_pos.pos(), 100);
    }

    #[test]
    fn test_genomic_position_into_seq_pos_unstranded_u64() {
        let gp = GenomicPosition::new("chrX".to_string(), 1_000_000u64);
        let bio_pos: SeqPosUnstranded = gp.into();
        assert_eq!(bio_pos.refid(), "chrX");
        assert_eq!(bio_pos.pos(), 1_000_000);
    }

    #[test]
    fn test_genomic_position_sub_same_seqname_ge() {
        let gp1 = GenomicPosition::new("chr1".to_string(), 100u32);
        let gp2 = GenomicPosition::new("chr1".to_string(), 50u32);
        let result = gp1.sub(gp2);
        assert!(result.is_some());
        let new_pos = result.unwrap();
        assert_eq!(new_pos.seqname(), "chr1");
        assert_eq!(new_pos.position(), 50);
    }

    #[test]
    fn test_genomic_position_sub_same_seqname_lt() {
        let gp1 = GenomicPosition::new("chr1".to_string(), 50u32);
        let gp2 = GenomicPosition::new("chr1".to_string(), 100u32);
        let result = gp1.sub(gp2);
        assert!(result.is_none());
    }

    #[test]
    fn test_genomic_position_sub_different_seqnames() {
        let gp1 = GenomicPosition::new("chr1".to_string(), 100u32);
        let gp2 = GenomicPosition::new("chr2".to_string(), 50u32);
        let result = gp1.sub(gp2);
        assert!(result.is_none());
    }

    #[test]
    fn test_genomic_position_add_same_seqname() {
        let gp1 = GenomicPosition::new("chr1".to_string(), 100u32);
        let gp2 = GenomicPosition::new("chr1".to_string(), 50u32);
        let result = gp1.add(gp2);
        assert!(result.is_some());
        let new_pos = result.unwrap();
        assert_eq!(new_pos.seqname(), "chr1");
        assert_eq!(new_pos.position(), 150);
    }

    #[test]
    fn test_genomic_position_add_different_seqnames() {
        let gp1 = GenomicPosition::new("chr1".to_string(), 100u32);
        let gp2 = GenomicPosition::new("chr2".to_string(), 50u32);
        let result = gp1.add(gp2);
        assert!(result.is_none());
    }

    #[test]
    #[should_panic(
        expected = "Cannot compare genomic positions with different seqnames"
    )]
    fn test_genomic_position_cmp_different_seqnames_panics() {
        let gp1 = GenomicPosition::new("chr1".to_string(), 100u32);
        let gp2 = GenomicPosition::new("chr2".to_string(), 100u32);
        // This should panic because partial_cmp returns None, and cmp expects
        // Some
        let _ = gp1.cmp(&gp2);
    }

    #[test]
    fn test_genomic_position_eq_and_ne() {
        let gp1a = GenomicPosition::new("chr1".to_string(), 100u32);
        let gp1b = GenomicPosition::new("chr1".to_string(), 100u32);
        let gp2 = GenomicPosition::new("chr1".to_string(), 200u32);
        let gp3 = GenomicPosition::new("chr2".to_string(), 100u32);

        assert!(gp1a.eq(&gp1b));
        assert_eq!(gp1a, gp1b); // Uses PartialEq
        assert!(!gp1a.ne(&gp1b));

        assert!(!gp1a.eq(&gp2));
        assert_ne!(gp1a, gp2);
        assert!(gp1a.ne(&gp2));

        assert!(!gp1a.eq(&gp3));
        assert_ne!(gp1a, gp3);
        assert!(gp1a.ne(&gp3));
    }

    #[test]
    fn test_genomic_position_display() {
        let gp = GenomicPosition::new("chr1".to_string(), 12345u32);
        assert_eq!(format!("{}", gp), "chr1:12345");
    }

    // --- Contig Tests ---

    #[test]
    #[should_panic(
        expected = "Start position must be less than or equal to end position"
    )]
    fn test_contig_new_invalid_range_panics() {
        Contig::new("chr1".to_string(), 100u32, 50u32, Strand::None);
    }

    #[test]
    fn test_contig_strand() {
        let contig = Contig::new("chr1".to_string(), 1u32, 10, Strand::Forward);
        assert_eq!(contig.strand(), Strand::Forward);
    }

    #[test]
    fn test_contig_length() {
        let contig = Contig::new("chr1".to_string(), 10u32, 100, Strand::None);
        assert_eq!(contig.length(), 90);
    }

    #[test]
    fn test_contig_extend_upstream() {
        let mut contig =
            Contig::new("chr1".to_string(), 100u32, 200u32, Strand::None);
        contig.extend_upstream(50u32);
        assert_eq!(contig.start(), 50);
        assert_eq!(contig.end(), 200);

        // Test saturating behavior
        let mut contig_saturate =
            Contig::new("chr1".to_string(), 10u32, 20u32, Strand::None);
        contig_saturate.extend_upstream(50u32);
        assert_eq!(contig_saturate.start(), 0); // Start is 0 due to saturating_sub
        assert_eq!(contig_saturate.end(), 20);
    }

    #[test]
    fn test_contig_extend_downstream() {
        let mut contig =
            Contig::new("chr1".to_string(), 100u32, 200u32, Strand::None);
        contig.extend_downstream(50u32);
        assert_eq!(contig.start(), 100);
        assert_eq!(contig.end(), 250);
    }

    #[test]
    fn test_contig_set_start() {
        let mut contig =
            Contig::new("chr1".to_string(), 100u32, 200u32, Strand::None);
        contig.set_start(150u32);
        assert_eq!(contig.start(), 150);
        assert_eq!(contig.end(), 200);
    }

    #[test]
    fn test_contig_set_end() {
        let mut contig =
            Contig::new("chr1".to_string(), 100u32, 200u32, Strand::None);
        contig.set_end(250u32);
        assert_eq!(contig.start(), 100);
        assert_eq!(contig.end(), 250);
    }

    #[test]
    fn test_contig_from_range_genomic_position() {
        let start = GenomicPosition::new("chr1".to_string(), 100u32);
        let end = GenomicPosition::new("chr1".to_string(), 200u32);
        let range = start..end;
        let contig: Contig<String, u32> = range.into();
        assert_eq!(contig.seqname(), "chr1");
        assert_eq!(contig.start(), 100);
        assert_eq!(contig.end(), 200);
        assert_eq!(contig.strand(), Strand::None);
    }

    #[test]
    #[should_panic(
        expected = "Start and end positions must have the same sequence name"
    )]
    fn test_contig_from_range_genomic_position_different_seqnames_panics() {
        let start = GenomicPosition::new("chr1".to_string(), 100u32);
        let end = GenomicPosition::new("chr2".to_string(), 200u32);
        let range = start..end;
        let _contig: Contig<String, u32> = range.into(); // This should panic
    }

    #[test]
    fn test_contig_from_bed_record() {
        let mut bed_record = BedRecord::new();
        bed_record.set_chrom("chr1");
        bed_record.set_start(100);
        bed_record.set_end(200);

        let contig: Contig<String, u64> = bed_record.into();
        assert_eq!(contig.seqname(), "chr1");
        assert_eq!(contig.start(), 100);
        assert_eq!(contig.end(), 200);
        assert_eq!(contig.strand(), Strand::None);

        let mut bed_record_no_strand = BedRecord::new();
        bed_record_no_strand.set_chrom("chr2");
        bed_record_no_strand.set_start(50);
        bed_record_no_strand.set_end(150);

        let contig_no_strand: Contig<String, u64> = bed_record_no_strand.into();
        assert_eq!(contig_no_strand.strand(), Strand::None);
    }

    #[test]
    fn test_contig_into_bed_record() {
        let contig =
            Contig::new("chr1".to_string(), 100u64, 200u64, Strand::Reverse);
        let bed_record: BedRecord = contig.into();
        assert_eq!(bed_record.chrom(), "chr1");
        assert_eq!(bed_record.start(), 100);
        assert_eq!(bed_record.end(), 200);

        let contig_no_strand =
            Contig::new("chr2".to_string(), 50u64, 150u64, Strand::None);
        let bed_record_no_strand: BedRecord = contig_no_strand.into();
        assert_eq!(bed_record_no_strand.strand(), None);
    }

    #[test]
    fn test_contig_from_bio_contig() {
        let bio_contig = BioContig::new(
            "chr1".to_string(),
            100,
            100 as usize,
            Some(ReqStrand::Reverse),
        ); // start 100, length 100 -> end 200
        let contig: Contig<String, u64> = bio_contig.into();
        assert_eq!(contig.seqname(), "chr1");
        assert_eq!(contig.start(), 100);
        assert_eq!(contig.end(), 200);
        assert_eq!(contig.strand(), Strand::Reverse);

        let bio_contig_no_strand = BioContig::new(
            "chr2".to_string(),
            50,
            50 as usize,
            None::<ReqStrand>,
        ); // start 50, length 50 -> end 100
        let contig_no_strand: Contig<String, u64> = bio_contig_no_strand.into();
        assert_eq!(contig_no_strand.seqname(), "chr2");
        assert_eq!(contig_no_strand.start(), 50);
        assert_eq!(contig_no_strand.end(), 100);
        assert_eq!(contig_no_strand.strand(), Strand::None);
    }

    #[test]
    fn test_contig_into_bio_contig() {
        let contig =
            Contig::new("chr1".to_string(), 100u64, 200u64, Strand::Forward); // start 100, end 200 -> length 100
        let bio_contig: BioContig<String, Option<ReqStrand>> = contig.into();
        assert_eq!(bio_contig.refid(), "chr1");
        assert_eq!(bio_contig.start(), 100);
        assert_eq!(bio_contig.length(), 100);
        assert_eq!(bio_contig.strand(), Some(ReqStrand::Forward));

        let contig_no_strand =
            Contig::new("chr2".to_string(), 50u64, 100u64, Strand::None); // start 50, end 100 -> length 50
        let bio_contig_no_strand: BioContig<String, Option<ReqStrand>> =
            contig_no_strand.into();
        assert_eq!(bio_contig_no_strand.strand(), None);
    }

    #[test]
    fn test_contig_partial_cmp() {
        let c1 = Contig::new("chr1".to_string(), 100u32, 200u32, Strand::None);
        let c2 = Contig::new("chr1".to_string(), 250u32, 300u32, Strand::None); // c1 before c2
        let c3 = Contig::new("chr1".to_string(), 50u32, 90u32, Strand::None); // c3 before c1
        let c4 = Contig::new("chr1".to_string(), 150u32, 250u32, Strand::None); // c4 intersects c1
        let c5 = Contig::new("chr2".to_string(), 100u32, 200u32, Strand::None); // Different chromosome

        assert_eq!(c1.partial_cmp(&c2), Some(Ordering::Less)); // c1.end <= c2.start is false, c1.start >= c2.end is false, but c1.end
                                                               // <= c2.start
        assert_eq!(c2.partial_cmp(&c1), Some(Ordering::Greater)); // c2.start >= c1.end

        assert_eq!(c1.partial_cmp(&c3), Some(Ordering::Greater)); // c1.start >= c3.end
        assert_eq!(c3.partial_cmp(&c1), Some(Ordering::Less)); // c3.end <= c1.start

        assert_eq!(c1.partial_cmp(&c4), None); // Intersects
        assert_eq!(c4.partial_cmp(&c1), None); // Intersects

        assert_eq!(c1.partial_cmp(&c5), None); // Different chromosome
        assert_eq!(c5.partial_cmp(&c1), None); // Different chromosome
    }

    #[test]
    fn test_contig_display() {
        let contig =
            Contig::new("chrX".to_string(), 1000u32, 2000u32, Strand::Forward);
        assert_eq!(format!("{}", contig), "chrX:1000-2000 (+)");
    }
}
