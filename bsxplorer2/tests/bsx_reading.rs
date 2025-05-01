mod common;
use std::path::PathBuf;
use std::fs::File;

use bsxplorer2::data_structs::batch::EncodedBsxBatch;
use bsxplorer2::{
    data_structs::{
        batch::BsxBatchMethods,
        coords::Contig,
        enums::Strand,
    },
    io::bsx::{BsxFileReader, RegionReader},
};
use rstest::{fixture, rstest};
use bio::io::gff::Reader as GffReader;

#[fixture]
fn test_bsxreader() -> BsxFileReader<File> {
    let reader = BsxFileReader::new(
        File::open(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/report.bsx")).expect("Error opening test report file")
    );
    reader
}

#[fixture]
fn test_gffreader() -> GffReader<File> {
    GffReader::from_file(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/annot.gff"),
        bio::io::gff::GffType::GFF3
    ).expect("Error opening test GFF file")
}

#[fixture]
fn test_regionreader(test_bsxreader: BsxFileReader<File>) -> RegionReader<File> {
    RegionReader::new(test_bsxreader)
}

/// Test that the reader will return error, if required batch has already been processed.
#[rstest]
#[should_panic]
fn unsorted_pos(
    mut test_regionreader: RegionReader<File>,
) {
    let _ = test_regionreader.query(Contig::new("NC_003070.9", 50_000, 50_899, Strand::None)).unwrap();
    let _ = test_regionreader.query(Contig::new("NC_003070.9", 5174, 5326, Strand::None)).unwrap();
}

/// Test that the reader will return error, if such chr does not exist
#[rstest]
#[should_panic]
fn unexistent_chr(
    mut test_regionreader: RegionReader<File>,
) {
    let _ = test_regionreader.query(Contig::new("NC_123456.9", 50_000, 50_899, Strand::None)).unwrap();
}

#[rstest]
fn basic_reading(
    mut test_regionreader: RegionReader<File>,
) -> anyhow::Result<()> {
    let test_contigs = [
        Contig::new("NC_003070.9".to_string(), 50_000, 51_000, Strand::None),
        Contig::new("NC_003070.9".to_string(), 52_000, 53_000, Strand::None),
        Contig::new("NC_003070.9".to_string(), 101_000, 150_000, Strand::None),
        Contig::new("NC_003070.9".to_string(), 151_000, 151_001, Strand::None)
    ];
    assert!(check_batch(&test_contigs[0], &test_regionreader.query(test_contigs[0].clone())?.unwrap())?);
    assert!(check_batch(&test_contigs[1], &test_regionreader.query(test_contigs[1].clone())?.unwrap())?);
    assert!(check_batch(&test_contigs[2], &test_regionreader.query(test_contigs[2].clone())?.unwrap())?);
    assert!(check_batch(&test_contigs[3], &test_regionreader.query(test_contigs[3].clone())?.unwrap())?);

    Ok(())
}

#[rstest]
fn sorted_reading(
    mut test_regionreader: RegionReader<File>,
) -> anyhow::Result<()> {
    let test_contigs = [
        Contig::new("NC_003070.9".to_string(), 52_000, 53_000, Strand::None),
        Contig::new("NC_003070.9".to_string(), 50_000, 51_000, Strand::None),
        Contig::new("NC_003070.9".to_string(), 151_000, 151_001, Strand::None),
        Contig::new("NC_003070.9".to_string(), 101_000, 150_000, Strand::None),
    ];

    let sorted = test_regionreader.index()?.sort(test_contigs);

    let expected_order = [
        Contig::new("NC_003070.9".to_string(), 50_000, 51_000, Strand::None),
        Contig::new("NC_003070.9".to_string(), 52_000, 53_000, Strand::None),
        Contig::new("NC_003070.9".to_string(), 101_000, 150_000, Strand::None),
        Contig::new("NC_003070.9".to_string(), 151_000, 151_001, Strand::None)
        ].to_vec();

    assert_eq!(sorted, expected_order);

    assert!(check_batch(&sorted[0], &test_regionreader.query(sorted[0].clone())?.unwrap())?);
    assert!(check_batch(&sorted[1], &test_regionreader.query(sorted[1].clone())?.unwrap())?);
    assert!(check_batch(&sorted[2], &test_regionreader.query(sorted[2].clone())?.unwrap())?);
    assert!(check_batch(&sorted[3], &test_regionreader.query(sorted[3].clone())?.unwrap())?);

    Ok(())
}

#[rstest]
fn different_chr_reading(
    mut test_regionreader: RegionReader<File>,
) -> anyhow::Result<()> {
    let test_contigs = [
        Contig::new("NC_003074.8".to_string(), 52_000, 53_000, Strand::None),
        Contig::new("NC_003071.7".to_string(), 50_000, 51_000, Strand::None),
        Contig::new("NC_003070.9".to_string(), 151_000, 151_001, Strand::None),
        Contig::new("NC_003070.9".to_string(), 101_000, 150_000, Strand::None),
    ];

    let sorted = test_regionreader.index()?.sort(test_contigs);

    let expected_order = [
        Contig::new("NC_003070.9".to_string(), 101_000, 150_000, Strand::None),
        Contig::new("NC_003070.9".to_string(), 151_000, 151_001, Strand::None),
        Contig::new("NC_003074.8".to_string(), 52_000, 53_000, Strand::None),
        Contig::new("NC_003071.7".to_string(), 50_000, 51_000, Strand::None),
        ].to_vec();

    assert_eq!(sorted, expected_order);

    assert!(check_batch(&sorted[0], &test_regionreader.query(sorted[0].clone())?.unwrap())?);
    assert!(check_batch(&sorted[1], &test_regionreader.query(sorted[1].clone())?.unwrap())?);
    assert!(check_batch(&sorted[2], &test_regionreader.query(sorted[2].clone())?.unwrap())?);
    assert!(check_batch(&sorted[3], &test_regionreader.query(sorted[3].clone())?.unwrap())?);

    Ok(())
}

fn check_batch(contig: &Contig<String, u32>, batch: &EncodedBsxBatch) -> anyhow::Result<bool> {
    if let Some(b_contig) = batch.as_contig()? {
        Ok(b_contig.is_in(contig))
    } else {
        Ok(true)
    }
}
