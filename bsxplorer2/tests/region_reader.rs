use std::fs::File;
use std::path::PathBuf;
use std::sync::Arc;

use bsxplorer2::data_structs::annotation::HcAnnotStore;
use bsxplorer2::data_structs::batch::BsxBatch;
use bsxplorer2::data_structs::coords::Contig;
use bsxplorer2::data_structs::Context;
use bsxplorer2::io::bsx::{
    BatchIndex,
    BsxFileReader,
    RegionReader,
};
use rstest::{
    fixture,
    rstest,
};

#[fixture]
fn bsx_reader() -> BsxFileReader {
    BsxFileReader::try_new(
        File::open(
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/report.bsx"),
        )
        .expect("Error opening test report file"),
    )
    .expect("Failed to create reader")
}

#[fixture]
#[once]
fn index(mut bsx_reader: BsxFileReader) -> BatchIndex {
    BatchIndex::from_reader(&mut bsx_reader).expect("Failed to create index")
}

#[fixture]
fn reader(
    bsx_reader: BsxFileReader,
    index: &BatchIndex,
) -> RegionReader {
    RegionReader::new(bsx_reader, index.clone(), None)
}

#[fixture]
#[once]
fn contigs(index: &BatchIndex) -> Vec<Contig> {
    let annot = HcAnnotStore::from_bed(
        File::open(
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/genes.bed"),
        )
        .expect("Error opening test contigs file"),
    )
    .expect("Error parsing test contigs file");

    let contigs_iter = annot.values().map(|e| e.contig().clone());
    index.sort(contigs_iter).collect()
}

#[rstest]
fn test_fixtures(
    _bsx_reader: BsxFileReader,
    _index: &BatchIndex,
    _reader: RegionReader,
    _contigs: &[Contig],
) {
}

#[rstest]
fn test_reader_basic(
    mut reader: RegionReader,
    contigs: &[Contig],
) {
    let mut count = 0;
    for batch in reader.iter_contigs(contigs) {
        assert!(batch.is_ok());
        count += 1;
    }
    let seqnames = reader.index().get_chr_order();
    assert_eq!(
        count,
        contigs
            .iter()
            .filter(|c| seqnames.contains(c.seqname()))
            .count()
    );
}

#[rstest]
fn test_reader_preprocess(
    mut reader: RegionReader,
    contigs: &[Contig],
) {
    let mut count = 0;
    reader.set_preprocess_fn(Some(Arc::new(|batch: BsxBatch| {
        batch
            .lazy()
            .filter_context(Context::CG)
            .collect()
            .map_err(|e| anyhow::anyhow!(e))
    })));
    for batch in reader.iter_contigs(contigs) {
        assert!(batch.is_ok());
        count += 1;
    }
    let seqnames = reader.index().get_chr_order();
    assert_eq!(
        count,
        contigs
            .iter()
            .filter(|c| seqnames.contains(c.seqname()))
            .count()
    );
}
