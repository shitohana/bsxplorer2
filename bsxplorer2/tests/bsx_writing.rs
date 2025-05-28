mod common;
use std::fs::File;
use std::path::PathBuf;

use bsxplorer2::data_structs::typedef::BsxSmallStr;
use bsxplorer2::io::bsx::{
    BatchIndex,
    BsxFileReader,
    BsxFileWriter,
};
use itertools::Itertools;
use rstest::{
    fixture,
    rstest,
};
use tempfile::NamedTempFile;

#[fixture]
fn test_bsxreader() -> BsxFileReader {
    let reader = BsxFileReader::try_new(
        File::open(
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/report.bsx"),
        )
        .expect("Error opening test report file"),
    )
    .expect("Failed to create reader");
    reader
}

#[rstest]
fn write_bsx(mut test_bsxreader: BsxFileReader) -> anyhow::Result<()> {
    let index = BatchIndex::from_reader(&mut test_bsxreader)?;
    let sink_tempfile = NamedTempFile::new()?;
    let mut writer = BsxFileWriter::try_new(
        sink_tempfile.reopen()?,
        index
            .get_chr_order()
            .into_iter()
            .map(BsxSmallStr::to_string)
            .collect_vec(),
        None,
        None,
    )?;

    for batch in test_bsxreader.iter() {
        writer.write_batch(batch?)?;
    }
    writer.close()?;

    let mut written_reader = BsxFileReader::try_new(sink_tempfile.reopen()?)?;

    for (orig_batch, read_batch) in test_bsxreader.iter().zip(written_reader.iter()) {
        let orig_batch = orig_batch.unwrap();
        let read_batch = read_batch.unwrap();

        assert_eq!(orig_batch, read_batch)
    }

    Ok(())
}
