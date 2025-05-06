mod common;
use std::collections::BTreeMap;
use std::io::Cursor;

use bsxplorer2::data_structs::batch::{BsxBatchBuilder, BsxBatchMethods};
use bsxplorer2::io::bsx::BsxFileReader;
use common::DemoReportBuilder;
use itertools::{izip, Itertools};
use rand::rngs::StdRng;
use rstest::{fixture, rstest};

#[fixture]
fn report_builder() -> DemoReportBuilder<StdRng> {
    DemoReportBuilder::new(123_456, 20, 15.0, 0.5, Some(42))
}

#[rstest]
#[case::small_batch(1_000usize, 3)]
#[case::medium_batch(10_000usize, 3)]
#[case::large_batch(200_000usize, 3)]
#[case::single_chr(1_000usize, 1)]
fn read_write(
    report_builder: DemoReportBuilder<StdRng>,
    #[case] batch_size: usize,
    #[case] num_chrs: usize,
) -> anyhow::Result<()> {
    // Writing
    let mut sink_buffer = Vec::<u8>::new();
    let orig_batches = report_builder
        .write_bsx(Cursor::new(&mut sink_buffer), num_chrs, batch_size)?
        .into_iter()
        .map(|batch| (batch.chr_val().to_owned(), batch))
        .collect::<BTreeMap<_, _>>();

    // Reading
    let reader = BsxFileReader::new(Cursor::new(&mut sink_buffer));
    let read_batches = reader
        .into_iter()
        .flatten()
        .map(|batch| BsxBatchBuilder::decode_batch(batch))
        .flatten()
        .into_group_map_by(|batch| batch.chr_val().to_string());
    let read_batches: BTreeMap<_, _> = read_batches
        .into_iter()
        .map(|(k, v)| (k, BsxBatchBuilder::concat(v).unwrap()))
        .collect();

    for ((chr1, batch1), (chr2, batch2)) in
        izip!(orig_batches.into_iter(), read_batches.into_iter())
    {
        assert_eq!(chr1, chr2);
        let orig_nodensity = batch1.data().drop("density")?;
        let read_nodensity = batch2.data().drop("density")?;
        assert_eq!(orig_nodensity, read_nodensity);
    }

    Ok(())
}
