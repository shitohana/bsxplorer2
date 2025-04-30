mod common;
use std::collections::HashSet;
use std::error::Error;
use std::io::Cursor;
use std::path::PathBuf;
use std::{collections::BTreeMap, fs::File};

use bsxplorer2::data_structs::batch::{BsxBatch, EncodedBsxBatch};
use bsxplorer2::io::bsx::BsxIpcWriter;
use bsxplorer2::io::report::{ReportReaderBuilder, ReportTypeSchema};
use bsxplorer2::{
    data_structs::{
        batch::{BsxBatchBuilder, BsxBatchMethods},
        coords::Contig,
        enums::Strand,
    },
    io::bsx::{BsxFileReader, RegionReader},
};
use common::DemoReportBuilder;
use itertools::{izip, Itertools};
use polars::prelude::IpcCompression;
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
        .map(|batch| (batch.chr_val().unwrap().to_owned(), batch))
        .collect::<BTreeMap<_, _>>();

    // Reading
    let reader = BsxFileReader::new(Cursor::new(&mut sink_buffer));
    let read_batches = reader
        .into_iter()
        .flatten()
        .map(|batch| BsxBatchBuilder::decode_batch(batch))
        .flatten()
        .into_group_map_by(|batch| batch.chr_val().unwrap().to_string());
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

#[test]
fn test_bsx_index() -> anyhow::Result<()> {
    const CHR_LEN: usize = 123_456;
    let demo_builder: DemoReportBuilder<StdRng> =
        DemoReportBuilder::new(CHR_LEN, 30, 15.0, 0.3, Some(42));

    const NUM_CHRS: usize = 10;
    const BATCH_SIZE: usize = 10_000;

    let mut sink_buffer = Vec::<u8>::new();
    let orig_batches = demo_builder
        .write_bsx(Cursor::new(&mut sink_buffer), NUM_CHRS, BATCH_SIZE)?
        .into_iter()
        .map(|batch| (batch.chr_val().unwrap().to_owned(), batch))
        .collect::<BTreeMap<_, _>>();

    let test_region = Contig::new(
        orig_batches
            .keys()
            .next()
            .unwrap()
            .to_owned(),
        5000 as u32,
        55_000 as u32,
        Strand::None,
    );

    let true_test_region = orig_batches
        .get(&test_region.seqname().to_string())
        .cloned()
        .unwrap()
        .lazy()
        .filter_pos_gt(test_region.start() - 1)
        .filter_pos_lt(test_region.end())
        .collect()?;

    let mut reader = BsxFileReader::new(Cursor::new(&sink_buffer));

    let query_res = reader.query(&test_region)?.unwrap();
    let read_test_region = BsxBatchBuilder::decode_batch(query_res)?;

    let orig_nodensity = true_test_region
        .data()
        .drop("density")?;
    let read_nodensity = read_test_region
        .data()
        .drop("density")?;
    assert_eq!(orig_nodensity, read_nodensity);

    let mut region_reader = RegionReader::new(reader);
    let query_res = region_reader
        .query(test_region)?
        .unwrap();
    let read_test_region = BsxBatchBuilder::decode_batch(query_res)?;

    let orig_nodensity = true_test_region
        .data()
        .drop("density")?;
    let read_nodensity = read_test_region
        .data()
        .drop("density")?;
    assert_eq!(orig_nodensity, read_nodensity);

    Ok(())
}

#[test]
fn convert_report() -> Result<(), Box<dyn Error>> {
    let file_handle =
        File::open("/Users/shitohana/Documents/CX_reports/old/A_thaliana.txt")?;
    let fai_path: PathBuf =
        "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa.fai".into();

    let reader = ReportReaderBuilder::<BsxBatch>::default()
        .with_report_type(ReportTypeSchema::Bismark)
        .with_chunk_size(10_000)
        .with_fai_path(fai_path.clone())
        .build_from_handle(Box::new(file_handle))?;

    let file_sink = File::create(
        "/Users/shitohana/Documents/CX_reports/old/A_thaliana.bsx",
    )?;
    let mut writer = BsxIpcWriter::try_from_sink_and_fai(
        file_sink,
        fai_path,
        Some(IpcCompression::ZSTD),
        None,
    )?;

    let mut selected_chroms = HashSet::new();
    for batch in reader {
        let batch = batch?;
        let chr = batch.chr_val()?.to_string();

        if !selected_chroms.contains(&chr) && selected_chroms.len() >= 3{
            break;
        } else {
            selected_chroms.insert(chr);
        }
        writer.write_batch(batch)?;
    }
    writer.close()?;
    println!("Conversion completed successfully.");
    println!("Selected chromosomes: {:?}", selected_chroms);

    Ok(())
}
