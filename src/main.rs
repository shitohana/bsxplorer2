#![allow(warnings)]

extern crate core;

use std::fs::File;
use std::io::BufWriter;
use log::info;
use std::thread;
use bsx_rs::io::bsx::writer::BsxIpcWriter;
use bsx_rs::io::report::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use bsx_rs::io::report::reader::ReportReaderBuilder;
use bsx_rs::io::report::schema::ReportTypeSchema;

extern crate pretty_env_logger;

struct ThreadTest {
    read_thread: thread::JoinHandle<()>,
}

/*
Writing
-------
1. Reader should read and process batches in separate threads
2. Reader should have method for processing batch into proper format
3. Reader should have internal cache for sorting batches (size > batches per read)
4. Batches should be marked whether they are end of chromosome or not
5. Alignment with reference genome should be obligatory for bedgraph or coverage
   and optional for other

Reading
-------
1. Reader should check consistency of batches among files
*/

fn main() {
    pretty_env_logger::init_timed();
    info!("Starting BSX RS Test");
    // ConvertReportOptions::default()
    //     .with_chunk_size(10_000)
    //     .with_batch_size(10_000)
    //     .with_batch_per_read(50)
    //     .finish(
    //         ReportType::BISMARK,
    //         "/Users/shitohana/Documents/CX_reports/old/A_thaliana.txt",
    //         "/Users/shitohana/Desktop/RustProjects/bsxplorer2_dev/bismark.ipc",
    //         "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa",
    //         "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa.fai",
    //     ).unwrap();
    // println!("Report written");
    
    let reader = ReportReaderBuilder::default()
        .with_report_type(ReportTypeSchema::Bismark)
        .try_finish(File::open("/Users/shitohana/Documents/CX_reports/old/A_thaliana.txt").unwrap())
        .unwrap();
    
    let mut writer = BsxIpcWriter::try_from_path_and_fai(
        BufWriter::new(File::create("/Users/shitohana/Desktop/RustProjects/bsxplorer2_dev/test.ipc").unwrap()),
        "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa.fai".into()
    ).unwrap();
    
    for batch in reader {
        let encoded = EncodedBsxBatch::encode(batch, writer.get_chr_dtype()).unwrap();
        let test = encoded.trim_region::<u32>(&(encoded.first_position().unwrap() + 3 >> encoded.last_position().unwrap() - 5).unwrap()).unwrap();
        println!("{}", encoded.data());
        println!("{}", test.data());
        writer.write_batch(encoded).unwrap()
    }
    
    writer.close().unwrap();
    
    // let reader = BSXReader::new(
    //     File::open("/Users/shitohana/Desktop/RustProjects/bsxplorer2_dev/bismark.ipc").unwrap(),
    //     None,
    // );
    // let annotation = AnnotationBuilder::from_gff(
    //     "/Users/shitohana/Documents/CX_reports/old/A_thaliana_genomic.gff",
    // )
    // .filter(col("type").eq(lit("gene")))
    // .finish()
    // .unwrap();
    // let region_reader = RegionsDataReader::new(reader, annotation);
    //
    // let start = Instant::now();
    //
    // for reg_data in {
    //     region_reader
    //         .into_iter()
    //         .with_batch_per_read(16)
    //         .with_filters(ReadFilters::new(Some(Context::CG), Some(Strand::None)))
    // } {
    //     // println!("{}", reg_data);
    // }
    // let duration = start.elapsed();
    // println!("Time elapsed is: {:?}", duration);
}
