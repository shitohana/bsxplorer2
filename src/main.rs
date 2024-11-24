extern crate core;

use bsx_rs::genome::AnnotationBuilder;
use bsx_rs::io::bsx::bsx_reader::BSXReader;
use bsx_rs::io::bsx::region_reader::{ReadFilters, RegionsDataReader};
use bsx_rs::utils::types::{Context, Strand};
use polars::prelude::*;
use std::fs::File;
use std::time::Instant;
use bsx_rs::io::bsx::write::ConvertReportOptions;
use bsx_rs::io::report::types::ReportType;

extern crate pretty_env_logger;

fn main() {
    pretty_env_logger::init_timed();
    ConvertReportOptions::default()
        .with_chunk_size(10_000)
        .with_batch_size(10_000)
        .with_batch_per_read(50)
        .finish(
            ReportType::BISMARK,
            "/Users/shitohana/Documents/CX_reports/old/A_thaliana.txt",
            "/Users/shitohana/Desktop/RustProjects/bsxplorer2_dev/bismark.ipc",
            "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa",
            "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa.fai",
        ).unwrap();
    println!("Report written");

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
