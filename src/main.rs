use std::fs::File;
use std::io::{BufReader, BufWriter};
use polars::prelude::*;
use bsx_rs::genome::AnnotationBuilder;
use bsx_rs::io::bsx::read::{df_to_regions, RegionReader};
use bsx_rs::io::bsx::writer::BsxIpcWriter;
use bsx_rs::io::report::bsx_batch::EncodedBsxBatch;
use bsx_rs::io::report::reader::ReportReaderBuilder;
use bsx_rs::io::report::schema::ReportTypeSchema;

fn main() {
    // let mut ipc_writer = BsxIpcWriter::try_from_sink_and_fai(
    //     BufWriter::new(File::create("/Users/shitohana/Desktop/RustProjects/bsxplorer2_dev/test.ipc").unwrap()),
    //     "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa.fai".into()
    // ).unwrap();
    // 
    // let report_reader = ReportReaderBuilder::default()
    //     .with_report_type(ReportTypeSchema::Bismark)
    //     .try_finish(
    //         BufReader::new(File::open("/Users/shitohana/Documents/CX_reports/old/A_thaliana.txt").unwrap())
    //     ).unwrap();
    // 
    // for batch in report_reader {
    //     ipc_writer.write_batch(EncodedBsxBatch::encode(batch, ipc_writer.get_chr_dtype()).unwrap())
    //         .expect("Can't write");
    // }
    // 
    // ipc_writer.close().unwrap();
    
    let annotation = AnnotationBuilder::from_gff("/Users/shitohana/Documents/CX_reports/old/A_thaliana_genomic.gff")
        .filter(col("type").eq(lit("gene")))
        .finish().unwrap();
    let regions_of_interest = df_to_regions(&annotation, None, None, None).unwrap();
    
    let regions_reader = RegionReader::try_new(
        BufReader::new(File::open("/Users/shitohana/Desktop/RustProjects/bsxplorer2_dev/test.ipc").unwrap()),
        &regions_of_interest,
    ).unwrap();
    
    for region_data in regions_reader.into_iter() {
        println!("Region: {:?}", region_data);
    }
    
    println!("{}", regions_of_interest.first().unwrap());
}