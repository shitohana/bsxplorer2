use bsx_rs::bsx_batch::BsxBatchMethods;
use bsx_rs::genome::AnnotationBuilder;
use bsx_rs::io::bsx::region_read::{df_to_regions, RegionReader};
use polars::prelude::*;
use std::fs::File;
use std::io::{BufReader, Read};
use std::time::Instant;

fn main() {
    // let mut ipc_writer = BsxIpcWriter::try_from_sink_and_fai(
    //     BufWriter::new(File::create("/Users/shitohana/Desktop/RustProjects/bsxplorer2_dev/test.ipc").unwrap()),
    //     "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa.fai".into()
    // ).unwrap();
    //
    // let report_reader = ReportReaderBuilder::default()
    //     .with_report_type(ReportTypeSchema::Bismark)
    //     .with_chunk_size(5000)
    //     .try_finish(
    //         BufReader::new(File::open("/Users/shitohana/Documents/CX_reports/old/A_thaliana.txt").unwrap())
    //     ).unwrap();
    //
    // for batch in report_reader {
    //     ipc_writer.write_batch(batch)
    //         .expect("Can't write");
    // }
    //
    // ipc_writer.close().unwrap();

    let annotation = AnnotationBuilder::from_gff(
        "/Users/shitohana/Documents/CX_reports/old/A_thaliana_genomic.gff",
    )
    .filter(col("type").eq(lit("gene")))
    .finish()
    .unwrap();
    let regions_of_interest = df_to_regions(&annotation, None, None, None).unwrap();

    for _ in 0..3 {
        let regions_reader = RegionReader::try_new(
            BufReader::new(
                File::open("/Users/shitohana/Desktop/RustProjects/bsxplorer2_dev/test.ipc")
                    .unwrap(),
            ),
            &regions_of_interest,
        )
        .unwrap();

        let start = Instant::now();
        for region_data in regions_reader.into_iter() {
            // println!("Region: {:?}", region_data.0);
        }
        println!("Elapsed time: {:?}", start.elapsed());
    }
}
