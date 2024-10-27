#![feature(btree_cursors)]

mod genome;
mod report;
mod sequence;
mod utils;

use genome::{get_index, Annotation};
use polars::prelude::{col, lit};
use std::time::Instant;
use utils::types::Context;

fn main() {
    // let start = Instant::now();
    // ReportType::BISMARK.convert_report(
    //     "/Users/shitohana/Documents/CX_reports/old/A_thaliana.txt",
    //     "/Users/shitohana/Desktop/RustProjects/bsxplorer2/res.ipc",
    //     "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa",
    //     Some(10000),
    //     Some(false),
    //     None,
    // );
    // let duration = start.elapsed();
    // println!("{:?}", duration);
    let annot =
        Annotation::from_gff("/Users/shitohana/Documents/CX_reports/old/A_thaliana_genomic.gff");
    let index_map = get_index("/Users/shitohana/Desktop/RustProjects/bsxplorer2/res.ipc").unwrap();
    let start = Instant::now();
    let annot_finished = annot
        .filter(col("type").eq(lit("CDS")))
        .align(
            "/Users/shitohana/Desktop/RustProjects/bsxplorer2/res.ipc",
            Context::CG,
        )
        .finish()
        .unwrap();

    println!("{:?}", annot_finished);
}
