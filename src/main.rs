#![feature(btree_cursors)]
#![feature(ascii_char)]

mod bsxipc;
mod genome;
mod python;
mod io;
mod utils;
mod bsx_reader;
mod bsx_writer;

use polars::prelude::*;
use rayon::prelude::*;

fn main() {
    let annotation = genome::AnnotationBuilder::from_gff(
        "/Users/shitohana/Documents/CX_reports/old/A_thaliana_genomic.gff"
    ).finish().unwrap();
    println!("{}", annotation);
}