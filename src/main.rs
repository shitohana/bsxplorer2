#![feature(btree_cursors)]

mod bsxipc;
mod genome;
mod python;
mod report;
mod sequence;
mod utils;
use bsxipc::*;
use genome::Annotation;
use itertools::Itertools;
use polars::prelude::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{Read, Seek};


fn main() {
    let annotation_path = "/Users/shitohana/Documents/CX_reports/old/A_thaliana_genomic.gff";
    let ipc_path = "/Users/shitohana/Desktop/RustProjects/bsxplorer2/res.ipc";
    let file = File::open(ipc_path).unwrap();
    let bsx_file = BSXFile::new(file, Default::default()).unwrap();

    let mut annotation = Annotation::from_gff(annotation_path)
        .filter(col("type").eq(lit("gene")))
        .finish()
        .unwrap();


    let annotation = bsx_file
        .add_index(&mut annotation, None)
        .unwrap();
    println!("{:?}", annotation);

    let index = annotation
        .column("index")
        .unwrap()
        .list()
        .unwrap()
        .into_iter()
        .filter_map(|s| match s {
            Some(series) => {
                let arr: ChunkedArray<UInt32Type> = series.u32().unwrap().to_owned();
                Some(
                    arr.into_iter()
                        .flatten()
                        .collect::<Vec<u32>>(),
                )
            }
            _ => None,
        })
        .collect::<Vec<Vec<u32>>>();
    let start = annotation
        .column("start")
        .unwrap()
        .u64()
        .unwrap()
        .into_iter()
        .flatten()
        .collect::<Vec<u64>>();
    let end = annotation
        .column("end")
        .unwrap()
        .u64()
        .unwrap()
        .into_iter()
        .flatten()
        .collect::<Vec<u64>>();

    let start_time = std::time::Instant::now();
    let mut bsxiter = bsx_file.into_iter();
    let total_n = annotation.len();

    let _ = itertools::izip!(&index, &start, &end)
        .enumerate()
        .for_each(|(n, (index_val, start_val, end_val))| {
            let df = bsxiter
                .get_region(index_val, start_val.clone(), end_val.clone())
                .expect("TODO: panic message");
            print!("{:?}/{:?}\r", n, total_n)
        });

    let duration = start_time.elapsed();
    println!("{:?}", duration)
}
