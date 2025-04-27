pub mod bsx;
#[cfg(feature = "compression")]
pub mod compression;
pub mod report;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::fasta::fai::io::Reader as FaiReader;
use noodles::fasta::fs::index as index_fasta;
use noodles::fasta::fai::Record;

pub(crate) fn read_chrom<P: AsRef<Path>>(
    path: P,
    is_index: bool
) -> anyhow::Result<Vec<String>> {
    let index = if is_index {
        FaiReader::new(BufReader::new(File::open(path)?))
            .read_index()?
    } else {
        index_fasta(path)?
    };
    let records: Vec<Record> = index.into();
    Ok(records.into_iter().map(|r| String::from_utf8_lossy(r.name()).to_string()).collect())
}