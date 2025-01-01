use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;
use std::ops::Div;
use std::sync::{Arc, Mutex};
use std::sync::mpmc::{Receiver, Sender};
use std::thread;
use std::thread::JoinHandle;
use itertools::Itertools;
use log::info;
use polars::io::mmap::MmapBytesReader;
use polars::prelude::*;
use polars_arrow::array::Utf8ViewArray;
use crate::io::report::read::OwnedBatchedCsvReader;
use crate::io::report::types::ReportType;
use crate::region::RegionCoordinates;
use crate::ubatch2::{BSXBatch, BSXCol, UBatch};
use crate::utils::types::BSXResult;

/*
Reading plan
------------

1. Get the raw batch.
2. Convert it to BSX.
3. Align data to reference CONDITIONED
4. Cast chromosome to Enum

*/

struct NumberedBatch {
    n: usize,
    data: DataFrame
}

impl NumberedBatch {
    pub fn new(n: usize, data: DataFrame) -> Self {
        Self { n, data }
    }
}

impl Eq for NumberedBatch {}

impl PartialEq<Self> for NumberedBatch {
    fn eq(&self, other: &Self) -> bool {
        self.n == other.n
    }
}

impl PartialOrd<Self> for NumberedBatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for NumberedBatch {
    fn cmp(&self, other: &Self) -> Ordering {
        self.n.cmp(&other.n)
    }
}

struct ReportReader {
    receiver: Receiver<NumberedBatch>,
    fasta_reader: Arc<bio::io::fasta::IndexedReader<Box<dyn MmapBytesReader>>>,
    reader_handle: JoinHandle<()>,
    batch_cache: UBatch,
    chunk_size: usize,
    report_type: ReportType,
    do_align: bool,
    seq_cache: Option<Vec<u8>>,
    seq_region: RegionCoordinates
}



impl Iterator for ReportReader {
    type Item = UBatch;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.batch_cache.length() < self.chunk_size {
            if !self.receiver.is_empty() {
                let queue_len = self.receiver.iter().count();
                let mut heap = BinaryHeap::new();
                for i in 0..queue_len {
                    heap.push(Reverse(self.receiver.recv().unwrap()));
                }
                while let Some(Reverse(batch)) = heap.pop() {
                    let converted = self.convert_batch(batch.data);
                    self.batch_cache.extend(&converted);
                }
                self.next()
            } else {
                None
            }
        } else {
            let (out, back) = self.batch_cache.split(self.chunk_size as i64);
            self.batch_cache = back;
            Some(out)
        }
    }
}

impl ReportReader {
    fn convert_batch(&self, batch: DataFrame) -> UBatch {
        todo!()
    }
}

#[derive(Clone)]
/// All it does is reads batches of data from
/// report with [OwnedBatchedCsvReader] and forwards
/// it to [Sender]
struct ReaderThread {
    /// How many batches will be read per read
    batch_per_read: usize,
    batched_reader: Arc<Mutex<OwnedBatchedCsvReader>>,
    sender: Arc<Sender<NumberedBatch>>
}

impl ReaderThread {
    fn new(
        batch_per_read: usize,
        batched_reader: Arc<Mutex<OwnedBatchedCsvReader>>, 
        sender: Arc<Sender<NumberedBatch>>
    ) -> Self {
        Self { batched_reader, sender, batch_per_read}
    }
    
    fn run(&mut self) -> JoinHandle<()> {
        let cloned_self = self.clone();
        thread::spawn(
            move || {
                let mut reader = cloned_self.batched_reader.lock().unwrap();
    
                info!("ReaderThread started");
                let count: usize = 0;
                while let Some(batches) = reader.next_batches(cloned_self.batch_per_read).unwrap() {
                    for batch in batches {
                        let item = NumberedBatch::new(count, batch);
                        cloned_self.sender.send(item).unwrap()
                    }
                }
            } 
        )
    }
}

/// This method converts vector of a strings into a [DataType::Categorical]
pub fn get_enum_chr(chr_names: &Vec<String>) -> DataType {
    let categories = DataType::Enum(
        Some(Arc::new(RevMapping::build_local({
            Utf8ViewArray::from_slice(chr_names.iter().map(|v| Some(v)).collect_vec().as_slice())
        }))),
        CategoricalOrdering::Physical,
    );
    categories
}

pub fn cast_chr_to_categorical(data_frame: DataFrame, data_type: &DataType) -> BSXResult<DataFrame> {
    match data_frame.lazy().cast(PlHashMap::from_iter([(BSXCol::Chr.as_str(), data_type.clone())]), true).collect() {
        Ok(data) => Ok(data),
        Err(e) => Err(Box::new(e)),
    }
}

/// Converts strand from String to Boolean
fn _strand_expr() -> Expr {
    when(col("strand").eq(lit("+")))
        .then(true)
        .when(col("strand").eq(lit("-")))
        .then(false)
        .otherwise(lit(NULL))
        .cast(BSXCol::Strand.data_type())
}

/// Converts contexts from String to Boolean
fn _context_expr() -> Expr {
    when(col("context").eq(lit("CG")))
        .then(lit(true))
        .when(col("context").eq(lit("CHG")))
        .then(lit(false))
        .otherwise(lit(NULL))
        .cast(DataType::Boolean)
}
/// Converts nuc field from [ReportType::CGMAP] to Boolean 
/// strand
fn _nuc_expr() -> Expr {
    when(col("nuc").eq(lit("C")))
        .then(true)
        .when(col("nuc").eq(lit("G")))
        .then(false)
        .otherwise(lit(NULL))
        .cast(DataType::Boolean)
}

/// Unifies [DataType] of [ReportType] into [BSXBatch] columns
/// 
/// 1. Converts "context" column to [DataType::Boolean]
/// 2. Converts "strand" column to [DataType::Boolean]
/// 3. Casts other columns to [BSXCol]
/// 
/// # Important
/// This method does not convert "chr" column to [DataType::Categorical]
pub fn report_type_to_bsx(report_type: &ReportType, df: DataFrame) -> BSXResult<DataFrame> {
    let lazy_frame = df.lazy();

    let res = match report_type {
        ReportType::BISMARK => lazy_frame
            // count_total
            .with_column(
                (col("count_m") + col("count_um"))
                    .cast(BSXCol::CountTotal.data_type())
                    .alias(BSXCol::CountTotal.as_str()),
            )
            // density
            .with_column(
                (col("count_m")
                    .cast(BSXCol::Density.data_type())
                    .div(BSXCol::CountTotal.into()))
                    .alias(BSXCol::Density.as_str()),
            )
            .with_columns([
                _strand_expr().alias(BSXCol::Strand.as_str()),
                _context_expr().alias(BSXCol::Context.as_str()),
            ]),
        ReportType::CGMAP => lazy_frame
            .with_columns([
                _nuc_expr().alias("strand"),
                _context_expr().alias("context"),
            ])
            .with_column(
                (col("count_m")
                    .cast(BSXCol::Density.data_type())
                    .div(BSXCol::CountTotal.into()))
                    .alias(BSXCol::Density.as_str())
            ),
        ReportType::BEDGRAPH => lazy_frame
            .rename(["start"], ["position"], true)
            .drop(["end"])
            .with_columns([
                lit(NULL).alias("strand").cast(BSXCol::Strand.data_type()),
                lit(NULL).alias("context").cast(BSXCol::Context.data_type()),
                lit(NULL).alias("count_m").cast(BSXCol::CountM.data_type()),
                lit(NULL).alias("count_total").cast(BSXCol::CountTotal.data_type()),
                col("density").div(lit(100)).alias("density").cast(BSXCol::Density.data_type()),
            ]),
        ReportType::COVERAGE => lazy_frame
            .rename(["start"], ["position"], true)
            .drop(["end"])
            .with_column((col("count_m") + col("count_um")).alias("count_total"))
            .with_columns([
                lit(NULL).alias("strand"),
                lit(NULL).alias("context"),
                col("density").div(lit(100)).alias("density"),
            ]),
    }.cast(
        PlHashMap::from_iter(vec![
            (BSXCol::Position.as_str(), BSXCol::Position.data_type()),
            (BSXCol::Strand.as_str(), BSXCol::Strand.data_type()),
            (BSXCol::Context.as_str(), BSXCol::Context.data_type()),
            (BSXCol::Density.as_str(), BSXCol::Density.data_type()),
            (BSXCol::CountTotal.as_str(), BSXCol::CountTotal.data_type()),
            (BSXCol::CountM.as_str(), BSXCol::CountM.data_type()),
        ]), false
    )
        .select(BSXCol::all_cols().into_iter().map(col).collect_vec())
        .collect()?;

    Ok(res)
}


#[cfg(test)]
mod tests {
    use super::*;
    
    fn get_dummy_bismark() -> DataFrame {
        df!(
            "chr" => &["chr1", "chr1", "chr2", "chr2", "chr2"],
            "position" => &[1, 2, 6, 7, 8],
            "strand" => &["+", "-", "-", "+", "+"],
            "count_m" => &[10, 20, 30, 0, 0],
            "count_um" => &[0, 0, 0, 40, 50],
            "context" => &["CG", "CHG", "CHH", "CG", "CHH"],
            "trinuc" => &["CGA", "CAG", "CAT", "CGT", "CTT"]
        ).unwrap()
    }

    fn cast_to_bsx_str(df: DataFrame) -> DataFrame {
        df.lazy().cast(
            PlHashMap::from_iter([
                (BSXCol::Chr.as_str(), DataType::String),
                (BSXCol::Strand.as_str(), BSXCol::Strand.data_type()),
                (BSXCol::Position.as_str(), BSXCol::Position.data_type()),
                (BSXCol::Context.as_str(), BSXCol::Context.data_type()),
                (BSXCol::CountM.as_str(), BSXCol::CountM.data_type()),
                (BSXCol::CountTotal.as_str(), BSXCol::CountTotal.data_type()),
                (BSXCol::Density.as_str(), BSXCol::Density.data_type()),
            ]), true
        ).collect().unwrap()
    }

    fn get_dummy_cgmap() -> DataFrame {
        df!(
            "chr" => &["chr1", "chr1", "chr2", "chr2", "chr2"],
            "nuc" => &["C", "G", "G", "C", "C"],
            "position" => &[1, 2, 6, 7, 8],
            "context" => &["CG", "CHG", "CHH", "CG", "CHH"],
            "dinuc" => &["GA", "AG", "AT", "GT", "TT"],
            "count_m" => &[10, 20, 30, 0, 0],
            "count_total" => &[10, 20, 30, 40, 50],
        ).unwrap()
    }

    fn get_dummy_coverage() -> DataFrame {
        df!(
            "chr" => &["chr1", "chr1", "chr2", "chr2", "chr2"],
            "start" => &[1, 2, 6, 7, 8],
            "end" => &[1, 2, 6, 7, 8],
            "density" => &[100., 100., 100., 0.0, 0.0],
            "count_m" => &[10, 20, 30, 0, 0],
            "count_um" => &[0, 0, 0, 40, 50],
        ).unwrap()
    }

    fn get_dummy_bedgraph() -> DataFrame {
        df!(
            "chr" => &["chr1", "chr1", "chr2", "chr2", "chr2"],
            "start" => &[1, 2, 6, 7, 8],
            "end" => &[1, 2, 6, 7, 8],
            "density" => &[100., 100., 100., 0.0, 0.0],
        ).unwrap()
    }
    
    #[test]
    fn check_enum_dtype_creation() {
        let dummy_chr_names: Vec<String> = vec!["chr1".to_string(), "chr2".to_string()];
        let data_type = get_enum_chr(&dummy_chr_names);
        assert!(data_type.is_enum());
        let dummy_series = Series::from_iter(dummy_chr_names);
        let dummy_cat = dummy_series.clone().cast(&data_type);
        assert!(dummy_cat.is_ok());
        let dummy_cat = dummy_cat.unwrap();
        assert_eq!(dummy_cat.null_count(), 0);
    }
    
    #[test]
    fn test_report_type_to_bsx() {
        let correct = df!(
            "chr" => &["chr1", "chr1", "chr2", "chr2", "chr2"],
            "strand" => &[Some(true), Some(false), Some(false), Some(true), Some(true)],
            "position" => &[1, 2, 6, 7, 8],
            "context" => &[Some(true), Some(false), None, Some(true), None],
            "count_m" => &[10, 20, 30, 0, 0],
            "count_total" => &[10, 20, 30, 40, 50],
            "density" => &[1.0, 1.0, 1.0, 0.0, 0.0],
        ).unwrap();
        let correct = cast_to_bsx_str(correct);

        let dummy_df = get_dummy_bismark();
        let converted = report_type_to_bsx(&ReportType::BISMARK, dummy_df);
        assert!(converted.is_ok(), "{:?}", converted);
        assert_eq!(correct, converted.unwrap());

        let dummy_df = get_dummy_cgmap();
        let converted = report_type_to_bsx(&ReportType::CGMAP, dummy_df);
        assert!(converted.is_ok(), "{:?}", converted);
        assert_eq!(correct, converted.unwrap());

        let correct = df!(
            "chr" => &["chr1", "chr1", "chr2", "chr2", "chr2"],
            "strand" => &[None::<bool>, None, None, None, None],
            "position" => &[1, 2, 6, 7, 8],
            "context" => &[None::<bool>, None, None, None, None],
            "count_m" => &[10, 20, 30, 0, 0],
            "count_total" => &[10, 20, 30, 40, 50],
            "density" => &[1.0, 1.0, 1.0, 0.0, 0.0],
        ).unwrap();
        let correct = cast_to_bsx_str(correct);
        let dummy_df = get_dummy_coverage();
        let converted = report_type_to_bsx(&ReportType::COVERAGE, dummy_df);
        assert!(converted.is_ok(), "{:?}", converted);
        assert_eq!(correct, converted.unwrap());

        let correct = df!(
            "chr" => &["chr1", "chr1", "chr2", "chr2", "chr2"],
            "strand" => &[None::<bool>, None, None, None, None],
            "position" => &[1, 2, 6, 7, 8],
            "context" => &[None::<bool>, None, None, None, None],
            "count_m" => &[None::<u32>, None, None, None, None],
            "count_total" => &[None::<u32>, None, None, None, None],
            "density" => &[1.0, 1.0, 1.0, 0.0, 0.0],
        ).unwrap();
        let correct = cast_to_bsx_str(correct);
        let dummy_df = get_dummy_bedgraph();
        let converted = report_type_to_bsx(&ReportType::BEDGRAPH, dummy_df);
        assert!(converted.is_ok(), "{:?}", converted);
        assert_eq!(correct, converted.unwrap());
    }
    
    #[test]
    fn test_categorical_cast() {
        let dummy_df = get_dummy_bismark();
        let dummy_df = report_type_to_bsx(&ReportType::BISMARK, dummy_df).unwrap();
        let data_type = get_enum_chr(&["chr1", "chr2"].into_iter().map(String::from).collect_vec());
        
        let new = cast_chr_to_categorical(dummy_df, &data_type);
        assert!(new.is_ok(), "{:?}", new);
        let new = new.unwrap();
        assert_eq!(new.column("chr").unwrap().dtype(), &data_type);
    }
}