#![feature(btree_cursors)]

mod genome;
mod ipc_reader;
use ipc_reader::FileReader as IPCFileReader;
use crate::genome::Annotation;
use polars::export::arrow::datatypes::IntegerType;
use polars::io::ipc::BatchedWriter;
use polars::io::ArrowReader;
use polars::prelude::{as_struct, col, lit, when, ArgAgg, BatchedCsvReader, CsvParseOptions, CsvReadOptions, CsvReader, CsvWriter, Expr, GetOutput, IntoLazy, IpcReader, IpcWriterOptions, SerWriter, NULL};
use polars_core::prelude::*;
use polars_core::utils::arrow::io::ipc::read::{read_file_metadata, FileReader};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::ops::{BitOr, Not};
use std::ops::Bound::Excluded;
use std::time::Instant;
use polars::export::arrow::io::ipc::read::FileMetadata;

/// Get schema of BSXplorer internal format
fn get_universal_schema(chr_names: Vec<&str>) -> Schema {
    let categories = DataType::Enum(
        Some({
            let mut cat_builder =
                CategoricalChunkedBuilder::new(Default::default(), 0, Default::default());
            for chr_name in chr_names {
                let _ = &cat_builder.append(Some(chr_name));
            }
            cat_builder.finish().get_rev_map().clone()
        }),
        CategoricalOrdering::Physical,
    );

    assert!(categories.contains_categoricals());

    let mut schema = Schema::default();
    schema.insert(PlSmallStr::from_str("chr"), categories);
    schema.insert(PlSmallStr::from_str("position"), DataType::UInt64);
    schema.insert(PlSmallStr::from_str("strand"), DataType::Boolean);
    schema.insert(PlSmallStr::from_str("context"), DataType::Boolean);
    schema.insert(PlSmallStr::from_str("count_m"), DataType::UInt16);
    schema.insert(PlSmallStr::from_str("count_total"), DataType::UInt16);
    schema.insert(PlSmallStr::from_str("density"), DataType::Float32);

    schema
}

/// Struct for sequence data
struct SeqData {
    /// Names of the records from FASTA
    chr_names: Vec<String>,
    /// Binary offsets of chroms
    chr_offsets: Vec<u64>,
    /// Number of G/C-nucleotides
    gc_content: u64,
    /// Number of non-G/C nucleotides
    non_gc_content: u64,
}

impl SeqData {
    fn new() -> Self {
        SeqData {
            chr_names: Vec::new(),
            chr_offsets: Vec::new(),
            gc_content: 0,
            non_gc_content: 0,
        }
    }
}

/// Scans fasta to determine
///
/// - Number of nucleotides
/// - G/C content count
/// - Chromosome names
/// - Sequence binary offsets
fn scan_fasta(path: &str) -> SeqData {
    let mut reader: BufReader<File> = BufReader::new(File::open(path).unwrap());
    let mut buf = [0_u8; 1];

    let mut reading_name = false;
    let mut reading_tag = false;
    let mut name = String::new();
    let mut seq_offset: u64;

    let mut seq_data = SeqData::new();

    while reader.read_exact(&mut buf).is_ok() {
        let c = u8::from_le_bytes(buf);

        // Start of tag
        if c == 0x3e {
            reading_name = true
        };
        // End of name
        if reading_name && c == 0x20 {
            reading_tag = true;
            reading_name = false;
            seq_data.chr_names.push(name.clone());
            name.clear()
        }
        // End of tag
        if reading_tag && c == 0x0a {
            reading_tag = false;
            // Save seq offset
            seq_offset = reader
                .stream_position()
                .expect("Could not get current position!");
            seq_data.chr_offsets.push(seq_offset);
        }

        // Push name char (if name continues)
        if reading_name && !reading_tag {
            if c != 0x3e {
                name.push(c as char)
            }
        }
        // Add GC content
        else if !reading_name && !reading_tag {
            // If not newline
            match c {
                // T|t|G|g
                0x54 | 0x74 | 0x47 | 0x67 => seq_data.gc_content += 1,
                // \n
                0x0a => {}
                // Other
                _ => seq_data.non_gc_content += 1,
            }
        }
    }
    seq_data
}

/// Type of methylation report
enum ReportType {
    BEDGRAPH,
    COVERAGE,
    CGMAP,
    BISMARK,
}

impl ReportType {
    /// Get corresponding polars schema for report type
    fn schema(&self) -> Schema {
        match self {
            ReportType::BISMARK => {
                let mut schema = Schema::default();
                schema.insert(PlSmallStr::from_str("chr"), DataType::String);
                schema.insert(PlSmallStr::from_str("position"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("strand"), DataType::String);
                schema.insert(PlSmallStr::from_str("count_m"), DataType::UInt32);
                schema.insert(PlSmallStr::from_str("count_um"), DataType::UInt32);
                schema.insert(PlSmallStr::from_str("context"), DataType::String);
                schema.insert(PlSmallStr::from_str("trinuc"), DataType::String);
                schema
            }
            ReportType::COVERAGE => {
                let mut schema = Schema::default();
                schema.insert(PlSmallStr::from_str("chr"), DataType::String);
                schema.insert(PlSmallStr::from_str("start"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("end"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("density"), DataType::Float32);
                schema.insert(PlSmallStr::from_str("count_m"), DataType::UInt32);
                schema.insert(PlSmallStr::from_str("count_um"), DataType::UInt32);
                schema
            }
            ReportType::CGMAP => {
                let mut schema = Schema::default();
                schema.insert(PlSmallStr::from_str("chr"), DataType::String);
                schema.insert(PlSmallStr::from_str("nuc"), DataType::String);
                schema.insert(PlSmallStr::from_str("position"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("context"), DataType::String);
                schema.insert(PlSmallStr::from_str("dinuc"), DataType::String);
                schema.insert(PlSmallStr::from_str("count_m"), DataType::UInt32);
                schema.insert(PlSmallStr::from_str("count_total"), DataType::UInt32);
                schema
            }
            ReportType::BEDGRAPH => {
                let mut schema = Schema::default();
                schema.insert(PlSmallStr::from_str("chr"), DataType::String);
                schema.insert(PlSmallStr::from_str("start"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("end"), DataType::UInt64);
                schema.insert(PlSmallStr::from_str("density"), DataType::Float32);
                schema
            }
        }
    }
    /// Get read options for report type
    fn read_options(
        &self,
        chunk_size: Option<usize>,
        low_memory: Option<bool>,
        n_threads: Option<usize>,
    ) -> CsvReadOptions {
        match self {
            ReportType::BISMARK
            | ReportType::COVERAGE
            | ReportType::CGMAP
            | ReportType::BEDGRAPH => CsvReadOptions::default()
                .with_has_header(false)
                .with_chunk_size(chunk_size.unwrap_or(2 << 20))
                .with_low_memory(low_memory.unwrap_or(false))
                .with_n_threads(n_threads)
                .with_schema(Some(SchemaRef::new(self.schema())))
                .with_parse_options({
                    CsvParseOptions::default()
                        .with_separator(b'\t')
                        .with_try_parse_dates(false)
                        .with_quote_char(Some(b'#'))
                }),
        }
    }

    /// Get mutation function to convert report format to BSX
    fn mutate_func(&self, data_frame: DataFrame, schema: &Schema) -> PolarsResult<DataFrame> {
        let lazy_frame = data_frame.lazy();
        let context_expr: Expr = when(col("context").eq(lit("CG")))
            .then(lit(true))
            .when(col("context").eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean);
        match self {
            ReportType::BEDGRAPH => todo!(),
            ReportType::COVERAGE => todo!(),
            ReportType::CGMAP => lazy_frame.with_columns([col("nuc").eq(lit("C")).alias("strand")]),
            ReportType::BISMARK => lazy_frame
                .with_column((col("count_m") + col("count_um")).alias("count_total"))
                .with_columns([
                    (col("count_m") / col("count_total")).alias("density"),
                    col("strand").eq(lit("+")).alias("strand"),
                ]),
        }
        .select([
            col("chr")
                .cast(schema.get_field("chr").unwrap().dtype)
                .cast(DataType::from_arrow(
                    &ArrowDataType::Dictionary(
                        IntegerType::UInt32,
                        Box::new(ArrowDataType::Utf8),
                        false,
                    ),
                    false,
                )),
            col("position").cast(schema.get_field("position").unwrap().dtype),
            col("strand").cast(schema.get_field("strand").unwrap().dtype),
            context_expr
                .alias("context")
                .cast(schema.get_field("context").unwrap().dtype),
            col("count_m").cast(schema.get_field("count_m").unwrap().dtype),
            col("count_total").cast(schema.get_field("count_total").unwrap().dtype),
            col("density").cast(schema.get_field("density").unwrap().dtype),
        ])
        .collect()
    }

    /// Master function to convert report into BSX schema
    fn convert_report(
        &self,
        input_path: &str,
        output_path: &str,
        fasta_path: &str,
        chunk_size: Option<usize>,
        low_memory: Option<bool>,
        n_threads: Option<usize>,
    ) {
        // Init chromosome names
        let seq_data: SeqData = scan_fasta(fasta_path);
        let chr_names: Vec<&str> = seq_data
            .chr_names
            .iter()
            .map(|s| &s as &str)
            .collect::<Vec<&str>>();

        // Init reader
        let mut binding: CsvReader<File> = self
            .read_options(chunk_size, low_memory, n_threads)
            .try_into_reader_with_file_path(Some(input_path.parse().unwrap()))
            .expect("could not load CSV");
        let mut reader: BatchedCsvReader =
            binding.batched_borrowed().expect("could not borrow CSV");

        // Init schema
        let schema: Schema = get_universal_schema(chr_names);

        // Init writer with categorical chromosome names
        let mut writer: BatchedWriter<File> = IpcWriterOptions::default()
            .to_writer(File::create(output_path).expect("could not create output file"))
            .batched(&schema)
            .expect("could not create ipc writer");

        // Read
        while let Some(batch_vec) = reader.next_batches(10).unwrap() {
            for batch in batch_vec {
                self.mutate_func(batch, &schema)
                    .expect("could not mutate record function")
                    .partition_by(["chr", "context", "strand"], true)
                    .expect("could not partition batch")
                    .iter_mut()
                    .try_for_each(|x| writer.write_batch(x.align_chunks()))
                    .expect("could not write batch");
            }
        }
        writer.finish().expect("could not finish write");
    }
}

#[derive(Eq, Hash, PartialEq, Copy, Clone)]
enum Context {
    CG,
    CHG,
    CHH,
    ALL,
}

impl Context {
    fn from_bool(value: Option<bool>) -> Self {
        match value {
            Some(true) => Context::CG,
            Some(false) => Context::CHG,
            None => Context::CHH,
        }
    }
}

#[derive(Eq, Hash, PartialEq, Copy, Clone)]
enum Strand {
    Forward,
    Reverse,
    None,
}

impl Strand {
    fn from_bool(value: Option<bool>) -> Self {
        match value {
            Some(true) => Strand::Forward,
            Some(false) => Strand::Reverse,
            None => Strand::None,
        }
    }

    fn from_str(value: &str) -> Self {
        match value {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::None,
        }
    }
}

type IpcIndex = HashMap<(String, Strand, Context), BTreeMap<u64, usize>>;
fn get_index(path: &str) -> Result<IpcIndex, PolarsError> {
    let mut file = File::open(path).expect("could not open ipc path");
    let metadata = read_file_metadata(&mut file).expect("could not read ipc metadata");
    let mut reader = FileReader::new(file, metadata.clone(), None, None);

    let mut num_batch: usize = 0;
    let mut index: HashMap<(String, Strand, Context), BTreeMap<u64, usize>> = HashMap::new();

    while let Some(batch) = reader.next_record_batch()? {
        let df = DataFrame::try_from((batch, metadata.schema.clone().as_ref()))
            .expect("could not create dataframe");

        let pos = df
            .column("position")?
            .u64()?
            .get(0)
            .expect("could not get position");
        let key: (String, Strand, Context) = {
            let chr_col = df.column("chr")?.categorical()?;
            (
                chr_col
                    .get_rev_map()
                    .get(chr_col.physical().get(0).unwrap())
                    .to_owned(),
                Strand::from_bool(df.column("strand")?.bool()?.get(0)),
                Context::from_bool(df.column("context")?.bool()?.get(0)),
            )
        };

        match index.get_mut(&key) {
            Some(btree) => {
                btree.insert(pos, num_batch);
            }
            None => {
                let mut btree: BTreeMap<u64, usize> = BTreeMap::new();
                btree.insert(pos, num_batch);
                index.insert(key, btree);
            }
        };
        num_batch += 1;
    }
    Ok(index)
}

fn add_index(annotation: Annotation, index: IpcIndex, context: Option<Context>) -> DataFrame {
    let df = annotation
        .finish()
        .expect("could not make annotation DataFrame");
    let context = context.unwrap_or(Context::CG);
    df.lazy()
        .with_column(
            as_struct(vec![col("chr"), col("strand"), col("start"), col("end")])
                .apply(
                    move |ser| series_to_index(ser, context, &index),
                    GetOutput::from_type(DataType::List(Box::new(DataType::UInt32))),
                )
                .alias("index"),
        )
        .collect()
        .expect("could not make index col")
}

fn series_to_index(
    series: Series,
    context: Context,
    index_map: &IpcIndex,
) -> Result<Option<Series>, PolarsError> {
    let s = series.struct_()?;
    // Bindings
    let chr_col = &s.field_by_name("chr")?;
    let strand_col = &s.field_by_name("strand")?;
    let start_col = &s.field_by_name("start")?;
    let end_col = &s.field_by_name("end")?;
    // Series
    let chr_arr = chr_col.str()?;
    let strand_arr = strand_col.str()?;
    let start_arr = start_col.u64()?;
    let end_arr = end_col.u64()?;

    let out: ChunkedArray<ListType> = itertools::izip!(
        chr_arr.into_iter(),
        strand_arr.into_iter(),
        start_arr.into_iter(),
        end_arr.into_iter()
    )
    .map(|(opt_chr, opt_strand, opt_start, opt_end)| {
        find_indices(opt_chr, opt_strand, opt_start, opt_end, context, &index_map)
    })
    .collect();
    Ok(Some(Series::from(out)))
}

fn find_indices(
    opt_chr: Option<&str>,
    opt_strand: Option<&str>,
    opt_start: Option<u64>,
    opt_end: Option<u64>,
    context: Context,
    index_map: &IpcIndex,
) -> Option<Series> {
    match (opt_chr, opt_strand, opt_start, opt_end) {
        (Some(chr), Some(strand), Some(start), Some(end)) => {
            let btree = index_map.get(&(chr.to_string(), Strand::from_str(strand), context));
            if let Some(btree) = btree {
                let mut root = btree.lower_bound(Excluded(&start));
                let mut batches = vec![root.peek_prev()];
                while let Some(new_batch) = root.next() {
                    if *new_batch.0 <= end {
                        batches.push(Some(new_batch));
                    } else {
                        break;
                    }
                }
                let idxs = batches
                    .iter()
                    .flatten()
                    .map(|x| *x.1 as u32)
                    .collect::<Vec<u32>>();

                let arr = UInt32Chunked::from_vec("index".into(), idxs);
                Some(arr.into_series())
            } else {
                None
            }
        }
        _ => None,
    }
}

fn series_to_offsets(series: Series, mut ipc_reader: IPCFileReader<File>) -> Result<Option<Series>, PolarsError> {
    let s = series.struct_()?;
    // Bindings
    let start_col = &s.field_by_name("start")?;
    let end_col = &s.field_by_name("end")?;
    let index_col = &s.field_by_name("index")?;
    // Series
    let start_arr = start_col.u64()?;
    let end_arr = end_col.u64()?;
    let index_arr = index_col.u32()?;
    
    let metadata = ipc_reader.metadata().clone();

    let offsets: ChunkedArray<ListType> = itertools::izip!(start_arr.into_iter(), end_arr.into_iter(), index_arr.into_iter())
        .map(|(opt_start, opt_end, opt_index)| 
            find_offsets(&mut ipc_reader, &metadata, opt_start, opt_end, opt_index)
        ).collect();
    Ok(Some(Series::from(offsets)))
}

fn find_offsets(
    ipc_reader: &mut IPCFileReader<File>,
    metadata: &FileMetadata,
    opt_start: Option<u64>,
    opt_end: Option<u64>,
    opt_index: Option<u32>,
) -> Option<Series> {
    match (opt_start, opt_end, opt_index) {
        (Some(start), Some(end), Some(index)) => {
            if let Some(block) = ipc_reader.read_at(index as usize) {
                let block = block.expect("could not read block");
                let df = DataFrame::try_from((block, metadata.schema.as_ref()))
                    .expect("could not create dataframe");

                let pos = df.column("position").unwrap();
                let bin_range = pos.lt(start).unwrap()
                    .bitor(pos.gt(end).unwrap())
                    .not();
                let (mut start_pos, mut end_pos) = (None, None);
                for (idx, val) in bin_range.into_iter().enumerate() {
                    if val.is_none() { continue; };
                    match val.unwrap() {
                        true => if start_pos.is_none() { start_pos = Some(idx); }
                        false => if start_pos.is_some() && end_pos.is_none() { end_pos = Some(idx); break; }
                    }
                }

                let offsets = vec![
                    start_pos.unwrap_or(0) as u32,
                    end_pos.unwrap_or(bin_range.len()) as u32,
                ];
                let arr = UInt32Chunked::from_vec("offset".into(), offsets);
                Some(arr.into_series())
            } else {
                None
            }
        }
        _ => None
    }
}

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
    let df_w_index = add_index(
        annot.filter(col("type").eq(lit("gene"))),
        index_map,
        Some(Context::CG),
    );
    let duration = start.elapsed();

    let ipc_path = "/Users/shitohana/Desktop/RustProjects/bsxplorer2/res.ipc";
    let mut file = File::open("/Users/shitohana/Desktop/RustProjects/bsxplorer2/res.ipc").expect("could not open ipc path");
    let metadata = read_file_metadata(&mut file).expect("could not read ipc metadata");


    let df_w_offsets = df_w_index
        .lazy()
        .explode(["index"])
        .with_column(
        as_struct(vec![col("start"), col("end"), col("index")])
            .apply(move |series: Series| {
                let s = series.struct_()?;
                let mut file = File::open(ipc_path).expect("could not open ipc path");
                let metadata = read_file_metadata(&mut file).expect("could not read ipc metadata");
                let mut ipc_reader = IPCFileReader::new(file, metadata.clone(), None, None);

                // Bindings
                let start_col = &s.field_by_name("start")?;
                let end_col = &s.field_by_name("end")?;
                let index_col = &s.field_by_name("index")?;
                // Series
                let start_arr = start_col.u64()?;
                let end_arr = end_col.u64()?;
                let index_arr = index_col.u32()?;

                let offsets: ChunkedArray<ListType> = itertools::izip!(start_arr.into_iter(), end_arr.into_iter(), index_arr.into_iter())
                        .map(|(opt_start, opt_end, opt_index)|
                            match (opt_start, opt_end, opt_index) {
                            (Some(start), Some(end), Some(index)) => {
                                if let Some(block) = ipc_reader.read_at(index as usize) {
                                    let block = block.expect("could not read block");
                                    let df = DataFrame::try_from((block, metadata.schema.as_ref()))
                                        .expect("could not create dataframe");

                                    let pos = df.column("position").unwrap();
                                    let bin_range = pos.lt(start).unwrap()
                                        .bitor(pos.gt(end).unwrap())
                                        .not();
                                    let (mut start_pos, mut end_pos) = (None, None);
                                    for (idx, val) in bin_range.into_iter().enumerate() {
                                        if val.is_none() { continue; };
                                        match val.unwrap() {
                                            true => if start_pos.is_none() { start_pos = Some(idx); }
                                            false => if start_pos.is_some() && end_pos.is_none() { end_pos = Some(idx); break; }
                                        }
                                    }

                                    let offsets = vec![
                                        start_pos.unwrap_or(0) as u32,
                                        end_pos.unwrap_or(bin_range.len()) as u32,
                                    ];
                                    let arr = UInt32Chunked::from_vec("offset".into(), offsets);
                                    Some(arr.into_series())
                                } else {
                                    None
                                }
                            }
                            _ => None
                        }).collect();
                Ok(Some(Series::from(offsets)))

                },
                GetOutput::from_type(DataType::Array(Box::new(DataType::UInt32), 2))
            ).alias("offset")
        ).collect().unwrap();

    // println!("{:?}", df_w_offsets);
    let test_csv = File::create("/Users/shitohana/Desktop/RustProjects/bsxplorer2/test.csv").unwrap();
    CsvWriter::new(test_csv).with_separator(b'\t').finish(&mut df_w_offsets.lazy().with_column(col("offset").cast(DataType::List(Box::new(DataType::String))).list().join(lit(", "), true)).select([col("chr"), col("start"), col("end"), col("index"), col("offset")]).collect().unwrap()).unwrap();
    println!("{:?}", duration);
}
