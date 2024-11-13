use crate::io::report::read::ReportReader;
use crate::io::report::types::ReportType;
use polars::io::ipc::BatchedWriter;
use polars::io::mmap::MmapBytesReader;
use polars::prelude::*;
use std::fs::File;
use std::io::Write;

pub struct BSXWriter<W: Write> {
    writer: BatchedWriter<W>,
    schema: Schema,
}

pub fn get_universal_schema(chr_names: Vec<String>) -> Schema {
    let categories = DataType::Enum(
        Some({
            let mut cat_builder =
                CategoricalChunkedBuilder::new(Default::default(), 0, Default::default());
            for chr_name in chr_names {
                let _ = &cat_builder.append(Some(chr_name.as_str()));
            }
            cat_builder
                .finish()
                .get_rev_map()
                .clone()
        }),
        CategoricalOrdering::Physical,
    );

    assert!(categories.contains_categoricals());

    let mut schema = Schema::default();
    schema.insert(PlSmallStr::from_str("chr"), categories);
    schema.insert(PlSmallStr::from_str("position"), DataType::UInt32);
    schema.insert(PlSmallStr::from_str("strand"), DataType::Boolean);
    schema.insert(PlSmallStr::from_str("context"), DataType::Boolean);
    schema.insert(PlSmallStr::from_str("count_m"), DataType::UInt16);
    schema.insert(PlSmallStr::from_str("count_total"), DataType::UInt16);
    schema.insert(PlSmallStr::from_str("density"), DataType::Float64);

    schema
}

impl<W: Write> BSXWriter<W> {
    pub fn new(writer: W, chr_names: Vec<String>) -> Self {
        let schema = get_universal_schema(chr_names);
        let writer = IpcWriter::new(writer)
            .with_compat_level(CompatLevel::newest())
            .with_compression(None)
            .batched(&schema)
            .unwrap();
        Self { writer, schema }
    }

    pub fn write_batch(&mut self, batch: &mut DataFrame) -> Result<(), PolarsError> {
        self.writer
            .write_batch(batch.align_chunks())
    }

    pub fn finish(mut self) -> Result<(), PolarsError> {
        self.writer.finish()
    }

    pub fn get_schema(&self) -> &Schema {
        &self.schema
    }
}


pub struct ConvertReportOptions {
    batch_per_read: Option<usize>,
    batch_size: Option<usize>,
    chunk_size: Option<usize>,
    low_memory: Option<bool>,
    n_threads: Option<usize>,
    rechunk: Option<bool>,
}

impl Default for ConvertReportOptions {
    fn default() -> Self {
        Self {
            batch_per_read: Some(usize::from(std::thread::available_parallelism().unwrap())),
            batch_size: Some(10_000),
            chunk_size: Some(10_000),
            low_memory: Some(false),
            n_threads: None,
            rechunk: Some(true),
        }
    }
}

impl ConvertReportOptions {
    pub fn with_low_memory(mut self, low_mem: bool) -> Self {
        self.low_memory = Some(low_mem);
        self
    }
    pub fn with_n_threads(mut self, n_threads: usize) -> Self {
        self.n_threads = Some(n_threads);
        self
    }
    pub fn with_batch_per_read(mut self, batch_per_read: usize) -> Self {
        self.batch_per_read = Some(batch_per_read);
        self
    }
    pub fn with_batch_size(mut self, batch_size: usize) -> Self {
        self.batch_size = Some(batch_size);
        self
    }
    pub fn with_chunk_size(mut self, chunk_size: usize) -> Self {
        self.chunk_size = Some(chunk_size);
        self
    }
    pub fn with_rechunk(mut self, rechunk: bool) -> Self {
        self.rechunk = Some(rechunk);
        self
    }
    pub fn finish(
        self,
        report_type: ReportType,
        report_path: &str,
        output_path: &str,
        fa_path: &str,
        fai_path: &str,
    ) -> Result<(), PolarsError> {
        let file = File::open(report_path).expect("could not open file");
        let csv_reader: CsvReader<Box<dyn MmapBytesReader>> = {
            report_type.get_reader(
                Box::new(file),
                self.chunk_size,
                self.low_memory,
                self.n_threads,
                self.rechunk,
            )
        };

        let mut reader = ReportReader::new(
            report_type,
            csv_reader,
            self.batch_per_read,
            self.batch_size,
            Some(fa_path.to_owned()),
            Some(fai_path.to_owned()),
        );

        let sink = File::create(output_path).expect("could not create output file");
        let mut writer = BSXWriter::new(sink, reader.get_chr_names().unwrap());
        let schema = writer.get_schema().clone();
        let schema_cast = PlHashMap::from_iter(
            schema
                .iter()
                .map(|(k, v)| (k.as_str(), v.clone())),
        );
        let schema_names = schema
            .iter_names()
            .map(|x| x.clone())
            .collect::<Vec<_>>();

        for mut batch in reader {
            batch = batch
                .lazy()
                .cast(schema_cast.clone(), true)
                .sort(
                    ["position"],
                    SortMultipleOptions::default()
                        .with_order_descending(false)
                        .with_multithreaded(true),
                )
                .collect()?
                .select(schema_names.clone())?;
            writer
                .write_batch(&mut batch)
                .expect("could not write batch");
        }
        writer
            .finish()
            .expect("could not finish write");
        Ok(())
    }
}
