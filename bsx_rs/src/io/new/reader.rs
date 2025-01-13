use crate::io::new::data_schema::ReportTypeSchema;
use polars::io::RowIndex;
use polars::prelude::*;
use std::error::Error;
use std::path::PathBuf;

pub struct ReportReaderBuilder {
    report_type: ReportTypeSchema,
    rechunk: bool,
    n_threads: Option<usize>,
    low_memory: bool,
    n_rows: Option<usize>,
    row_index: Option<RowIndex>,
    chunk_size: usize,
    skip_rows_after_header: usize,
    fasta_path: Option<PathBuf>,
    fai_path: Option<PathBuf>,
    batch_size: usize,
}

impl ReportReaderBuilder {
    pub fn with_rechunk(mut self, rechunk: bool) -> Self {
        self.rechunk = rechunk;
        self
    }
    pub fn with_n_threads(mut self, n_threads: usize) -> Self {
        self.n_threads = Some(n_threads);
        self
    }
    pub fn with_skip_rows_after_header(mut self, skip_rows_after_header: usize) -> Self {
        self.skip_rows_after_header = skip_rows_after_header;
        self
    }
    pub fn with_low_memory(mut self, low_memory: bool) -> Self {
        self.low_memory = low_memory;
        self
    }
    pub fn with_n_rows(mut self, n_rows: usize) -> Self {
        self.n_rows = Some(n_rows);
        self
    }
    pub fn with_report_type(mut self, report_type: ReportTypeSchema) -> Self {
        self.report_type = report_type;
        self
    }
    pub fn with_row_index(mut self, row_index: RowIndex) -> Self {
        self.row_index = Some(row_index);
        self
    }
    pub fn with_chunk_size(mut self, chunk_size: usize) -> Self {
        self.chunk_size = chunk_size;
        self
    }
    pub fn with_batch_size(mut self, batch_size: usize) -> Self {
        self.batch_size = batch_size;
        self
    }

    fn check_fasta_path(&self) -> Result<(), Box<dyn Error>> {
        if self.fasta_path.is_some() && self.fai_path.is_some() {
            Ok(())
        } else {
            Err(Box::from(
                "Both FASTA and FASTA index \
            paths should be specified",
            ))
        }
    }
    fn check_type_and_fasta(&self) -> Result<bool, Box<dyn Error>> {
        if self.report_type.need_align() || self.fasta_path.is_some() {
            match self.check_fasta_path() {
                Ok(_) => Ok(true),
                Err(e) => Err(e),
            }
        } else {
            Ok(false)
        }
    }
    //
    //     pub fn finish<F>(self, handle: F) -> ReportReader
    //     where
    //         F: Read + Seek + MmapBytesReader + 'static,
    //     {
    //         let reader = self.report_type.read_options()
    //             .with_low_memory(self.low_memory)
    //             .with_n_rows(self.n_rows)
    //             .with_skip_rows_after_header(self.skip_rows_after_header)
    //             .with_row_index(self.row_index.clone())
    //             .with_chunk_size(self.chunk_size)
    //             .with_low_memory(self.low_memory)
    //             .with_n_threads(self.n_threads)
    //             .into_reader_with_file_handle(handle);
    //
    //         let fasta_reader = match self.check_type_and_fasta() {
    //             Ok(true) => {
    //                 let fasta_reader = File::open(&self.fasta_path.unwrap())
    //                     .expect("Could not open FASTA file");
    //                 let index_reader = File::open(&self.fai_path.unwrap())
    //                     .expect("Could not open FASTA index file");
    //                 let indexed_reader = IndexedReader::new(fasta_reader, index_reader)
    //                     .expect("Could not create indexed reader");
    //                 Some(indexed_reader)
    //             }
    //             Ok(false) => None,
    //             Err(e) => panic!("Could not read FASTA file: {:?}", e),
    //         };
    //
    //         ReportReader::new(
    //             reader,
    //             fasta_reader,
    //             self.report_type,
    //             self.batch_size,
    //             self.chunk_size,
    //         )
    //     }
}
//
// pub(crate) fn get_reader<F>(handle: F, report_type: ReportTypeSchema) -> CsvReader<F>
// where
//     F: Read + Seek + MmapBytesReader,
// {
//     let options: CsvReadOptions = report_type.read_options();
//     options.into_reader_with_file_handle(handle)
// }
//
// pub struct OwnedBatchedCsvReader<F>
// where F: Read + Seek + MmapBytesReader{
//     #[allow(dead_code)]
//     // this exists because we need to keep ownership
//     pub schema: SchemaRef,
//     pub batched_reader: BatchedCsvReader<'static>,
//     // keep ownership
//     pub _reader: CsvReader<F>,
// }
//
// impl<F> OwnedBatchedCsvReader<F>
// where F: Read + Seek + MmapBytesReader + 'static {
//     pub(crate) fn new(mut reader: CsvReader<F>, schema: Arc<Schema>) -> Self {
//         let batched_reader = reader.batched_borrowed()
//             .expect("Could not create batched CSV reader.");
//         let batched_reader: BatchedCsvReader<'static> =
//             unsafe { std::mem::transmute(batched_reader) };
//         Self {
//             batched_reader,
//             _reader: reader,
//             schema,
//         }
//     }
// }
//
// impl<F> OwnedBatchedCsvReader<F>
// where F: Read + Seek + MmapBytesReader + 'static {
//     pub fn next_batches(&mut self, n: usize) -> PolarsResult<Option<Vec<DataFrame>>> {
//         self.batched_reader.next_batches(n)
//     }
// }
//
// pub struct ReportReader
// {
//     // Reader thread
//     join_handle: JoinHandle<()>,
//     receiver: Receiver<DataFrame>,
//     report_schema: ReportSchema,
//     // Fasta
//     fasta_reader: Option<IndexedReader<File>>,
//
//     // Config
//     chunk_size: usize,
// }
//
// impl ReportReader
// {
//     pub(crate) fn new<F>(
//         reader: CsvReader<F>,
//         fasta_reader: Option<IndexedReader<File>>,
//         report_schema: ReportSchema,
//         batch_per_read: usize,
//         chunk_size: usize
//     ) -> Self
//     where F: MmapBytesReader + 'static {
//         let (sc, receiver) = mpsc::channel();
//
//         let join_handle = thread::spawn(move || {
//             let mut owned_batched = OwnedBatchedCsvReader::new(
//                 reader,
//                 Arc::from(report_schema.get().schema()),
//             );
//             while let Ok(Some(data)) = owned_batched.next_batches(batch_per_read) {
//                 data.into_iter()
//                     .for_each(|frame| {sc.send(frame).expect("Could not send data to main thread.")});
//             }
//         });
//
//         Self {
//             receiver,
//             chunk_size,
//             report_schema,
//             fasta_reader,
//             join_handle,
//         }
//     }
//
//     fn get_context_data(&mut self, region: RegionCoordinates) -> Result<ContextData, Box<dyn Error>> {
//         let mut seq = bio::utils::Text::new();
//         if let Some(reader) = &mut self.fasta_reader {
//             reader.fetch(region.chr.as_str(), region.start as u64, region.end as u64)?;
//             reader.read(&mut seq)?;
//
//             // TODO Change position type in RegionCoordinates, ContextData,
//             //      GenomicPosition to something unified.
//             Ok(ContextData::from_sequence(String::from_utf8(seq)?, GenomicPosition::new(region.chr, region.start as usize)))
//         } else {
//             Err(Box::from("FASTA reader not provided"))
//         }
//     }
//
//     fn process_batch(&mut self, batch: DataFrame) -> PolarsResult<BsxBatch>  {
//         let report_schema = self.report_schema.get();
//         let batch = report_schema.bsx_schema_mutate(batch.lazy()).collect()?;
//         let do_align = report_schema.need_align();
//
//
//
//         todo!()
//     }
// }
//
// impl Iterator for ReportReader {
//
//     type Item = Result<BsxBatch, Box<dyn Error>>;
//
//     fn next(&mut self) -> Option<Self::Item> {
//         let next_batches = self.receiver.recv();
//         if next_batches.is_err() {
//             return None;
//         }
//         let next_batches = next_batches.unwrap();
//         todo!()
//     }
// }
