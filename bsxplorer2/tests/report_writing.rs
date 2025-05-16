// #![allow(unused)]
// use std::io::Write;

// use bio::io::fasta::Writer as FastaWriter;
// use bsxplorer2::data_structs::batch::{BsxBatch, BsxBatchMethods};
// #[cfg(feature = "compression")]
// use bsxplorer2::io::compression::Compression;
// use bsxplorer2::io::report::{ReportReaderBuilder,
//                              ReportType,
//                              ReportWriter};
// use polars::prelude::*;
// use rand::rngs::StdRng;
// use rstest::{fixture, rstest};

// mod common;
// use common::DemoReportBuilder;
// #[fixture]
// fn report_data() -> DemoReportBuilder<StdRng> {
//     DemoReportBuilder::new(123_456, 20, 15.0, 0.5, Some(42))
// }

// const N_CHR: usize = 3;
// const CHUNK_SIZE: usize = 10000;

// #[rstest]
// #[case::bismark_batched(ReportType::Bismark, true)]
// #[case::bismark(ReportType::Bismark, false)]
// #[case::cgmap(ReportType::CgMap, false)]
// #[case::bedgraph(ReportType::BedGraph, false)]
// #[case::coverage(ReportType::Coverage, false)]
// fn test_report_writing_reading_roundtrip(
//     #[case] report_type: ReportType,
//     #[case] write_batch: bool,
//     mut report_data: DemoReportBuilder<StdRng>,
// ) -> anyhow::Result<()> {
//     // 1. Generate data

//     use bsxplorer2::data_structs::batch::LazyBsxBatch;
//     use common::compare_batches;
//     use tempfile::NamedTempFile;

//     let n_batches = N_CHR; // Use N_CHR for consistency with reading test
//     let mut original_batches: Vec<BsxBatch> = Vec::new();
//     let mut original_dfs_in_report_format: Vec<DataFrame> = Vec::new();

//     // Need sequence file for alignment-required types (BedGraph, Coverage)
//     // during reading
//     let sequence_file = tempfile::NamedTempFile::new()?;
//     let mut sequence_writer = FastaWriter::new(sequence_file.reopen()?);

//     for _ in 0..n_batches {
//         let (record, original_batch) = report_data.next().unwrap();
//         original_batches.push(original_batch.clone());

//         // Generate report format DF from the batch
//         let lazy_batch = LazyBsxBatch::from(original_batch);
//         let report_df = lazy_batch.into_report(&report_type)?;
//         original_dfs_in_report_format.push(report_df);

//         // Write sequence for alignment
//         sequence_writer.write_record(&record)?;
//     }
//     sequence_writer.flush()?;
//     drop(sequence_writer); // Close sequence file handle

//     // 2. Write data to buffer
//     let buffer = NamedTempFile::new()?;

//     // Handle BedGraph header manually for writing test if needed
//     if report_type == ReportType::BedGraph {
//         writeln!(buffer.reopen()?, "track type=bedGraph")?;
//     }

//     #[cfg(not(feature = "compression"))]
//     let mut writer = ReportWriter::try_new(buffer.reopen()?, report_type, 8)?;
//     #[cfg(feature = "compression")]
//     let mut writer = ReportWriter::try_new(
//         buffer.reopen()?,
//         report_type,
//         8,
//         Compression::None,
//         None,
//     )?;

//     if !write_batch {
//         for df in &original_dfs_in_report_format {
//             writer.write_df(df)?;
//         }
//     }
//     else {
//         for batch in original_batches.clone() {
//             writer.write_batch(batch)?;
//         }
//     }

//     // Writer is dropped here, flushing happens

//     // 3. Read data back from buffer

//     let mut report_reader_builder: ReportReaderBuilder<BsxBatch> =
//         ReportReaderBuilder::default()
//         .with_chunk_size(CHUNK_SIZE) // Use same chunk size as reading test
//         .with_report_type(report_type);

//     // Add fasta path if alignment is needed by the reader for this report type
//     if report_type.need_align() {
//         report_reader_builder = report_reader_builder
//             .with_fasta_path(Some(sequence_file.path().to_path_buf()));
//     }
//     // Note: We are reading into `BsxBatch` (decoded), which doesn't require
//     // chr_dtype from fasta for non-aligning types. The builder handles this.

//     let report_reader =
//         report_reader_builder.build_from_handle(Box::new(buffer.reopen()?))?;

//     // 4. Collect read batches
//     let read_batches: Vec<BsxBatch> = report_reader
//         .into_iter()
//         .collect::<Result<_, _>>()?;

//     // 5. Combine original and read batches by chromosome (as reader outputs one
//     //    batch per chr after final chunk)
//     // Note: The reader combines chunks of the same chromosome. So we need to
//     // combine original batches too.
//     let mut combined_original_batches: Vec<BsxBatch> = Vec::new();
//     let mut current_original_chr = String::new();
//     let mut current_original_batch: Option<BsxBatch> = None;

//     // Original batches are ordered by generation, which is per chromosome
//     for batch in original_batches {
//         let batch_chr = batch.chr_val().to_string();
//         if batch_chr != current_original_chr {
//             if let Some(completed_batch) = current_original_batch.take() {
//                 combined_original_batches.push(completed_batch);
//             }
//             current_original_chr = batch_chr;
//             current_original_batch = Some(batch);
//         }
//         else {
//             if let Some(ref mut cb) = current_original_batch {
//                 cb.extend(&batch)?;
//             }
//             else {
//                 current_original_batch = Some(batch);
//             }
//         }
//     }
//     if let Some(completed_batch) = current_original_batch {
//         combined_original_batches.push(completed_batch);
//     }

//     // The reader already returns combined batches per chromosome (at least for
//     // the final batches) The logic in the reading test combines chunks from
//     // the reader into per-chromosome batches. Let's reuse that combination
//     // logic here for the read batches.

//     let mut combined_read_batches: Vec<BsxBatch> = Vec::new();
//     let mut current_read_chr = String::new();
//     let mut current_read_batch = None;

//     for batch in read_batches {
//         let batch_chr = batch.chr_val().to_string();

//         if batch_chr != current_read_chr {
//             // New chromosome encountered
//             if let Some(completed_batch) = current_read_batch.take() {
//                 combined_read_batches.push(completed_batch);
//             }
//             current_read_chr = batch_chr;
//             current_read_batch = Some(batch);
//         }
//         else {
//             // Same chromosome, extend the current batch
//             if let Some(ref mut cb) = current_read_batch {
//                 cb.extend(&batch)?;
//             }
//             else {
//                 current_read_batch = Some(batch);
//             }
//         }
//     }

//     // Add the last batch
//     if let Some(completed_batch) = current_read_batch {
//         combined_read_batches.push(completed_batch);
//     }

//     // 6. Compare original vs read combined batches
//     assert_eq!(
//         combined_read_batches.len(),
//         combined_original_batches.len(),
//         "Number of chromosomes read back does not match original"
//     );

//     for (original, read) in combined_original_batches
//         .iter()
//         .zip(combined_read_batches.iter())
//     {
//         println!("{read:?}");
//         println!("{original:?}");
//         compare_batches(original, read, &report_type)?;
//     }

//     Ok(())
// }
