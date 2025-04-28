use std::io::Cursor;
use std::ops::BitOr;

use anyhow::bail;
use bio::io::fasta::Writer as FastaWriter;
use bsxplorer2::data_structs::batch::colnames::*;
use bsxplorer2::data_structs::batch::LazyBsxBatch;
use bsxplorer2::data_structs::batch::{BsxBatch, BsxBatchMethods};
use bsxplorer2::io::report::ReportReaderBuilder;
use bsxplorer2::io::report::ReportTypeSchema;
use rand::rngs::StdRng;
use rstest::*;

use polars::prelude::*;
use std::io::Write;

mod common;
use common::DemoReportBuilder;
#[fixture]
fn report_data() -> DemoReportBuilder<StdRng> {
    DemoReportBuilder::new(123_456, 20, 15.0, 0.5, Some(42))
}

fn compare_batches(
    original: &BsxBatch,
    read: &BsxBatch,
    report_type: &ReportTypeSchema,
) -> anyhow::Result<()> {
    if original.chr_val()? != read.chr_val()? {
        bail!(
            "Chromosomes differ for original: {}\nread: {}",
            original.data(),
            read.data()
        );
    }
    if original.height() != read.height() {
        bail!(
            "Heights differ for original: {}\nread: {}",
            original.data(),
            read.data()
        );
    }
    let (original_df, read_df) =
        if matches!(report_type, ReportTypeSchema::BedGraph) {
            (
                original
                    .data()
                    .drop_many([COUNT_M_NAME, COUNT_TOTAL_NAME]),
                read.data()
                    .drop_many([COUNT_M_NAME, COUNT_TOTAL_NAME]),
            )
        } else {
            (original.data().clone(), read.data().clone())
        };
    if !original_df.equals_missing(&read_df) {
        let diff_mask = original_df
            .materialized_column_iter()
            .zip(read_df.materialized_column_iter())
            .map(|(orig, read)| !orig.equal(read).unwrap())
            .reduce(|acc, new| acc.bitor(new))
            .unwrap();
        let orig_diff = original_df.filter(&diff_mask)?;
        let read_diff = read_df.filter(&diff_mask)?;
        bail!(
            "Data differs for original: {}\nread: {}",
            orig_diff,
            read_diff
        );
    }

    Ok(())
}

const N_CHR: usize = 3;
const CHUNK_SIZE: usize = 10000;

/// Tests the report reading functionality with different report types and alignment options.
///
/// This test verifies that the ReportReader correctly processes different report formats
/// and handles alignment with reference sequences when requested.
///
/// # Test Cases
/// - Tests various report types (Bismark, CgMap, BedGraph, Coverage)
/// - Tests with and without alignment to reference sequences
///
/// # Parameters
/// - `report_type`: The schema/format of the report being tested
/// - `do_alignment`: Whether to perform alignment with reference sequences
/// - `report_data`: Fixture providing demo methylation data
#[rstest]
#[case::bismark_align(ReportTypeSchema::Bismark, true)]
#[case::cgmap_align(ReportTypeSchema::CgMap, true)]
#[case::bedgraph_align(ReportTypeSchema::BedGraph, true)]
#[case::coverage_align(ReportTypeSchema::Coverage, true)]
#[case::bismark_no_align(ReportTypeSchema::Bismark, false)]
#[case::cgmap_no_align(ReportTypeSchema::CgMap, false)]
fn test_report_reading_with_alignment(
    #[case] report_type: ReportTypeSchema,
    #[case] do_alignment: bool,
    report_data: DemoReportBuilder<StdRng>,
) -> anyhow::Result<()> {
    // Initialize buffers for report data and sequence data
    let mut report_buffer = Vec::new();
    let sequence_file = tempfile::NamedTempFile::new()?;

    // Add header for BedGraph format if needed
    if report_type == ReportTypeSchema::BedGraph {
        writeln!(&mut report_buffer, "track type=bedGraph")?;
    }

    // Set up writers for the report and sequence data
    let mut csv_writer = CsvWriter::new(&mut report_buffer)
        .with_separator(b'\t')
        .include_header(false)
        .batched(&report_type.schema())?;
    let mut sequence_writer = FastaWriter::new(sequence_file.reopen()?);

    // Store original batches for later comparison
    let mut original_batches = Vec::new();

    // Generate test data for multiple chromosomes
    let mut demo_report_iter = report_data;
    for _chr in 0..N_CHR {
        // Get a record and batch for a chromosome
        let (record, original_batch) = demo_report_iter.next().unwrap();

        // Save the original batch for later comparison
        original_batches.push(original_batch.clone());

        // Convert to the appropriate report format
        let lazy_batch = LazyBsxBatch::from(original_batch.clone());
        let report_df = lazy_batch.into_report(&report_type)?;

        // Write the data to our buffers
        csv_writer.write_batch(&report_df)?;
        sequence_writer.write_record(&record)?;
    }

    // Finalize and close the writers
    csv_writer.finish()?;
    sequence_writer.flush()?;
    drop(sequence_writer);

    // Create a reader for the generated report data
    let report_handle = Cursor::new(report_buffer);
    let mut report_reader_builder = ReportReaderBuilder::<BsxBatch>::default()
        .with_chunk_size(CHUNK_SIZE)
        .with_report_type(report_type);

    // Configure alignment if requested
    if do_alignment {
        report_reader_builder = report_reader_builder
            .with_fasta_path(sequence_file.path().to_path_buf());
    }

    // Build the report reader and create an iterator
    let report_reader =
        report_reader_builder.build_from_handle(Box::new(report_handle))?;
    let report_reader_iter = report_reader.into_iter();

    // Process the batches and count final batches
    // A final batch is identified by having fewer rows than CHUNK_SIZE
    let mut final_batch_num = 0;
    let mut read_batches = Vec::new();

    for batch in report_reader_iter {
        let batch = batch?;
        read_batches.push(batch.clone());
        if batch.height() != CHUNK_SIZE {
            final_batch_num += 1;
        }
    }

    // Verify we have the expected number of final batches (one per chromosome)
    assert_eq!(final_batch_num, N_CHR);

    // Compare read data with original data
    // Note: We may need to combine read batches if they were split due to chunk size
    let mut combined_read_batches = Vec::new();
    let mut current_chr = String::new();
    let mut current_batch = None;

    for batch in read_batches {
        let batch_chr = batch.chr_val()?.to_string();

        if batch_chr != current_chr {
            // New chromosome encountered
            if let Some(completed_batch) = current_batch.take() {
                combined_read_batches.push(completed_batch);
            }
            current_chr = batch_chr;
            current_batch = Some(batch);
        } else {
            // Same chromosome, extend the current batch
            if let Some(ref mut cb) = current_batch {
                cb.extend(&batch)?;
            } else {
                current_batch = Some(batch);
            }
        }
    }

    // Add the last batch
    if let Some(completed_batch) = current_batch {
        combined_read_batches.push(completed_batch);
    }

    // Verify we have the same number of chromosomes
    assert_eq!(combined_read_batches.len(), original_batches.len());

    // Compare each chromosome's data
    for (original, read) in original_batches
        .iter()
        .zip(combined_read_batches.iter())
    {
        // Compare key columns - this may need adjustment based on what's preserved in the report format
        compare_batches(original, read, &report_type)?;
    }

    Ok(())
}
