use _lib::data_structs::bsx_batch::BsxBatchMethods;
use _lib::io::bsx::read::BsxFileReader;
use _lib::io::bsx::write::{BsxIpcWriter, PolarsIpcCompression};
use _lib::io::report::read::ReportReaderBuilder;
use _lib::io::report::schema::ReportTypeSchema;
use _lib::io::report::write::ReportWriter;
use clap::{Args, ValueEnum};
use console::style;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

#[derive(Args, Debug, Clone)]
pub(crate) struct ReportArgs {
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads to use.",
        help_heading = "REPORT ARGS"
    )]
    n_threads: Option<usize>,
    #[arg(
        long,
        default_value_t = false,
        help = "Use less RAM, but elongate computation.",
        help_heading = "REPORT ARGS"
    )]
    low_memory: bool,
    #[arg(
        short,
        long = "chunk",
        default_value_t = 10_000,
        help = "Number of rows in the output batches (Important when converting to bsx format).",
        help_heading = "REPORT ARGS"
    )]
    chunk_size: usize,
    #[arg(
        long = "fa",
        help = "Path to the reference sequence file. Obligatory when converting BedGraph or Coverage.",
        help_heading = "REPORT ARGS"
    )]
    fasta_path: Option<PathBuf>,
    #[arg(
        long = "fai",
        help = "Path to the fasta index. Obligatory when converting BedGraph or Coverage.",
        help_heading = "REPORT ARGS"
    )]
    fai_path: Option<PathBuf>,
    #[arg(
        long,
        default_value_t = 8,
        help = "Number of batches to read simultaneously. Affects RAM usage.",
        help_heading = "REPORT ARGS"
    )]
    batch_per_read: usize,
    #[arg(long, default_value_t = 2 << 20, help = "Size of raw batches.")]
    batch_size: usize,
}

#[derive(Debug, Clone, ValueEnum, Eq, PartialEq)]
pub(crate) enum ConvertReportType {
    Bsx,
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
}

#[derive(Debug, Clone, ValueEnum)]
pub(crate) enum IpcCompression {
    LZ4,
    ZSTD,
    None,
}

pub fn run(
    input: PathBuf,
    output: PathBuf,
    from_type: ConvertReportType,
    into_type: ConvertReportType,
    ipc_compression: IpcCompression,
    report: ReportArgs,
) {
    if from_type == into_type {
        println!(
            "{}. Nothing to be done.",
            style("From report type matches into report type.")
        )
    }
    if matches!(
        from_type,
        ConvertReportType::BedGraph | ConvertReportType::Coverage
    ) && !(report.fasta_path.is_some() && report.fai_path.is_some())
    {
        eprintln!(
            "For {:?} report type conversion, both {}",
            from_type,
            style("fasta_path and fai_path must be set").red()
        )
    }
    if matches!(into_type, ConvertReportType::Bsx)
        && !(report.fasta_path.is_some() || report.fai_path.is_some())
    {
        eprintln!(
            "For {:?} report type conversion, either of {}",
            into_type,
            style("fasta_path and fai_path must be set").red()
        )
    }
    if !input.exists() {
        eprintln!("Path {} does not exist.", style(input.display()).red());
    }
    if !input.is_file() {
        eprintln!("Path {} is not a file.", style(input.display()).red());
    }
    if output.is_dir() {
        eprintln!(
            "Output path {} is a directory.",
            style(output.display()).red()
        );
    }

    let sink = BufWriter::new(File::create(output.clone()).unwrap_or_else(|e| {
        panic!(
            "Could not create output file {}: {}",
            output.to_string_lossy(),
            e
        )
    }));

    match from_type {
        ConvertReportType::Bismark
        | ConvertReportType::BedGraph
        | ConvertReportType::CgMap
        | ConvertReportType::Coverage => {
            let spinner = ProgressBar::new_spinner();
            spinner.set_style(
                ProgressStyle::default_spinner()
                    .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                    .template("{spinner} {msg}")
                    .expect("Failed to set spinner template"),
            );
            spinner.set_message("Processing...");

            let report_schema = match from_type {
                ConvertReportType::Bismark => ReportTypeSchema::Bismark,
                ConvertReportType::Coverage => ReportTypeSchema::Coverage,
                ConvertReportType::BedGraph => ReportTypeSchema::BedGraph,
                ConvertReportType::CgMap => ReportTypeSchema::CgMap,
                ConvertReportType::Bsx => unreachable!(),
            };
            let (fasta_path, fai_path) = match (report.fasta_path.clone(), report.fai_path.clone())
            {
                (Some(fasta_path), Some(fai_path)) => (Some(fasta_path), Some(fai_path)),
                _ => (None, None),
            };

            let report_reader = ReportReaderBuilder {
                report_type: report_schema,
                rechunk: true,
                n_threads: report.n_threads,
                low_memory: report.low_memory,
                n_rows: None,
                row_index: None,
                chunk_size: report.chunk_size,
                skip_rows_after_header: 0,
                fasta_path,
                fai_path,
                batch_per_read: report.batch_per_read,
                batch_size: report.batch_size,
            }
            .try_finish(File::open(input.clone()).unwrap_or_else(|e| {
                panic!("Could not open file {}: {}", input.to_string_lossy(), e)
            }))
            .unwrap_or_else(|e| {
                panic!("Could not create reader {}: {}", input.to_string_lossy(), e)
            });

            if into_type == ConvertReportType::Bsx {
                let compression = match ipc_compression {
                    IpcCompression::LZ4 => Some(PolarsIpcCompression::LZ4),
                    IpcCompression::ZSTD => Some(PolarsIpcCompression::ZSTD),
                    IpcCompression::None => None,
                };

                let mut writer = if report.fai_path.is_some() {
                    BsxIpcWriter::try_from_sink_and_fai(
                        sink,
                        report.fai_path.clone().unwrap(),
                        compression,
                        None,
                    )
                    .unwrap_or_else(|e| {
                        panic!(
                            "Could not create output file {}: {}",
                            output.to_string_lossy(),
                            e
                        )
                    })
                } else if report.fasta_path.is_some() {
                    BsxIpcWriter::try_from_sink_and_fasta(
                        sink,
                        report.fasta_path.clone().unwrap(),
                        compression,
                        None,
                    )
                    .unwrap_or_else(|e| {
                        panic!(
                            "Could not create output file {}: {}",
                            output.to_string_lossy(),
                            e
                        )
                    })
                } else {
                    unreachable!()
                };

                for (idx, batch) in report_reader.enumerate() {
                    let position = batch.last_position().expect("Error getting the position");
                    spinner.set_message(format!("Processing at {}", position));
                    if idx % 100 == 0 {
                        spinner.inc(1);
                    }
                    writer.write_batch(batch).expect("Could not write batch");
                }
                writer.close().expect("Could not close writer");
            } else {
                let into_report_schema = match into_type {
                    ConvertReportType::Bismark => ReportTypeSchema::Bismark,
                    ConvertReportType::Coverage => ReportTypeSchema::Coverage,
                    ConvertReportType::BedGraph => ReportTypeSchema::BedGraph,
                    ConvertReportType::CgMap => ReportTypeSchema::CgMap,
                    ConvertReportType::Bsx => unreachable!(),
                };

                let mut writer = ReportWriter::try_new(
                    sink,
                    into_report_schema,
                    report.n_threads.unwrap_or_default(),
                )
                .unwrap_or_else(|e| {
                    panic!(
                        "Could not create output writer {}: {}",
                        output.to_string_lossy(),
                        e
                    )
                });

                for (idx, batch) in report_reader.enumerate() {
                    let position = batch.last_position().expect("Error getting the position");
                    spinner.set_message(format!("Processing at {}", position));
                    if idx % 100 == 0 {
                        spinner.inc(1);
                    }
                    writer.write_batch(batch).expect("Could not write batch");
                }
            }

            spinner.finish_with_message(format!("Done: {}", style(output.display()).green()));
        }
        ConvertReportType::Bsx => {
            let reader = BsxFileReader::new(File::open(input.clone()).unwrap_or_else(|e| {
                panic!("Could not open file {}: {}", input.to_string_lossy(), e)
            }));

            let pb = ProgressBar::new(reader.blocks_total() as u64);

            // Set a custom style for the progress bar
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("{msg} [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                    .expect("Failed to set progress bar template")
                    .progress_chars("#>-"),
            );

            // Set the initial message
            pb.set_message("Processing...");

            let into_report_schema = match into_type {
                ConvertReportType::Bismark => ReportTypeSchema::Bismark,
                ConvertReportType::Coverage => ReportTypeSchema::Coverage,
                ConvertReportType::BedGraph => ReportTypeSchema::BedGraph,
                ConvertReportType::CgMap => ReportTypeSchema::CgMap,
                ConvertReportType::Bsx => unreachable!(),
            };

            let mut writer =
                ReportWriter::try_new(sink, into_report_schema, report.n_threads.unwrap_or(1))
                    .unwrap_or_else(|e| {
                        panic!(
                            "Could not create output writer {}: {}",
                            output.to_string_lossy(),
                            e
                        )
                    });

            for batch in reader {
                let batch = batch.expect("Could not decode ipc batch");

                pb.inc(1);

                let decoded_batch = batch.decode().expect("Could not decode ipc batch");
                writer
                    .write_batch(decoded_batch)
                    .expect("Could not write batch");
            }

            pb.finish_with_message(format!("Done: {}", style(output.display()).green()));
        }
    }
}
