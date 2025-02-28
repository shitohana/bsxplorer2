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

use crate::{init_logger, init_rayon_threads, UtilsArgs};

#[derive(Args, Debug, Clone)]
pub struct ReportArgs {
    #[arg(help = "Path of the input file.")]
    input: PathBuf,
    #[arg(
        short = 'o',
        long,
        required = true,
        help = "Path for the generated output file."
    )]
    output: PathBuf,
    #[clap(short='f', long = "from", required = true, value_enum, default_value_t = ConvertReportType::Bismark)]
    from_type: ConvertReportType,
    #[clap(short='i', long = "into", required = true, value_enum, default_value_t = ConvertReportType::Bsx)]
    into_type: ConvertReportType,
    #[clap(short='C', long = "compression", required = false, value_enum, default_value_t = IpcCompression::ZSTD)]
    ipc_compression: IpcCompression,
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
pub enum ConvertReportType {
    Bsx,
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
}

#[derive(Debug, Clone, ValueEnum, Eq, PartialEq)]
pub enum IpcCompression {
    LZ4,
    ZSTD,
    None,
}

pub(crate) fn run(args: ReportArgs, utils: UtilsArgs) {
    init_rayon_threads(utils.threads).expect("Failed to initialyze global thread pool");
    init_logger(utils.verbose).expect("Failed to initialize logger");

    if args.from_type == args.into_type {
        println!(
            "{}. Nothing to be done.",
            style("From report type matches into report type.")
        )
    }
    if matches!(
        args.from_type,
        ConvertReportType::BedGraph | ConvertReportType::Coverage
    ) && !(args.fasta_path.is_some() && args.fai_path.is_some())
    {
        eprintln!(
            "For {:?} report type conversion, both {}",
            args.from_type,
            style("fasta_path and fai_path must be set").red()
        )
    }
    if matches!(args.into_type, ConvertReportType::Bsx)
        && !(args.fasta_path.is_some() || args.fai_path.is_some())
    {
        eprintln!(
            "For {:?} report type conversion, either of {}",
            args.into_type,
            style("fasta_path and fai_path must be set").red()
        )
    }
    if !args.input.exists() {
        eprintln!("Path {} does not exist.", style(args.input.display()).red());
    }
    if !args.input.is_file() {
        eprintln!("Path {} is not a file.", style(args.input.display()).red());
    }
    if args.output.is_dir() {
        eprintln!(
            "Output path {} is a directory.",
            style(args.output.display()).red()
        );
    }

    let sink = BufWriter::new(File::create(args.output.clone()).unwrap_or_else(|e| {
        panic!(
            "Could not create output file {}: {}",
            args.output.to_string_lossy(),
            e
        )
    }));

    match args.from_type {
        ConvertReportType::Bismark
        | ConvertReportType::BedGraph
        | ConvertReportType::CgMap
        | ConvertReportType::Coverage => {
            let spinner = if utils.progress {
                ProgressBar::new_spinner()
            } else {
                ProgressBar::hidden()
            };
            spinner.set_style(
                ProgressStyle::default_spinner()
                    .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                    .template("{spinner} {msg}")
                    .expect("Failed to set spinner template"),
            );
            spinner.set_message("Processing...");

            let report_schema = match args.from_type {
                ConvertReportType::Bismark => ReportTypeSchema::Bismark,
                ConvertReportType::Coverage => ReportTypeSchema::Coverage,
                ConvertReportType::BedGraph => ReportTypeSchema::BedGraph,
                ConvertReportType::CgMap => ReportTypeSchema::CgMap,
                ConvertReportType::Bsx => unreachable!(),
            };
            let (fasta_path, fai_path) = match (args.fasta_path.clone(), args.fai_path.clone()) {
                (Some(fasta_path), Some(fai_path)) => (Some(fasta_path), Some(fai_path)),
                _ => (None, None),
            };

            let report_reader = ReportReaderBuilder {
                report_type: report_schema,
                rechunk: true,
                n_threads: Some(utils.threads),
                low_memory: args.low_memory,
                n_rows: None,
                row_index: None,
                chunk_size: args.chunk_size,
                skip_rows_after_header: 0,
                fasta_path,
                fai_path,
                batch_per_read: args.batch_per_read,
                batch_size: args.batch_size,
            }
            .try_finish(File::open(args.input.clone()).unwrap_or_else(|e| {
                panic!(
                    "Could not open file {}: {}",
                    args.input.to_string_lossy(),
                    e
                )
            }))
            .unwrap_or_else(|e| {
                panic!(
                    "Could not create reader {}: {}",
                    args.input.to_string_lossy(),
                    e
                )
            });

            if args.into_type == ConvertReportType::Bsx {
                let compression = match args.ipc_compression {
                    IpcCompression::LZ4 => Some(PolarsIpcCompression::LZ4),
                    IpcCompression::ZSTD => Some(PolarsIpcCompression::ZSTD),
                    IpcCompression::None => None,
                };

                let mut writer = if args.fai_path.is_some() {
                    BsxIpcWriter::try_from_sink_and_fai(
                        sink,
                        args.fai_path.clone().unwrap(),
                        compression,
                        None,
                    )
                    .unwrap_or_else(|e| {
                        panic!(
                            "Could not create output file {}: {}",
                            args.output.to_string_lossy(),
                            e
                        )
                    })
                } else if args.fasta_path.is_some() {
                    BsxIpcWriter::try_from_sink_and_fasta(
                        sink,
                        args.fasta_path.clone().unwrap(),
                        compression,
                        None,
                    )
                    .unwrap_or_else(|e| {
                        panic!(
                            "Could not create output file {}: {}",
                            args.output.to_string_lossy(),
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
                let into_report_schema = match args.into_type {
                    ConvertReportType::Bismark => ReportTypeSchema::Bismark,
                    ConvertReportType::Coverage => ReportTypeSchema::Coverage,
                    ConvertReportType::BedGraph => ReportTypeSchema::BedGraph,
                    ConvertReportType::CgMap => ReportTypeSchema::CgMap,
                    ConvertReportType::Bsx => unreachable!(),
                };

                let mut writer = ReportWriter::try_new(sink, into_report_schema, utils.threads)
                    .unwrap_or_else(|e| {
                        panic!(
                            "Could not create output writer {}: {}",
                            args.output.to_string_lossy(),
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

            spinner.finish_with_message(format!("Done: {}", style(args.output.display()).green()));
        }
        ConvertReportType::Bsx => {
            let reader = BsxFileReader::new(File::open(args.input.clone()).unwrap_or_else(|e| {
                panic!(
                    "Could not open file {}: {}",
                    args.input.to_string_lossy(),
                    e
                )
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

            let into_report_schema = match args.into_type {
                ConvertReportType::Bismark => ReportTypeSchema::Bismark,
                ConvertReportType::Coverage => ReportTypeSchema::Coverage,
                ConvertReportType::BedGraph => ReportTypeSchema::BedGraph,
                ConvertReportType::CgMap => ReportTypeSchema::CgMap,
                ConvertReportType::Bsx => unreachable!(),
            };

            let mut writer = ReportWriter::try_new(sink, into_report_schema, utils.threads)
                .unwrap_or_else(|e| {
                    panic!(
                        "Could not create output writer {}: {}",
                        args.output.to_string_lossy(),
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

            pb.finish_with_message(format!("Done: {}", style(args.output.display()).green()));
        }
    }
}
