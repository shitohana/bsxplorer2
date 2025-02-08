use std::io::Write;
use clap::ColorChoice;
use std::cmp::PartialEq;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use clap::{Args, Parser, Subcommand, ValueEnum};
use std::iter::repeat_n;
// For parsing command-line arguments.
use console::style;
// For colored text output.
use dialoguer::Confirm;
use _lib::tools::dmr::{MethyLassoConfig, MethylLassoRunConfig, RegionType, SegmentModel};
use _lib::utils::types::{Context, IPCEncodedEnum};
// For interactive user confirmations.
use glob::glob;
use std::path::PathBuf;
use indicatif::{ProgressBar, ProgressStyle};
use wild::ArgsOs;
use _lib::data_structs::bsx_batch::BsxBatchMethods;
use _lib::io::bsx::region_read::BsxFileReader;
use _lib::io::bsx::write::{BsxIpcWriter, PolarsIpcCompression};
use _lib::io::report::read::ReportReaderBuilder;
use _lib::io::report::schema::ReportTypeSchema;
use _lib::io::report::write::ReportWriter;

#[derive(Parser, Debug)]
#[command(
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = env!("CARGO_PKG_DESCRIPTION"),
    long_about = None,
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

const DMR_ABOUT: &'static str = "BSXplorer DMR identification algorithm";
const CONVERT_ABOUT: &'static str = "BSXplorer report type conversion tool";
#[derive(Subcommand, Debug)]
enum Commands {
    #[command(
        about = DMR_ABOUT, 
        name = "dmr",
        after_help = include_str!("strings/dmr_ahelp.txt"),
    )]
    Dmr {
        #[arg(value_parser, short='A', long, required = true, help = "Paths to BSX files of the first sample group.")]
        group_a: Vec<String>,
        #[arg(value_parser, short='B', long, required = true, help = "Paths to BSX files of the second sample group.")]
        group_b: Vec<String>,
        #[arg(short='o', long, required = true, help = "Path for the generated output file.")]
        output: PathBuf,
        #[clap(flatten)]
        filters: FilterArgs,
        #[arg(short, long, required = false, default_value_t = false, help = "Automatically confirm selected paths.")]
        force: bool,
        #[clap(flatten)]
        segmentation: SegmentationArgs,
    },

    #[command(
        name = "convert",
        about = CONVERT_ABOUT,
        after_help = include_str!("strings/convert_ahelp.txt"),
    )]
    Convert {
        #[arg(help = "Path of the input file.")]
        input: PathBuf,
        #[arg(short='o', long, required = true, help = "Path for the generated output file.")]
        output: PathBuf,
        #[clap(short='f', long = "from", required = true, value_enum, default_value_t = ConvertReportType::Bismark)]
        from_type: ConvertReportType,
        #[clap(short='i', long = "into", required = true, value_enum, default_value_t = ConvertReportType::Bsx)]
        into_type: ConvertReportType,
        #[clap(short='C', long = "compression", required = false, value_enum, default_value_t = IpcCompression::ZSTD)]
        ipc_compression: IpcCompression,
        #[clap(flatten)]
        report: ReportArgs,
    }
}

#[derive(Debug, Clone, ValueEnum, Eq, PartialEq)]
enum ConvertReportType {
    Bsx,
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
}

#[derive(Debug, Clone, ValueEnum)]
enum IpcCompression {
    LZ4,
    ZSTD,
    None
}

#[derive(Debug, Clone, ValueEnum)]
enum DmrContext {
    CG, CHG, CHH
}

#[derive(Args, Debug , Clone)]
struct ReportArgs {
    #[arg(short='t', long="threads", help = "Number of threads to use.", help_heading="REPORT ARGS")]
    n_threads: Option<usize>,
    #[arg(long, default_value_t = false, help = "Use less RAM, but elongate computation.", help_heading="REPORT ARGS")]
    low_memory: bool,
    #[arg(short, long = "chunk", default_value_t = 10_000, help = "Number of rows in the output batches (Important when converting to bsx format).", help_heading="REPORT ARGS")]
    chunk_size: usize,
    #[arg(long = "fa", help = "Path to the reference sequence file. Obligatory when converting BedGraph or Coverage.", help_heading="REPORT ARGS")]
    fasta_path: Option<PathBuf>,
    #[arg(long = "fai", help = "Path to the fasta index. Obligatory when converting BedGraph or Coverage.", help_heading="REPORT ARGS")]
    fai_path: Option<PathBuf>,
    #[arg(long, default_value_t = 8, help = "Number of batches to read simultaneously. Affects RAM usage.", help_heading="REPORT ARGS")]
    batch_per_read: usize,
    #[arg(long, default_value_t = 2 << 20, help = "Size of raw batches.")]
    batch_size: usize,
}


#[derive(Args, Debug, Clone)]
struct FilterArgs {
    #[clap(short, long, value_enum, default_value_t = DmrContext::CG, help_heading="FILTER ARGS", help = "Select cytosine methylation context.")]
    context: DmrContext,
    #[arg(short, long, default_value_t = 0, help_heading="FILTER ARGS", help = "Set missing values threshold. Cytosines with no data in > n_samples will be discarded.")]
    n_missing: usize,
    #[arg(short='v', long, default_value_t = 5, help_heading="FILTER ARGS", help = "Cytosines with coverage below threshold will be discarded.")]
    min_coverage: i16,
    #[arg(short, long, default_value_t = 0.05, help_heading="FILTER ARGS", help = "P-value for Welch t-Test for DMR identification. Regions with p-value above threshold will be discarded.")]
    p_value: f64,
    #[arg(short, long, default_value_t = 10, help_heading="FILTER ARGS", help = "DMRs with cytosine count below threshold will be discarded.")]
    min_cytosines: usize,
    #[arg(short, long, default_value_t = 0.05, help_heading="FILTER ARGS", help = "DMRs with difference between samples below threshold will be discarded.")]
    diff_threshold: f64,
}

#[derive(Args, Debug, Clone)]
struct SegmentationArgs {
    #[arg(long, default_value_t = 100, help_heading="SEGMENTATION ARGS", help = "Number of optimal lambda approximation iterations.")]
    nlambda: usize,
    #[arg(long, default_value_t = 1e-3, help_heading="SEGMENTATION ARGS", help = "Initial lambda value. Affects segments length. Must be at least 10 times less than expected methylation density range.")]
    initial_l: f64,
    #[arg(long, default_value_t = 1e-6, help_heading="SEGMENTATION ARGS", help = "Tolerance of segment identifying after denoising step (should be low).")]
    tolerance: f64,
    #[arg(long, default_value_t = 1e5, help_heading="SEGMENTATION ARGS", help = "Segment length penalty in optimal lambda approximation. Longer segments -> less penalty.")]
    length_penalty: f64,
    #[arg(long, default_value_t = 1e4, help_heading="SEGMENTATION ARGS", help = "Segments count penalty in optimal lambda approximation. Fewer segments -> less penalty.")]
    count_penalty: f64,
    #[arg(long, default_value_t = 1e-2, help_heading="SEGMENTATION ARGS", help = "Welch t-Test P-value for segments merging step. Smaller p-value -> more iterations, less falsely merged segments.")]
    merge_p: f64,
}


fn main() {
    let args: ArgsOs = wild::args_os();
    let cli = Cli::parse_from(args);
    // Dispatch the command based on the provided subcommand.
    match cli.command {
        Commands::Dmr { 
            filters, 
            segmentation, 
            group_a, 
            group_b, 
            output ,
            force,
        } => {
            let a_paths = expand_wildcards(group_a);
            let b_paths = expand_wildcards(group_b);

            if !force {
                let prompt = format!("Do you want to proceed with the following paths?\n\nGroup A: {:?}\nGroup B: {:?}\nOutput:{:?}", a_paths, b_paths, output);
                let confirmed = Confirm::new()
                    .with_prompt(prompt)
                    .default(true)
                    .interact()
                    .unwrap_or(false);

                if !confirmed {
                    println!("{}", style("Process aborted by the user.").red());
                    return;
                }
            }
            
            for path in a_paths.iter().chain(b_paths.iter()) {
                if !path.exists() {
                    eprintln!("Path {} does not exist.", style(path.display()).red());
                }
                if !path.is_file() {
                    eprintln!("Path {} is not a file.", style(path.display()).red());
                }
            }
            
            if output.is_dir() {
                eprintln!("Output path {} is a directory.", style(output.display()).red());
            }
            
            let context = match filters.context {
                DmrContext::CG => Context::CG,
                DmrContext::CHG => Context::CHG,
                DmrContext::CHH => Context::CHH,
            };
            let left_labels = repeat_n("A".to_string(), a_paths.len());
            let right_labels = repeat_n("B".to_string(), b_paths.len());

            let sample_paths = a_paths.into_iter().chain(b_paths.into_iter()).map(|p| p.to_string_lossy().to_string()).collect::<Vec<_>>();
            let sample_labels = left_labels.into_iter().chain(right_labels.into_iter()).collect::<Vec<_>>();

            let run_config = MethylLassoRunConfig {
                analysis_config: MethyLassoConfig {
                    context,
                    n_missing: filters.n_missing,
                    min_coverage: filters.min_coverage,
                    diff_threshold: filters.diff_threshold,
                    p_value: filters.p_value,
                    type_density: Default::default(),
                    min_cpgs: filters.min_cytosines,
                    segment_model: SegmentModel {
                        min_seg_length: filters.min_cytosines,
                        num_lambdas: segmentation.nlambda,
                        lambda_low: segmentation.initial_l,
                        penalty_weight: segmentation.length_penalty,
                        seg_count_weight: segmentation.count_penalty,
                        merge_pvalue: segmentation.merge_p,
                        segmentation_tol: segmentation.tolerance,
                    },
                },

                sample_paths,
                sample_labels,
                selected_regions: vec![RegionType::DMR],
                output: output.to_string_lossy().to_string(),
            };

            let files = run_config.sample_paths
                .iter()
                .map(|path| File::open(path).expect(format!("Could not open file {}", path).as_str()))
                .map(BufReader::new)
                .collect::<Vec<_>>();

            let iterator = run_config.analysis_config.clone().finish(
                run_config.sample_labels
                    .iter()
                    .cloned()
                    .zip(files.into_iter())
                    .collect::<Vec<_>>()
            ).expect("Failed to initialize sample reader");

            let spinner = ProgressBar::new_spinner();
            spinner.set_style(ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                .template("{spinner} {msg}")
                .expect("Failed to set spinner template"));
            spinner.set_message("Processing...");

            let mut sink = BufWriter::new(File::create(&run_config.output).expect("Failed to open output file"));
            
            let mut dmr_count = 0;
            for region in iterator {
                if run_config.selected_regions.contains(&region.region_type) {
                    spinner.set_message(format!("Found {} DMRs", style(dmr_count).green().bold()));
                    dmr_count += 1;
                    
                    if dmr_count % 10 == 0 {
                        spinner.inc(1);
                    }
                    
                    let meth_diff = region.mean_right - region.mean_left;
                    let overall = (region.mean_left + region.mean_right) / 2.0;
                    let n_cytosines = region.pairwise_sites.len();
                    let ref_id = region.pairwise_sites.left.chr.as_str();
                    let start = region.pairwise_sites.left.positions.first().unwrap();
                    let end = region.pairwise_sites.left.positions.last().unwrap();
                    let region_type = region.region_type;
                    let mean_a = region.mean_left;
                    let mean_b = region.mean_right;
                    let pvalue = region.p_value;

                    let line = format!("{ref_id}\t{start}\t{end}\t{n_cytosines}\t{region_type:?}\t{pvalue:.3e}\t{mean_a:.5}\t{mean_b:.5}\t{meth_diff:.5}\t{overall}");
                    writeln!(sink, "{line}").expect("Failed to write line");
                }
            }
        },

        Commands::Convert {
            input,
            output,
            from_type,
            into_type,
            ipc_compression,
            report
        } => {
            if from_type == into_type {
                println!("{}. Nothing to be done.", style("From report type matches into report type."))
            }
            if matches!(from_type, ConvertReportType::BedGraph | ConvertReportType::Coverage) &&
                !(report.fasta_path.is_some() && report.fai_path.is_some()) {
                eprintln!("For {:?} report type conversion, both {}", from_type, style("fasta_path and fai_path must be set").red())
            }
            if matches!(into_type, ConvertReportType::Bsx) &&
                !(report.fasta_path.is_some() || report.fai_path.is_some()) {
                eprintln!("For {:?} report type conversion, either of {}", into_type, style("fasta_path and fai_path must be set").red())
            }
            if !input.exists() {
                eprintln!("Path {} does not exist.", style(input.display()).red());
            }
            if !input.is_file() {
                eprintln!("Path {} is not a file.", style(input.display()).red());
            }
            if output.is_dir() {
                eprintln!("Output path {} is a directory.", style(output.display()).red());
            }

            let sink = BufWriter::new(
                File::create(output.clone()).unwrap_or_else(|e| panic!("Could not create output file {}: {}", output.to_string_lossy(), e)),
            );

            match from_type {
                ConvertReportType::Bismark |
                ConvertReportType::BedGraph |
                ConvertReportType::CgMap |
                ConvertReportType::Coverage => {
                    let spinner = ProgressBar::new_spinner();
                    spinner.set_style(ProgressStyle::default_spinner()
                        .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                        .template("{spinner} {msg}")
                        .expect("Failed to set spinner template"));
                    spinner.set_message("Processing...");

                    let report_schema = match from_type {
                        ConvertReportType::Bismark => ReportTypeSchema::Bismark,
                        ConvertReportType::Coverage => ReportTypeSchema::Coverage,
                        ConvertReportType::BedGraph => ReportTypeSchema::BedGraph,
                        ConvertReportType::CgMap => ReportTypeSchema::CgMap,
                        ConvertReportType::Bsx => unreachable!()
                    };
                    let (fasta_path, fai_path) = match (report.fasta_path.clone(), report.fai_path.clone()) {
                        (Some(fasta_path), Some(fai_path)) => (Some(fasta_path), Some(fai_path)),
                        _ => (None, None)
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
                    }.try_finish(
                        File::open(input.clone())
                            .unwrap_or_else(|e| panic!("Could not open file {}: {}", input.to_string_lossy(), e)),
                    ).unwrap_or_else(
                        |e| panic!("Could not create reader {}: {}", input.to_string_lossy(), e)
                    );

                    if into_type == ConvertReportType::Bsx {
                        let compression = match ipc_compression {
                            IpcCompression::LZ4 => Some(PolarsIpcCompression::LZ4),
                            IpcCompression::ZSTD => Some(PolarsIpcCompression::ZSTD),
                            IpcCompression::None => None
                        };

                        let mut writer = if report.fai_path.is_some() {
                            BsxIpcWriter::try_from_sink_and_fai(
                                sink,
                                report.fai_path.clone().unwrap(),
                                compression,
                                None
                            ).unwrap_or_else(|e| panic!("Could not create output file {}: {}", output.to_string_lossy(), e))
                        } else if report.fasta_path.is_some() {
                            BsxIpcWriter::try_from_sink_and_fasta(
                                sink,
                                report.fasta_path.clone().unwrap(),
                                compression,
                                None
                            ).unwrap_or_else(|e| panic!("Could not create output file {}: {}", output.to_string_lossy(), e))
                        } else {
                            unreachable!()
                        };

                        for (idx, batch) in report_reader.enumerate() {
                            let position = batch.last_position()
                                .expect("Error getting the position");
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
                            ConvertReportType::Bsx => unreachable!()
                        };

                        let mut writer = ReportWriter::try_new(
                            sink, into_report_schema, report.n_threads.unwrap_or_default()
                        ).unwrap_or_else(|e| panic!("Could not create output writer {}: {}", output.to_string_lossy(), e));

                        for (idx, batch) in report_reader.enumerate() {
                            let position = batch.last_position()
                                .expect("Error getting the position");
                            spinner.set_message(format!("Processing at {}", position));
                            if idx % 100 == 0 {
                                spinner.inc(1);
                            }
                            writer.write_batch(batch).expect("Could not write batch");
                        }
                    }

                    spinner.finish_with_message(format!("Done: {}", style(output.display()).green()));
                },
                ConvertReportType::Bsx => {
                    let reader = BsxFileReader::new(
                        File::open(input.clone()).unwrap_or_else(|e| panic!("Could not open file {}: {}", input.to_string_lossy(), e)),
                    );

                    let pb = ProgressBar::new(reader.blocks_total() as u64);

                    // Set a custom style for the progress bar
                    pb.set_style(ProgressStyle::default_bar()
                        .template("{msg} [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                        .expect("Failed to set progress bar template")
                        .progress_chars("#>-"));

                    // Set the initial message
                    pb.set_message("Processing...");

                    let into_report_schema = match into_type {
                        ConvertReportType::Bismark => ReportTypeSchema::Bismark,
                        ConvertReportType::Coverage => ReportTypeSchema::Coverage,
                        ConvertReportType::BedGraph => ReportTypeSchema::BedGraph,
                        ConvertReportType::CgMap => ReportTypeSchema::CgMap,
                        ConvertReportType::Bsx => unreachable!()
                    };

                    let mut writer = ReportWriter::try_new(
                        sink, into_report_schema, report.n_threads.unwrap_or(1)
                    ).unwrap_or_else(|e| panic!("Could not create output writer {}: {}", output.to_string_lossy(), e));

                    for batch in reader {
                        let batch = batch
                            .expect("Could not decode ipc batch");

                        pb.inc(1);

                        let decoded_batch =
                            batch.decode()
                            .expect("Could not decode ipc batch");
                        writer.write_batch(decoded_batch).expect("Could not write batch");
                    }

                    pb.finish_with_message(format!("Done: {}", style(output.display()).green()));
                },
            }
        }
    }
}

fn expand_wildcards(paths: Vec<String>) -> Vec<PathBuf> {
    let mut expanded_paths = Vec::new();

    for path in paths {
        if path.contains('*') || path.contains('?') {
            // Expand wildcard using glob
            match glob(&path) {
                Ok(matches) => {
                    for entry in matches.filter_map(Result::ok) {
                        expanded_paths.push(entry);
                    }
                }
                Err(e) => eprintln!("Error processing wildcard '{}': {}", path, e),
            }
        } else {
            // If not a wildcard, push the path as-is
            expanded_paths.push(PathBuf::from(path));
        }
    }

    expanded_paths
}
