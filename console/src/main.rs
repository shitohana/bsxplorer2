#![feature(path_add_extension)]

mod dmr_binom;
mod dmr_fast;
mod convert;
mod stats;

use crate::dmr_fast::{FilterArgs, SegmentationArgs};
use crate::dmr_binom::{SegmentationArgs2};
use clap::{Parser, Subcommand, ValueEnum};
use glob::glob;
use std::path::PathBuf;
use wild::ArgsOs;
use crate::convert::{ConvertReportType, IpcCompression, ReportArgs};

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

const DMR_ABOUT_FAST: &'static str = "BSXplorer DMR fast identification algorithm.";
const DMR_ABOUT: &'static str = "BSXplorer DMR Beta-binomial identification algorithm.";
const CONVERT_ABOUT: &'static str = "BSXplorer report type conversion tool";
#[derive(Subcommand, Debug)]
enum Commands {
    #[command(
        about = DMR_ABOUT_FAST,
        name = "dmrfast",
        after_help = include_str!("strings/dmr_fast_ahelp.txt"),
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
        #[arg(long, required = false, default_value_t = true, help = "Display progress bar (Disable if you need clean pipeline logs).")]
        progress: bool,
        #[arg(long, required = false, default_value_t = 1, help = "Number of threads to use.")]
        threads: usize,
        #[clap(flatten)]
        segmentation: SegmentationArgs,
    },

    #[command(
        about = DMR_ABOUT, 
        name = "dmr",
        after_help = include_str!("strings/dmr_binom_ahelp.txt"),
    )]
    Dmr2 {
        #[arg(value_parser, short='A', long, required = true, help = "Paths to BSX files of the first sample group.")]
        group_a: Vec<String>,
        #[arg(value_parser, short='B', long, required = true, help = "Paths to BSX files of the second sample group.")]
        group_b: Vec<String>,
        #[arg(short='o', long, required = true, help = "Prefix for the generated output files.")]
        output: PathBuf,
        #[clap(flatten)]
        filters: FilterArgs,
        #[arg(short, long, required = false, default_value_t = false, help = "Automatically confirm selected paths.")]
        force: bool,
        #[arg(long, required = false, default_value_t = true, help = "Display progress bar (Disable if you need clean pipeline logs).")]
        progress: bool,
        #[arg(long, required = false, default_value_t = 1, help = "Number of threads to use.")]
        threads: usize,
        #[clap(flatten)]
        segmentation: SegmentationArgs2,
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
    },

    #[command(
        name = "stats",
        about = "Compute methylation statistics.",
        after_help = include_str!("strings/stats_ahelp.txt"),
    )]
    Stats {
        #[clap(flatten)]
        stats: stats::StatsArgs,
    }
}

#[derive(Debug, Clone, ValueEnum)]
enum DmrContext {
    CG, CHG, CHH
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
            threads,
            progress,
        } => dmr_fast::run(filters, segmentation, group_a, group_b, output, force, threads, progress),

        Commands::Dmr2 {
            filters,
            segmentation,
            group_a,
            group_b,
            output ,
            force,
            threads,
            progress,
        } => dmr_binom::run(filters, segmentation, group_a, group_b, output, force, threads, progress),

        Commands::Convert {
            input,
            output,
            from_type,
            into_type,
            ipc_compression,
            report
        } => convert::run(input, output, from_type, into_type, ipc_compression, report),

        Commands::Stats { stats } => stats::run(stats),
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
