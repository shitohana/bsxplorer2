#![feature(path_add_extension)]
pub mod convert;
mod dmr;
pub mod stats;
pub mod utils;
use crate::convert::ReportArgs;
use clap::{Parser, Subcommand};
use serde::Serialize;
pub(crate) use utils::*;
use wild::ArgsOs;

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

const DMR_ABOUT: &'static str = "BSXplorer DMR identification algorithm.";
const CONVERT_ABOUT: &'static str = "BSXplorer report type conversion tool";
#[derive(Subcommand, Debug)]
enum Commands {
    #[command(
        about = DMR_ABOUT,
        name = "dmr",
    )]
    Dmr {
        #[clap(flatten)]
        utils: UtilsArgs,
        #[clap(flatten)]
        args: dmr::DmrArgs,
    },

    #[command(
        name = "convert",
        about = CONVERT_ABOUT,
        after_help = include_str!("strings/convert_ahelp.txt"),
    )]
    Convert {
        #[clap(flatten)]
        args: ReportArgs,
        #[clap(flatten)]
        utils: UtilsArgs,
    },

    #[command(
        name = "stats",
        about = "Compute methylation statistics.",
        after_help = include_str!("strings/stats_ahelp.txt"),
    )]
    Stats {
        #[clap(flatten)]
        stats: stats::StatsArgs,
        #[clap(flatten)]
        utils: UtilsArgs,
    },
}

fn main() {
    let args: ArgsOs = wild::args_os();
    let cli = Cli::parse_from(args);
    match cli.command {
        Commands::Dmr { args, utils } => dmr::run(args, utils),

        Commands::Convert { args, utils } => convert::run(args, utils),

        Commands::Stats { stats, utils } => stats::run(stats, utils),
    }
}
