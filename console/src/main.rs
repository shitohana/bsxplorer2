mod convert;
mod dimred;
mod dmr;
mod sort;
mod stats;
mod strings;
mod utils;
mod validate;

use clap::{
    Parser,
    Subcommand,
};
use convert::{
    FromBsxConvert,
    R2RConvert,
    ToBsxConvert,
};
use dimred::{
    DimRedArgs,
    MergeArgs,
};
use sort::SortArgs;
use validate::ValidateArgs;
use wild::ArgsOs;

pub(crate) trait PipelineCommand {
    fn run(&self) -> anyhow::Result<()>;
}

#[derive(Parser, Debug)]
#[command(
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = env!("CARGO_PKG_DESCRIPTION"),
    long_about = None,)]
struct Cli {
    #[command(subcommand)]
    command: MainMenu,
}

#[derive(Subcommand, Debug)]
enum MainMenu {
    #[command(
        subcommand,
        name = "convert",
        about = "Convert between different file formats",
        long_about = include_str!("strings/convert_help.txt")
    )]
    Convert(ConvertMenu),

    #[command(
        name = "dmr",
        about = "Perform DMR identification",
        long_about = include_str!("strings/dmr_help.txt")
    )]
    Dmr {
        #[clap(flatten)]
        args: dmr::DmrArgs,
    },

    #[command(subcommand, name = "dimred", about = "Dimensionality reduction tools")]
    Dimred(DimRedMenu),

    #[command(
        name = "validate",
        about = "Validate the integrity of a datasets",
        long_about = include_str!("strings/validate_help.txt")
    )]
    Validate {
        #[clap(flatten)]
        args: ValidateArgs,
    },

    #[command(
        name = "sort",
        about = "Sort a dataset",
        long_about = include_str!("strings/sort_help.txt")
    )]
    Sort {
        #[clap(flatten)]
        args: SortArgs,
    },
}

#[derive(Subcommand, Debug)]
enum ConvertMenu {
    #[command(name = "tobsx", about = "Convert a methylation report to BSX format")]
    ToBsx {
        #[clap(flatten)]
        args: ToBsxConvert,
    },
    #[command(name = "fbsx", about = "Convert a BSX dataset to a methylation report")]
    FromBsx {
        #[clap(flatten)]
        args: FromBsxConvert,
    },
    #[command(name = "r2r", about = "Convert a methylation report to another format")]
    R2R {
        #[clap(flatten)]
        args: R2RConvert,
    },
}

#[derive(Subcommand, Debug)]
enum DimRedMenu {
    #[command(name = "shrink", about = "Perform dimensionality reduction", long_about = include_str!("strings/dimred_help.txt"))]
    Shrink {
        #[clap(flatten)]
        args: DimRedArgs,
    },
    #[command(name = "merge", about = "Merge multiple dimred files")]
    Merge {
        #[clap(flatten)]
        args: MergeArgs,
    },
}

impl MainMenu {
    fn get_command(self) -> Box<dyn PipelineCommand> {
        match self {
            MainMenu::Convert(ConvertMenu::ToBsx { args }) => Box::new(args),
            MainMenu::Convert(ConvertMenu::FromBsx { args }) => Box::new(args),
            MainMenu::Convert(ConvertMenu::R2R { args }) => Box::new(args),
            MainMenu::Dmr { args } => Box::new(args),
            MainMenu::Dimred(DimRedMenu::Shrink { args }) => Box::new(args),
            MainMenu::Dimred(DimRedMenu::Merge { args }) => Box::new(args),
            MainMenu::Validate { args } => Box::new(args),
            MainMenu::Sort { args } => Box::new(args),
        }
    }
}

fn main() -> anyhow::Result<()> {
    let args: ArgsOs = wild::args_os();
    Cli::parse_from(args).command.get_command().run()?;
    Ok(())
}
