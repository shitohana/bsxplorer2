mod convert;
mod dimred;
mod dmr;
mod sort;
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
use sort::SortArgs;
use validate::ValidateArgs;
use wild::ArgsOs;

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

    #[command(
        name = "dimred",
        about = "Perform dimensionality reduction",
        long_about = include_str!("strings/dimred_help.txt")
    )]
    Dimred {
        #[clap(flatten)]
        args: dimred::DimRedArgs,
    },

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

fn main() -> anyhow::Result<()> {
    let args: ArgsOs = wild::args_os();
    let cli = Cli::parse_from(args);

    match cli.command {
        MainMenu::Convert(ConvertMenu::FromBsx { args }) => {
            args.run()?;
        },
        MainMenu::Convert(ConvertMenu::ToBsx { args }) => {
            args.run()?;
        },
        MainMenu::Convert(ConvertMenu::R2R { args }) => {
            args.run()?;
        },
        MainMenu::Dimred { args } => {
            args.run()?;
        },
        MainMenu::Dmr { args } => {
            args.run()?;
        },
        MainMenu::Validate { args } => args.run()?,
        MainMenu::Sort { args } => args.run()?,
    }
    Ok(())
}
