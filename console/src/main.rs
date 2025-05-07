pub mod convert;
pub mod dmr;
pub mod utils;

use clap::{Parser, Subcommand};
use convert::{FromBsxConvert, R2RConvert, ToBsxConvert};
use utils::UtilsArgs;
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
    #[command(subcommand, name = "convert")]
    Convert(ConvertMenu),

    Dmr {
        #[clap(flatten)]
        utils: UtilsArgs,
        #[clap(flatten)]
        args:  dmr::DmrArgs,
    },
}

#[derive(Subcommand, Debug)]
enum ConvertMenu {
    #[command(name = "to-bsx")]
    ToBsx {
        #[clap(flatten)]
        utils: UtilsArgs,
        #[clap(flatten)]
        args:  ToBsxConvert,
    },
    #[command(name = "from-bsx")]
    FromBsx {
        #[clap(flatten)]
        utils: UtilsArgs,
        #[clap(flatten)]
        args:  FromBsxConvert,
    },
    #[command(name = "r2r")]
    R2R {
        #[clap(flatten)]
        utils: UtilsArgs,
        #[clap(flatten)]
        args:  R2RConvert,
    },
}

fn main() -> anyhow::Result<()> {
    let args: ArgsOs = wild::args_os();
    let cli = Cli::parse_from(args);

    match cli.command {
        MainMenu::Convert(ConvertMenu::FromBsx { utils, args }) => {
            utils.setup()?;
            args.run(&utils)?;
        },
        MainMenu::Convert(ConvertMenu::ToBsx { utils, args }) => {
            utils.setup()?;
            args.run(&utils)?;
        },
        MainMenu::Convert(ConvertMenu::R2R { utils, args }) => {
            utils.setup()?;
            args.run(&utils)?;
        },
        MainMenu::Dmr { utils, args } => {
            utils.setup()?;
            args.run(&utils)?;
        },
    }
    Ok(())
}
