use std::fs::File;
use std::io::IsTerminal;
use std::path::PathBuf;

use bsxplorer2::prelude::*;
use clap::{
    Args,
    ValueEnum,
};
use indicatif::ProgressBar;

use crate::utils::init_pbar;

#[derive(Debug, Clone, ValueEnum, Eq, PartialEq)]
pub enum ConvertReportType {
    Bsx,
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
}

#[derive(Debug, Clone, Args)]
pub struct FromBsxConvert {
    #[arg(help = "Path of the input file.")]
    input: PathBuf,

    #[arg(
        short = 'o',
        long,
        required = true,
        help = "Path for the generated output file."
    )]
    output: PathBuf,

    #[clap(short='t', long = "to", required = true, value_enum, default_value_t = ReportType::Bismark)]
    to_type: ReportType,

    #[clap(short='T', long = "to-compression", required = true, value_enum, default_value_t = Compression::None)]
    to_compression: Compression,

    #[clap(short='L', long = "level", required = false, value_enum, default_value = None)]
    compression_level: Option<u32>,
}

impl FromBsxConvert {
    pub fn run(&self) -> anyhow::Result<()> {
        let input_handle = File::open(self.input.clone())?;
        let bsx_reader = BsxFileReader::try_new(input_handle)?;

        let output_handle = File::create(self.output.clone())?;
        let mut writer = ReportWriter::try_new(
            output_handle,
            self.to_type,
            bsxplorer2::utils::n_threads(),
            self.to_compression.clone(),
            self.compression_level,
        )?;

        let pbar = if std::io::stdin().is_terminal() {
            init_pbar(0)?
        }
        else {
            ProgressBar::hidden()
        };

        for batch in bsx_reader.into_iter() {
            writer.write_batch(batch?)?;
            pbar.inc(1);
        }
        writer.finish()?;
        Ok(())
    }
}
