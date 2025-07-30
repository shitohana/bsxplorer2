use std::fs::File;
use std::path::PathBuf;

use bsxplorer2::prelude::*;
use clap::{
    Args,
    ValueEnum,
};
use spipe::spipe;

use crate::utils::{
    init_progress,
    validate_input,
    validate_output,
};
use crate::PipelineCommand;

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

    #[clap(short='t', long = "to", required = true, default_value_t = ReportType::Bismark)]
    to_type: ReportType,

    #[clap(short='T', long = "to-compression", required = true, default_value_t = Compression::None)]
    to_compression: Compression,

    #[clap(short='L', long = "level", required = false, value_enum, default_value = None)]
    compression_level: Option<u32>,
}

impl PipelineCommand for FromBsxConvert {
    fn run(&self) -> anyhow::Result<()> {
        let reader = spipe!(
            &self.input =>
            validate_input =>?
            File::open =>?
            BsxFileReader::try_new =>?
            ...
        );

        let mut writer = spipe!(
            &self.output =>
            validate_output =>?
            File::create =>?
            ReportWriter::try_new(
                self.to_type,
                bsxplorer2::utils::n_threads(),
                self.to_compression,
                self.compression_level
            ) =>? ...
        );

        let pbar = init_progress(Some(reader.blocks_total()))?;

        for batch in pbar.wrap_iter(reader.into_iter()) {
            writer.write_batch(batch?)?;
        }
        writer.finish()?;
        Ok(())
    }
}
