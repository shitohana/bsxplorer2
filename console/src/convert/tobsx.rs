use std::fs::File;
use std::path::PathBuf;

use anyhow::{anyhow, ensure};
use bsxplorer2::prelude::*;
use clap::{
    Args,
    ValueEnum,
};
use spipe::spipe;

use super::FromReportArgs;
use crate::utils::{
    init_progress, validate_output, CliIpcCompression
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
pub struct ToBsxConvert {
    #[arg(help = "Path of the input file.")]
    input: PathBuf,

    #[arg(
        short = 'o',
        long,
        required = true,
        help = "Path for the generated output file."
    )]
    output: PathBuf,

    #[clap(short='T', long = "to-compression", required = true, value_enum, default_value_t = CliIpcCompression::None)]
    to_compression: CliIpcCompression,

    #[clap(
        short,
        long,
        value_enum,
        default_value_t = Context::CG,
        help = "Select cytosine methylation context."
    )]
    context: Context,

    #[arg(
        short,
        long = "chunk",
        default_value_t = 10_000,
        help = "Number of rows in the output batches."
    )]
    chunk_size: usize,

    #[arg(long = "fai", help = "Path to the fasta index.")]
    fai_path: Option<PathBuf>,

    #[clap(flatten)]
    from_report: FromReportArgs,
}

impl PipelineCommand for ToBsxConvert {
    fn run(&self) -> anyhow::Result<()> {
        {
            let need_align = self.from_report.from_type.need_align();
            let fasta_defined = self.from_report.fasta_path.is_some();
            let fai_defined = self.fai_path.is_some();

            ensure!(
                !need_align || (fasta_defined && fai_defined),
                "For {:?} report type conversion, both fasta path and fai path must be set",
                self.from_report.from_type,
            )
        }

        let mut report_reader_builder = ReportReaderBuilder::default()
            .with_batch_size(self.from_report.batch_size)
            .with_chunk_size(self.chunk_size)
            .with_low_memory(self.from_report.low_memory)
            .with_report_type(self.from_report.from_type)
            .with_compression(Some(self.from_report.from_compression));

        if let Some(fasta_path) = &self.from_report.fasta_path {
            report_reader_builder =
                report_reader_builder.with_fasta_path(Some(fasta_path.clone()));
        }
        if let Some(fai_path) = &self.fai_path {
            report_reader_builder =
                report_reader_builder.with_fai_path(Some(fai_path.clone()));
        }

        let report_reader = report_reader_builder.build(self.input.clone())?;

        let mut bsx_writer = if let Some(fai_path) = &self.fai_path {
            spipe!(
                &self.output =>
                validate_output =>?
                File::create =>?
                BsxFileWriter::try_from_sink_and_fai(
                    fai_path.clone(),
                    self.to_compression.into()
                ) =>? ...
            )
        }
        else if let Some(fasta_path) = &self.from_report.fasta_path {
            spipe!(
                &self.output =>
                validate_output =>?
                File::create =>?
                BsxFileWriter::try_from_sink_and_fasta(
                    fasta_path.clone(),
                    self.to_compression.into()
                ) =>? ...
            )
        }
        else {
            return Err(anyhow!("Either fasta_path or fai_path must be set"))
        };

        let pbar = init_progress(None)?;
        for batch in pbar.wrap_iter(report_reader) {
            let batch = batch?;
            let cur_pos: GenomicPosition = batch.last_genomic_pos().unwrap_or_default();

            bsx_writer.write_batch(batch)?;
            pbar.set_message(format!("Reading: {}", cur_pos));
        }

        bsx_writer.close()?;
        Ok(())
    }
}
