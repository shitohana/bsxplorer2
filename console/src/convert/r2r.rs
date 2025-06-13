use std::fs::File;
use std::io::IsTerminal;
use std::path::PathBuf;

use bsxplorer2::prelude::*;
use clap::{
    Args,
    ValueEnum,
};

use super::FromReportArgs;
use crate::utils::{
    init_hidden,
    init_spinner,
};

#[derive(Debug, Clone, ValueEnum, Eq, PartialEq)]
pub enum ConvertReportType {
    Bsx,
    Bismark,
    CgMap,
    BedGraph,
    Coverage,
}

#[derive(Debug, Clone, Args)]
pub struct R2RConvert {
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

    #[clap(flatten)]
    from_report: FromReportArgs,
}

impl R2RConvert {
    pub fn run(&self) -> anyhow::Result<()> {
        let pbar = if std::io::stdin().is_terminal() {
            init_spinner()?
        }
        else {
            init_hidden()?
        };

        let mut report_reader_builder = ReportReaderBuilder::default()
            .with_batch_size(self.from_report.batch_size)
            .with_low_memory(self.from_report.low_memory)
            .with_report_type(self.from_report.from_type)
            .with_compression(Some(self.from_report.from_compression));

        if let Some(fasta_path) = &self.from_report.fasta_path {
            report_reader_builder =
                report_reader_builder.with_fasta_path(Some(fasta_path.clone()));
        }
        let report_reader = report_reader_builder.build(self.input.clone())?;

        let output_handle = File::create(self.output.clone())?;
        let mut writer = ReportWriter::try_new(
            output_handle,
            self.to_type,
            bsxplorer2::utils::n_threads(),
            self.to_compression.clone(),
            self.compression_level,
        )?;

        for batch in report_reader {
            let batch = batch?;
            let cur_pos: GenomicPosition = batch.last_genomic_pos().unwrap_or_default();
            pbar.set_message(format!("Reading: {}", cur_pos));
            pbar.inc(1);

            writer.write_batch(batch)?;
        }
        writer.finish()?;
        Ok(())
    }
}
