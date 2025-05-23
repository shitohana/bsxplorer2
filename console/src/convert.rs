use std::fs::File;
use std::path::PathBuf;
use std::process::exit;

use bsxplorer2::data_structs::coords::GenomicPosition;
use bsxplorer2::io::bsx::{BsxFileReader, BsxFileWriter};
use bsxplorer2::io::compression::Compression;
use bsxplorer2::io::report::{ReportReaderBuilder, ReportType, ReportWriter};
use clap::{Args, ValueEnum};
use console::style;

use crate::utils::{init_hidden,
                   init_pbar,
                   init_spinner,
                   CliIpcCompression,
                   UtilsArgs};

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

    #[clap(short='f', long = "from", required = true, value_enum, default_value_t = ReportType::Bismark)]
    from_type: ReportType,

    #[clap(short='F', long = "from-compression", required = true, value_enum, default_value_t = Compression::None)]
    from_compression: Compression,

    #[clap(short='T', long = "to-compression", required = true, value_enum, default_value_t = CliIpcCompression::None)]
    to_compression: CliIpcCompression,

    #[arg(
        long,
        default_value_t = false,
        help = "Use less RAM, but elongate computation.",
        help_heading = "REPORT ARGS"
    )]
    low_memory: bool,

    #[arg(
        short,
        long = "chunk",
        default_value_t = 10_000,
        help = "Number of rows in the output batches.",
        help_heading = "REPORT ARGS"
    )]
    chunk_size: usize,

    #[arg(
        long = "fa",
        help = "Path to the reference sequence file.",
        help_heading = "REPORT ARGS"
    )]
    fasta_path: Option<PathBuf>,

    #[arg(
        long = "fai",
        help = "Path to the fasta index.",
        help_heading = "REPORT ARGS"
    )]
    fai_path: Option<PathBuf>,

    #[arg(long, default_value_t = 2 << 20, help = "Size of raw batches.")]
    batch_size: usize,
}

impl ToBsxConvert {
    pub fn run(
        &self,
        utils_args: &UtilsArgs,
    ) -> anyhow::Result<()> {
        if matches!(self.from_type, ReportType::BedGraph | ReportType::Coverage)
            && !(self.fasta_path.is_some() && self.fai_path.is_some())
        {
            eprintln!(
                "For {:?} report type conversion, both {}",
                self.from_type,
                style("fasta_path and fai_path must be set").red()
            )
        }

        let pbar = if utils_args.progress {
            init_spinner()?
        }
        else {
            init_hidden()?
        };

        let mut report_reader_builder = ReportReaderBuilder::default()
            .with_n_threads(Some(utils_args.threads))
            .with_batch_size(self.batch_size)
            .with_chunk_size(self.chunk_size)
            .with_low_memory(self.low_memory)
            .with_report_type(self.from_type)
            .with_compression(Some(self.from_compression));

        if let Some(fasta_path) = &self.fasta_path {
            report_reader_builder =
                report_reader_builder.with_fasta_path(Some(fasta_path.clone()));
        }
        if let Some(fai_path) = &self.fai_path {
            report_reader_builder =
                report_reader_builder.with_fai_path(Some(fai_path.clone()));
        }

        let report_reader = report_reader_builder.build(self.input.clone())?;

        let sink = File::create(self.output.clone())?;
        let mut bsx_writer = if let Some(fai_path) = &self.fai_path {
            BsxFileWriter::try_from_sink_and_fai(
                sink,
                fai_path.clone(),
                self.to_compression.into(),
                None,
            )?
        }
        else if let Some(fasta_path) = &self.fasta_path {
            BsxFileWriter::try_from_sink_and_fasta(
                sink,
                fasta_path.clone(),
                self.to_compression.into(),
                None,
            )?
        }
        else {
            eprintln!(
                "{}",
                style("Either fasta_path or fai_path must be set").red()
            );
            exit(0x0100)
        };

        for batch in report_reader {
            let batch = batch?;
            let cur_pos: GenomicPosition = batch.last_genomic_pos().unwrap_or_default();
            pbar.set_message(format!("Reading: {}", cur_pos));
            pbar.inc(1);

            bsx_writer.write_batch(batch)?;
        }
        bsx_writer.close()?;
        Ok(())
    }
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
    pub fn run(
        &self,
        utils_args: &UtilsArgs,
    ) -> anyhow::Result<()> {
        let input_handle = File::open(self.input.clone())?;
        let bsx_reader = BsxFileReader::try_new(input_handle)?;

        let output_handle = File::create(self.output.clone())?;
        let mut writer = ReportWriter::try_new(
            output_handle,
            self.to_type,
            utils_args.threads,
            self.to_compression.clone(),
            self.compression_level,
        )?;

        let pbar = if utils_args.progress {
            init_pbar(bsx_reader.blocks_total())?
        }
        else {
            init_hidden()?
        };

        for batch in bsx_reader.into_iter() {
            writer.write_batch(batch?)?;
            pbar.inc(1);
        }
        writer.finish()?;
        Ok(())
    }
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

    #[clap(short='f', long = "from", required = true, value_enum, default_value_t = ReportType::Bismark)]
    from_type: ReportType,

    #[clap(short='t', long = "to", required = true, value_enum, default_value_t = ReportType::Bismark)]
    to_type: ReportType,

    #[clap(short='F', long = "from-compression", required = true, value_enum, default_value_t = Compression::None)]
    from_compression: Compression,

    #[clap(short='T', long = "to-compression", required = true, value_enum, default_value_t = Compression::None)]
    to_compression: Compression,

    #[clap(short='L', long = "level", required = false, value_enum, default_value = None)]
    compression_level: Option<u32>,

    #[arg(
        long,
        default_value_t = false,
        help = "Use less RAM, but elongate computation.",
        help_heading = "REPORT ARGS"
    )]
    low_memory: bool,

    #[arg(
        long = "fa",
        help = "Path to the reference sequence file.",
        help_heading = "REPORT ARGS"
    )]
    fasta_path: Option<PathBuf>,

    #[arg(long, default_value_t = 2 << 20, help = "Size of raw batches.")]
    batch_size: usize,
}

impl R2RConvert {
    pub fn run(
        &self,
        utils_args: &UtilsArgs,
    ) -> anyhow::Result<()> {
        let pbar = if utils_args.progress {
            init_spinner()?
        }
        else {
            init_hidden()?
        };

        let mut report_reader_builder = ReportReaderBuilder::default()
            .with_n_threads(Some(utils_args.threads))
            .with_batch_size(self.batch_size)
            .with_low_memory(self.low_memory)
            .with_report_type(self.from_type)
            .with_compression(Some(self.from_compression));

        if let Some(fasta_path) = &self.fasta_path {
            report_reader_builder =
                report_reader_builder.with_fasta_path(Some(fasta_path.clone()));
        }
        let report_reader = report_reader_builder.build(self.input.clone())?;

        let output_handle = File::create(self.output.clone())?;
        let mut writer = ReportWriter::try_new(
            output_handle,
            self.to_type,
            utils_args.threads,
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
