use std::fs::File;
use std::path::PathBuf;

use anyhow::anyhow;
use bsxplorer2::prelude::*;
use clap::Args;
use itertools::Itertools;
use log::info;
use spipe::spipe;

use crate::utils::{
    init_progress, validate_input, validate_output, CliIpcCompression
};
use crate::PipelineCommand;

#[derive(Args, Debug, Clone)]
pub(crate) struct SortArgs {
    #[arg(required = true, help = "Path to BSX file")]
    file: PathBuf,

    #[arg(required = true, help = "Path to sorted file", short, long)]
    output: PathBuf,

    #[clap(short='T', long = "to-compression", required = true, value_enum, default_value_t = CliIpcCompression::None)]
    to_compression: CliIpcCompression,

    #[arg(
        short = 'O', long,
        num_args=2..,
        help = "Order of the chromosomes",
    )]
    order: Option<Vec<String>>,
}

impl PipelineCommand for SortArgs {
    fn run(&self) -> anyhow::Result<()> {
        validate_input(&self.file)?;

        let mut reader = spipe!(
            &self.file =>
            validate_input =>?
            File::open =>?
            BsxFileReader::try_new =>?
            ...
        );
        info!("Created reader at {}", self.file.display());

        let index = BatchIndex::from_reader(&mut reader)?;
        let chr_order = self.order.clone().unwrap_or(
            index
                .get_chr_order()
                .iter()
                .sorted()
                .map(|v| v.to_string())
                .collect_vec(),
        );

        let mut writer = spipe!(
            &self.output =>
            validate_output =>?
            File::create =>?
            BsxFileWriter::try_new(&chr_order, self.to_compression.into()) =>?
            ...
        );
        info!("Created writer at {}", self.output.display());

        let progress_bar = init_progress(Some(reader.blocks_total()))?;

        for chr in chr_order {
            let batches = index.chr_indices(&chr)
                .ok_or(anyhow!("Chromosome {} not present in file", chr))?;

            for batch_idx in batches {
                spipe!(
                    reader.get_batch(batch_idx) =>
                    .unwrap =>?
                    |x| writer.write_batch(x) =>? ...
                );
                progress_bar.inc(1);
            }
        }

        progress_bar.finish_and_clear();
        writer.close()?;
        Ok(())
    }
}
