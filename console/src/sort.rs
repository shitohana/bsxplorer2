use std::fs::File;
use std::path::PathBuf;

use bsxplorer2::prelude::*;
use clap::Args;
use indicatif::ProgressBar;
use itertools::Itertools;

use crate::{assert_or_exit, init_progress, utils::{
    init_pbar,
    CliIpcCompression,
}};


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
        help = "Order of the chromosomes"
    )]
    order: Option<Vec<String>>,
}

impl SortArgs {
    pub fn run(&self) -> anyhow::Result<()> {
        assert_or_exit!(self.file.exists(), "File {} does not exist!", self.file.display());
        assert_or_exit!(!self.output.is_dir(), "{} is a directory!", self.output.display());

        let mut reader = BsxFileReader::try_new(File::open(self.file.clone())?)?;

        let index = BatchIndex::from_reader(&mut reader)?;
        let chr_order = self.order.clone().unwrap_or(
            index
                .get_chr_order()
                .iter()
                .sorted()
                .map(|v| v.to_string())
                .collect_vec(),
        );

        let mut writer = BsxFileWriter::try_new(
            File::create(self.output.clone())?,
            chr_order.clone(),
            self.to_compression.into(),
        )?;


        let progress_bar = init_progress!(reader.blocks_total())?;

        for chr in chr_order {
            let batch_indices =
                index.find(&Contig::new(chr.clone().into(), 0, u32::MAX, Strand::None));

            assert_or_exit!(batch_indices.is_some(), "Chromosome {} not found in file", chr);
            let batches = batch_indices.unwrap();
            for batch_idx in batches {
                let batch = reader.get_batch(batch_idx).unwrap()?;
                writer.write_batch(batch)?;
                progress_bar.inc(1);
            }
        }

        progress_bar.finish_and_clear();
        writer.close()?;
        Ok(())
    }
}
