use std::fs::File;
use std::path::PathBuf;
use std::process::exit;

use bsxplorer2::data_structs::coords::Contig;
use bsxplorer2::data_structs::enums::Strand;
use bsxplorer2::io::bsx::{BatchIndex, BsxFileReader, BsxFileWriter};
use clap::Args;
use console::style;
use indicatif::ProgressBar;
use itertools::Itertools;

use crate::utils::{init_pbar, CliIpcCompression, UtilsArgs};


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
    pub fn run(
        &self,
        utils: &UtilsArgs,
    ) -> anyhow::Result<()> {
        if !self.file.exists() {
            eprintln!(
                "File {} does not exist!",
                style(self.file.to_string_lossy()).red()
            );
            exit(-1)
        }
        if self.output.is_dir() {
            eprintln!(
                "{} is a directory!",
                style(self.output.to_string_lossy()).red()
            );
            exit(-1)
        }

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
            None,
        )?;


        let progress_bar = if utils.progress {
            let progress_bar =
                init_pbar(reader.blocks_total()).expect("Failed to initialize progress bar");
            progress_bar
        }
        else {
            ProgressBar::hidden()
        };

        for chr in chr_order {
            let batch_indices =
                index.find(&Contig::new(chr.clone().into(), 0, u32::MAX, Strand::None));

            if batch_indices.is_none() {
                eprintln!("Chromosome {} not found in file", style(chr).red());
                continue;
            }
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
