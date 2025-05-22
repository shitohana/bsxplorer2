use std::collections::HashMap;
use std::fs::File;
use std::process::exit;
use std::rc::Rc;

use bsxplorer2::io::bsx::BsxFileReader;
use clap::Args;
use console::style;
use indicatif::ProgressBar;
use itertools::Itertools;

use crate::utils::{expand_wildcards_single, init_pbar, UtilsArgs};

#[derive(Args, Debug, Clone)]
pub(crate) struct ValidateArgs {
    #[arg(
        value_parser,
        num_args=2..,
        required = true,
        help = "Paths to BSX files"
    )]
    files: Vec<String>,
}

impl ValidateArgs {
    pub fn run(
        &self,
        utils: &UtilsArgs,
    ) -> anyhow::Result<()> {
        let paths = self
            .files
            .iter()
            .map(|path| expand_wildcards_single(path))
            .flatten()
            .collect::<Vec<_>>();

        for path in paths.iter() {
            if !path.exists() {
                eprintln!("Path {} does not exist.", style(path.display()).red());
                exit(-1);
            }
            if !path.is_file() {
                eprintln!("Path {} is not a file.", style(path.display()).red());
                exit(-1);
            }
        }

        let mut readers = paths
            .iter()
            .map(|path| File::open(path).map(|f| (path.clone(), f)))
            .collect::<Result<Vec<_>, _>>()
            .unwrap_or_else(|e| {
                panic!("Error when opening files: {}", style(e.to_string()).red())
            })
            .into_iter()
            .map(|(path, file)| (path, BsxFileReader::new(file)))
            .collect::<HashMap<_, _>>();

        let batch_num = readers
            .iter()
            .map(|(path, reader)| (reader.blocks_total(), path.clone()))
            .into_group_map();

        if batch_num.len() > 1 {
            let max_batch_num = batch_num
                .iter()
                .sorted_by_key(|(_k, v)| v.len())
                .map(|(k, _v)| k)
                .max()
                .cloned()
                .unwrap();

            eprintln!("Detected {}", style("inconsistency in batch number").red());
            eprintln!("Most files have {} batches", style(max_batch_num).green());
            for (k, v) in batch_num.iter() {
                if *k == max_batch_num {
                    continue;
                }
                eprintln!("Following files have {} batches:", style(k).red());

                for path in v {
                    eprintln!("\t{}", path.to_string_lossy())
                }
                exit(-1)
            }
        }
        let common_batch_num = {
            let keys = batch_num.keys().collect_vec();
            *keys[0]
        };
        println!(
            "[{}] All files have {} batches",
            style("V").green(),
            style(common_batch_num).green()
        );

        let mut iterators = readers
            .iter_mut()
            .map(|(k, v)| (Rc::new(k.clone()), v.iter()))
            .collect::<HashMap<_, _>>();

        let progress_bar = if utils.progress {
            let progress_bar =
                init_pbar(common_batch_num).expect("Failed to initialize progress bar");
            progress_bar
        }
        else {
            ProgressBar::hidden()
        };

        let mut batch_count = 0;
        while let Some(new_batches) = iterators
            .iter_mut()
            .map(|(k, v)| {
                v.next().map(|batch_res| {
                    let batch = batch_res.unwrap_or_else(|e| {
                        eprintln!(
                            "Error when reading batch from {}",
                            k.to_string_lossy()
                        );
                        panic!("{e}")
                    });
                    (k.clone(), batch)
                })
            })
            .collect::<Option<Vec<_>>>()
        {
            let batch_size = new_batches
                .iter()
                .map(|(path, batch)| (batch.len(), path))
                .into_group_map();

            if batch_size.len() > 1 {
                let max_size = batch_size
                    .keys()
                    .max_by_key(|k| batch_size.get(k).unwrap().len())
                    .unwrap();
                eprintln!(
                    "Detected {}. Batch №{}",
                    style("inconsistency in batch size").red(),
                    batch_count
                );
                eprintln!("Most files have {} rows", style(max_size).green());

                for (size, paths) in batch_size
                    .iter()
                    .sorted_by_key(|(_s, paths)| paths.len())
                    .take(batch_size.len() - 1)
                {
                    eprintln!("Following files have {} rows:", style(size).red());
                    for path in paths {
                        eprintln!("\t{}", path.to_string_lossy())
                    }
                }

                exit(-1)
            }

            let contigs = new_batches
                .iter()
                .map(|(path, batch)| (batch.as_contig().unwrap(), path))
                .into_group_map();

            if contigs.len() > 1 {
                let common_contig = contigs
                    .keys()
                    .max_by_key(|k| contigs.get(k).unwrap().len())
                    .unwrap();
                eprintln!(
                    "Detected {}. Batch №{}",
                    style("inconsistency in batch values range").red(),
                    batch_count
                );
                eprintln!("Most files cover {} region", style(common_contig).green());

                for (contig, paths) in contigs
                    .iter()
                    .sorted_by_key(|(_s, paths)| paths.len())
                    .take(contigs.len() - 1)
                {
                    eprintln!("Following files cover {} region:", style(contig).red());
                    for path in paths {
                        eprintln!("\t{}", path.to_string_lossy())
                    }
                }
                exit(-1)
            }

            batch_count += 1;
            progress_bar.inc(1);
        }
        progress_bar.finish_and_clear();

        println!("[{}] All files have equivalent batches", style("V").green(),);
        println!("{}", style("Files are valid").green(),);
        Ok(())
    }
}
