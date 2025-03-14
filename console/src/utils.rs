/// ***********************************************************************
/// *****
/// * Copyright (c) 2025
/// The Prosperity Public License 3.0.0
///
/// Contributor: [shitohana](https://github.com/shitohana)
///
/// Source Code: https://github.com/shitohana/BSXplorer
/// ***********************************************************************
/// ****

/// ***********************************************************************
/// *****
/// * Copyright (c) 2025
/// ***********************************************************************
/// ****
use std::path::PathBuf;

use bsxplorer2::exports::{anyhow, rayon};
use bsxplorer2::utils::types::Context;
use clap::{Args, ValueEnum};
use glob::glob;
use indicatif::{ProgressBar, ProgressStyle};

pub fn init_pbar(total: usize) -> anyhow::Result<ProgressBar> {
    let progress_bar = ProgressBar::new(total as u64);
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}, ETA: {eta}] \
                 [{bar:40.cyan/blue}] {pos:>5.green}/{len:5} {msg}",
            )?
            .progress_chars("#>-"),
    );
    progress_bar.set_message("Processing...");
    Ok(progress_bar)
}

#[derive(Debug, Clone, ValueEnum)]
pub(crate) enum DmrContext {
    CG,
    CHG,
    CHH,
}

impl DmrContext {
    pub fn tobsxplorer2(&self) -> Context {
        match self {
            DmrContext::CG => Context::CG,
            DmrContext::CHG => Context::CHG,
            DmrContext::CHH => Context::CHH,
        }
    }
}

pub(crate) fn expand_wildcards(paths: Vec<String>) -> Vec<PathBuf> {
    let mut expanded_paths = Vec::new();

    for path in paths {
        if path.contains('*') || path.contains('?') {
            // Expand wildcard using glob
            match glob(&path) {
                Ok(matches) => {
                    for entry in matches.filter_map(Result::ok) {
                        expanded_paths.push(entry);
                    }
                },
                Err(e) => {
                    eprintln!("Error processing wildcard '{}': {}", path, e)
                },
            }
        }
        else {
            // If not a wildcard, push the path as-is
            expanded_paths.push(PathBuf::from(path));
        }
    }

    expanded_paths
}

#[derive(Args, Debug)]
pub(crate) struct UtilsArgs {
    #[arg(
        long,
        required = false,
        default_value_t = true,
        help = "Display progress bar (Disable if you need clean pipeline \
                logs)."
    )]
    pub(crate) progress: bool,
    #[arg(
        long,
        required = false,
        default_value_t = 1,
        help = "Number of threads to use."
    )]
    pub(crate) threads:  usize,
    #[arg(
        long,
        required = false,
        default_value_t = false,
        help = "Verbose output."
    )]
    pub(crate) verbose:  bool,
}

#[inline]
pub(crate) fn init_logger(logger: bool) -> anyhow::Result<()> {
    if logger {
        Ok(bsxplorer2::exports::pretty_env_logger::try_init()?)
    }
    else {
        Ok(())
    }
}

#[inline]
pub(crate) fn init_rayon_threads(threads: usize) -> anyhow::Result<()> {
    Ok(rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?)
}
