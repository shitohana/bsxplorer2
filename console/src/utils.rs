use std::path::PathBuf;

use clap::ValueEnum;
use glob::glob;
use indicatif::{
    ProgressBar,
    ProgressStyle,
};
use polars::prelude::IpcCompression;

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

pub fn init_spinner() -> anyhow::Result<ProgressBar> {
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
            .template("{spinner} {msg}")
            .expect("Failed to set spinner template"),
    );
    Ok(spinner)
}

#[macro_export]
macro_rules! init_progress {
    () => { init_progress!(@inner init_spinner())  };
    ($total:expr) => { init_progress!(@inner init_pbar($total)) };
    (@inner $init_with:expr ) => {{
        use std::io::IsTerminal;
        if std::io::stdin().is_terminal() {
            $init_with
        }
        else {
            Ok(ProgressBar::hidden())
        }
    }}
}

#[macro_export]
macro_rules! exit_with_msg {
    ($fmt:literal$(, $value:expr)*) => {
        eprintln!($fmt$(, console::style($value).red())*);
        std::process::exit(-1);
    };
}

#[macro_export]
macro_rules! assert_or_exit {
    ($expr:expr, $fmt:literal$(, $value:expr)*) => {
        if !$expr { crate::exit_with_msg!($fmt $(, $value)*); }
    };
}


pub fn init_hidden() -> anyhow::Result<ProgressBar> {
    Ok(ProgressBar::hidden())
}

pub fn expand_wildcards(paths: Vec<String>) -> Vec<PathBuf> {
    paths
        .iter()
        .map(|path| expand_wildcards_single(path))
        .flatten()
        .collect()
}

pub fn expand_wildcards_single(path: &String) -> Vec<PathBuf> {
    let mut expanded_paths = vec![];
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
    expanded_paths
}

#[derive(Debug, Clone, ValueEnum, Eq, PartialEq, Copy)]
pub enum CliIpcCompression {
    LZ4,
    ZSTD,
    None,
}

impl From<CliIpcCompression> for Option<IpcCompression> {
    fn from(compression: CliIpcCompression) -> Self {
        match compression {
            CliIpcCompression::LZ4 => Some(IpcCompression::LZ4),
            CliIpcCompression::ZSTD => Some(IpcCompression::ZSTD),
            CliIpcCompression::None => None,
        }
    }
}
