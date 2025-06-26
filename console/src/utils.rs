use std::error::Error;
use std::fmt::Display;
use std::fs::File;
use std::path::{
    Path,
    PathBuf,
};

use anyhow::{
    ensure,
    Context,
    Result,
};
use bsxplorer2::prelude::BsxFileReader;
use clap::ValueEnum;
use indicatif::{
    ProgressBar,
    ProgressStyle,
};
use log::{
    debug,
    info,
};
use polars::prelude::IpcCompression;
use spipe::spipe;

#[derive(Debug)]
#[allow(unused)]
pub(crate) enum CliError {
    PathExpand(String),
    IsDirectory(PathBuf),
    CantOpenFile(PathBuf, std::io::Error),
    CantOpenBsx(anyhow::Error),
}

impl Display for CliError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            CliError::PathExpand(path) => write!(f, "Can't expand path {}", path),
            CliError::IsDirectory(path_buf) => {
                write!(f, "{} is a directory!", path_buf.display())
            },
            CliError::CantOpenFile(path_buf, error) => {
                write!(f, "Failed to open file {}: {}", path_buf.display(), error)
            },
            CliError::CantOpenBsx(polars_error) => {
                write!(
                    f,
                    "Failed to open bsx file. Ensure the file is valid. Error: {}",
                    polars_error
                )
            },
        }
    }
}

impl Error for CliError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::CantOpenFile(_, e) => Some(e),
            _ => None,
        }
    }

    fn description(&self) -> &str {
        "description() is deprecated; use Display"
    }

    fn cause(&self) -> Option<&dyn Error> {
        self.source()
    }
}

pub fn init_progress(
    max: Option<usize>
) -> std::result::Result<ProgressBar, anyhow::Error> {
    use std::io::IsTerminal;
    if std::io::stdin().is_terminal() {
        if let Some(max) = max {
            spipe!(
                ProgressBar::new(max as u64)
                    =>$ .set_style(
                        ProgressStyle::default_bar()
                            .template(
                                "{spinner:.green} [{elapsed_precise}, ETA: {eta}] \
                                [{bar:40.cyan/blue}] {pos:>5.green}/{len:5} {msg}",
                            )?
                            .progress_chars("#>-")
                    )
                    => Ok
            )
        }
        else {
            spipe!(
                ProgressBar::new_spinner()
                    =>$ .set_style(ProgressStyle::default_spinner()
                        .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                        .template("{spinner} {msg}")?)
                    => Ok
            )
        }
    }
    else {
        Ok(ProgressBar::hidden())
    }
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

pub fn expand_wildcard<S: AsRef<str>>(path: S) -> Result<Vec<PathBuf>, CliError> {
    if path.as_ref().contains('*') || path.as_ref().contains('?') {
        Ok(vec![PathBuf::from(path.as_ref())])
    }
    else {
        Ok(glob::glob(path.as_ref())
            .map_err(|_| CliError::PathExpand(path.as_ref().to_string()))?
            .map(|p| p.expect("Failed to access path"))
            .collect())
    }
}

#[macro_export]
macro_rules! log_on_success {
    ($op:expr, $log_level:ident, $err_msg:literal$(, $subst:expr)*) => {
        if $op.is_ok() {
            log::$log_level!($err_msg$(, $subst)*);
        }
    };
}

/// For every path:
///     1. Expand wildcards
///     2. Ensure path exists
///     3. Ensure path is not empty
///     4. Ensure path is file
pub fn check_validate_paths<R: AsRef<str>>(
    paths: &[R]
) -> anyhow::Result<Vec<PathBuf>> {
    let paths = paths
        .iter()
        .map(expand_wildcard)
        .collect::<anyhow::Result<Vec<_>, _>>()
        .context("Failed to expand wildcards in paths")?
        .concat();
    ensure!(paths.len() > 0, "Paths cannot be empty");

    paths.iter().try_for_each(|x| { validate_input(x).map(|_| {}) })?;
    info!("Paths {:?} valid", paths);
    Ok(paths)
}

pub fn validate_input<P: AsRef<Path>>(path: P) -> anyhow::Result<P> {
    ensure!(path.as_ref().exists(), "Path {} does not exist", path.as_ref().display());
    ensure!(path.as_ref().is_file(), "Path {} is not a file", path.as_ref().display());
    debug!("Path {} valid", path.as_ref().display());
    Ok(path)
}

pub fn validate_output<P: AsRef<Path>>(path: P) -> anyhow::Result<P> {
    ensure!(!path.as_ref().is_dir(), "Path {} is a directory", path.as_ref().display());
    if path.as_ref().exists() {
        info!("Overwriting {}", path.as_ref().display())
    }
    Ok(path)
}

/// For every path
///     1. Open file
///     2. Open [BsxFileReader]
pub fn init_readers<I, S>(iter: I) -> Result<Vec<BsxFileReader>>
where
    I: IntoIterator<Item = S>,
    S: AsRef<Path> + Clone, {
    iter.into_iter()
        .map(|p| spipe!{
            p
                =>  .as_ref
                =>  File::open
                =>  .context(p.as_ref().display().to_string())
                =>& BsxFileReader::try_new
                =>  .context(p.as_ref().display().to_string())
                =># |res: &Result<_>| log_on_success!(
                    &res,
                    debug,
                    "Opened {} for {}", stringify!(BsxFileReader), p.as_ref().display())

        })
        .collect::<Result<Vec<_>, _>>()
}

// We need this because IpcCompression is defined in polars crate
// and we cannot just add FromStr to it
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
