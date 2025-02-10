#![feature(path_add_extension)]

use clap::ValueEnum;
use glob::glob;
use std::path::PathBuf;

pub mod convert;
pub mod dmr_binom;
pub mod dmr_fast;
pub mod stats;

#[derive(Debug, Clone, ValueEnum)]
pub(crate) enum DmrContext {
    CG,
    CHG,
    CHH,
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
                }
                Err(e) => eprintln!("Error processing wildcard '{}': {}", path, e),
            }
        } else {
            // If not a wildcard, push the path as-is
            expanded_paths.push(PathBuf::from(path));
        }
    }

    expanded_paths
}
