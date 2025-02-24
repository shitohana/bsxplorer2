use crate::dmr_fast::{init_pbar, run_dmr, FilterArgs};
use crate::{expand_wildcards, DmrContext};
use clap::Args;
use std::fs::File;
use std::path::PathBuf;

#[derive(Args, Debug, Clone)]
pub(crate) struct SegmentationArgs {
    #[arg(
        short = 'D',
        long,
        default_value_t = 100,
        help_heading = "SEGMENTATION ARGS",
        help = "Maximum distance between adjacent cytosines in a segment."
    )]
    max_dist: u32,
    #[arg(
        short = 'T',
        long,
        default_value_t = 0.5,
        help_heading = "SEGMENTATION ARGS",
        help = "Trend threshold"
    )]
    trend_threshold: f64,
    #[arg(
        short = 'V',
        long,
        default_value_t = 0.7,
        help_heading = "SEGMENTATION ARGS",
        help = "Value of a valley filter"
    )]
    valley: f64,
}

pub fn run(
    filters: FilterArgs,
    segmentation: SegmentationArgs,
    group_a: Vec<String>,
    group_b: Vec<String>,
    output: PathBuf,
    force: bool,
    threads: usize,
    progress: bool,
) {
    let (sample_paths, sample_labels) = if let Ok(result) =
        crate::dmr_fast::init(group_a, group_b, output.clone(), force, threads)
    {
        result
    } else {
        return;
    };
    let context = filters.context.to_lib();

    let run_config = _lib::tools::metilene::DmrConfig {
        context,
        n_missing: filters.n_missing,
        min_coverage: filters.min_coverage,
        diff_threshold: filters.diff_threshold,
        min_cpgs: filters.min_cytosines,
        valley: segmentation.valley,
        trend_threshold: segmentation.trend_threshold,
        max_dist: segmentation.max_dist as u64,
    };

    let files = sample_paths
        .iter()
        .map(|path| File::open(path).expect(format!("Could not open file {}", path).as_str()));
    let labeled = sample_labels.into_iter().zip(files).collect::<Vec<_>>();

    let dmr_iterator = run_config
        .try_finish(labeled)
        .expect("Failed to initialize DMR iterator");

    let progress_bar = if progress {
        let progress_bar =
            init_pbar(dmr_iterator.blocks_total()).expect("Failed to initialize progress bar");

        Some(progress_bar)
    } else {
        None
    };

    let mut filters = filters;
    filters.min_cytosines = 0;
    run_dmr(Box::new(dmr_iterator), progress_bar, filters, output);
}
