use crate::dmr_fast::{init_pbar, run_dmr, FilterArgs};
use crate::DmrContext;
use _lib::tools::dmr::dmr_binom::DmrConfig;
use _lib::tools::dmr::sure_segment::SureSegmentModelConfig;
use _lib::utils::types::Context;
use clap::Args;
use std::fs::File;
use std::path::PathBuf;

#[derive(Args, Debug, Clone)]
pub(crate) struct SegmentationArgs2 {
    #[arg(
        long,
        default_value_t = 1e-3,
        help_heading = "SEGMENTATION ARGS",
        help = "Tolerance of segment identifying after denoising step (should be low)."
    )]
    tolerance: f64,
    #[arg(
        long,
        default_value_t = 1e-2,
        help_heading = "SEGMENTATION ARGS",
        help = "Welch t-Test P-value for segments merging step. Smaller p-value -> more iterations, less falsely merged segments."
    )]
    merge_p: f64,
    #[arg(
        long,
        default_value_t = 20,
        help_heading = "SEGMENTATION ARGS",
        help = "Block size for optimal lambda approximation."
    )]
    block_size: usize,
    #[arg(
        long,
        default_value_t = 20,
        help_heading = "SEGMENTATION ARGS",
        help = "Number of optimal lambda approximation iterations for coarse grid search."
    )]
    coarse_steps: usize,
    #[arg(
        long,
        default_value_t = 30,
        help_heading = "SEGMENTATION ARGS",
        help = "Number of optimal lambda approximation iterations for fine grid search."
    )]
    refined_steps: usize,
    #[arg(
        long,
        default_value_t = 1e-3,
        help_heading = "SEGMENTATION ARGS",
        help = "Initial lambda value. Affects segments length. Must be at least 10 times less than expected methylation density range."
    )]
    lambda_min: f64,
    #[arg(
        long,
        default_value_t = 2.0,
        help_heading = "SEGMENTATION ARGS",
        help = "Maximum lambda value."
    )]
    lambda_max: f64,
    #[arg(
        short = 'D',
        long,
        default_value_t = 100,
        help_heading = "SEGMENTATION ARGS",
        help = "Maximum distance between adjacent cytosines in a segment."
    )]
    max_dist: u32,
    #[arg(
        short = 'u',
        long,
        default_value_t = 0.7,
        help_heading = "SEGMENTATION ARGS",
        help = "Segment intersection ratio threshold (intersection length / full length)."
    )]
    union_threshold: f64,
    #[arg(
        long,
        default_value_t = 1000,
        help_heading = "SEGMENTATION ARGS",
        help = "Beta-binomial distribution model fitting interations number."
    )]
    model_iters: u64,
}

pub fn run(
    filters: FilterArgs,
    segmentation: SegmentationArgs2,
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
    let context = match filters.context {
        DmrContext::CG => Context::CG,
        DmrContext::CHG => Context::CHG,
        DmrContext::CHH => Context::CHH,
    };

    let run_config = DmrConfig {
        context,
        n_missing: filters.n_missing,
        min_coverage: filters.min_coverage,
        diff_threshold: filters.diff_threshold,
        min_cpgs: filters.min_cytosines,
        segment_model: SureSegmentModelConfig {
            block_size: segmentation.block_size,
            coarse_steps: segmentation.coarse_steps,
            refined_steps: segmentation.refined_steps,
            lambda_min: segmentation.lambda_min,
            lambda_max: segmentation.lambda_max,
            seg_tolerance: segmentation.tolerance,
            max_dist: segmentation.max_dist,
            union_threshold: segmentation.union_threshold,
            max_iters: segmentation.model_iters,
            ..Default::default()
        },
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

    run_dmr(Box::new(dmr_iterator), progress_bar, filters, output);
}
