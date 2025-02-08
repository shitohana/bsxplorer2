use crate::{expand_wildcards, DmrContext};
use _lib::tools::dmr2::{DmrConfig, SegmentModelConfig};
use _lib::utils::types::Context;
use clap::Args;
use console::style;
use dialoguer::Confirm;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::iter::repeat_n;
use std::path::PathBuf;

#[derive(Args, Debug, Clone)]
pub(crate) struct FilterArgs2 {
    #[clap(short='c', long, value_enum, default_value_t = DmrContext::CG, help_heading="FILTER ARGS", help = "Select cytosine methylation context.")]
    context: DmrContext,
    #[arg(
        short = 'm',
        long,
        default_value_t = 1,
        help_heading = "FILTER ARGS",
        help = "Set missing values threshold. Cytosines with no data in > n_samples will be discarded."
    )]
    n_missing: usize,
    #[arg(
        short = 'v',
        long,
        default_value_t = 5,
        help_heading = "FILTER ARGS",
        help = "Cytosines with coverage below threshold will be discarded."
    )]
    min_coverage: i16,
    #[arg(
        short = 'q',
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "Q-value for DMR identification."
    )]
    q_value: f64,
    #[arg(
        short = 'C',
        long,
        default_value_t = 10,
        help_heading = "FILTER ARGS",
        help = "DMRs with cytosine count below threshold will be discarded."
    )]
    min_cytosines: usize,
    #[arg(
        short = 'd',
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "DMRs with difference between samples below threshold will be discarded."
    )]
    diff_threshold: f64,
}

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
    filters: FilterArgs2,
    segmentation: SegmentationArgs2,
    group_a: Vec<String>,
    group_b: Vec<String>,
    output: PathBuf,
    force: bool,
    progress: bool,
) {
    let a_paths = expand_wildcards(group_a);
    let b_paths = expand_wildcards(group_b);

    if !force {
        let prompt = format!("Do you want to proceed with the following paths?\n\nGroup A: {:?}\nGroup B: {:?}\nOutput:{:?}", a_paths, b_paths, output);
        let confirmed = Confirm::new()
            .with_prompt(prompt)
            .default(true)
            .interact()
            .unwrap_or(false);

        if !confirmed {
            println!("{}", style("Process aborted by the user.").red());
            return;
        }
    }

    for path in a_paths.iter().chain(b_paths.iter()) {
        if !path.exists() {
            eprintln!("Path {} does not exist.", style(path.display()).red());
        }
        if !path.is_file() {
            eprintln!("Path {} is not a file.", style(path.display()).red());
        }
    }

    if output.is_dir() {
        eprintln!(
            "Output path {} is a directory.",
            style(output.display()).red()
        );
    }

    let context = match filters.context {
        DmrContext::CG => Context::CG,
        DmrContext::CHG => Context::CHG,
        DmrContext::CHH => Context::CHH,
    };
    let left_labels = repeat_n("A".to_string(), a_paths.len());
    let right_labels = repeat_n("B".to_string(), b_paths.len());

    let sample_paths = a_paths
        .into_iter()
        .chain(b_paths.into_iter())
        .map(|p| p.to_string_lossy().to_string())
        .collect::<Vec<_>>();
    let sample_labels = left_labels
        .into_iter()
        .chain(right_labels.into_iter())
        .collect::<Vec<_>>();

    let run_config = DmrConfig {
        context,
        n_missing: filters.n_missing,
        min_coverage: filters.min_coverage,
        diff_threshold: filters.diff_threshold,
        min_cpgs: filters.min_cytosines,
        segment_model: SegmentModelConfig {
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

    let mut progress_bar = if progress {
        let progress_bar = ProgressBar::new(dmr_iterator.blocks_total() as u64);
        progress_bar.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}, ETA: {eta}] [{bar:40.cyan/blue}] {pos:>5.green}/{len:5} {msg}").expect("Failed to create progress bar style")
                .progress_chars("#>-"),
        );
        progress_bar.set_message("Processing...");
        Some(progress_bar)
    } else {
        None
    };

    let mut dmrs = Vec::new();
    let mut p_values = Vec::new();
    let mut last_batch_idx = 0;
    let mut count = 0;

    let mut writer = BufWriter::new(
        File::create(
            output
                .with_added_extension("segments")
                .with_added_extension("tsv"),
        )
        .expect("Failed to open output file"),
    );
    for (batch_idx, dmr) in dmr_iterator {
        if let Some(progress_bar) = progress_bar.as_mut() {
            if batch_idx != last_batch_idx {
                progress_bar.inc(batch_idx as u64 - last_batch_idx as u64);
                last_batch_idx = batch_idx;
            }

            progress_bar.set_message(format!(
                "{}{}",
                style(format!("{}:", dmr.chr.as_str())).blue(),
                style(format!("{}-{}", dmr.start, dmr.end)).green(),
            ));
        }

        count += 1;

        if dmr.n_cytosines >= filters.min_cytosines && dmr.meth_diff() >= filters.diff_threshold {
            p_values.push(dmr.p_value);

            writeln!(
                &mut writer,
                "{}\t{}\t{}\t{}\t{}",
                dmr.chr.as_str(),
                dmr.start,
                dmr.end,
                dmr.p_value,
                dmr.meth_diff()
            )
            .unwrap();

            dmrs.push(dmr);
        }
    }

    // Calculate adjusted p-values using Benjamini-Hochberg procedure
    let m = p_values.len();
    let mut sorted_p_values: Vec<(usize, f64)> = p_values.iter().cloned().enumerate().collect();
    sorted_p_values.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut adjusted_p_values = vec![0.0; m];
    for (i, (index, p_value)) in sorted_p_values.iter().enumerate() {
        adjusted_p_values[*index] = (*p_value * m as f64 / (i + 1) as f64).min(1.0);
    }

    // Ensure the adjusted p-values are monotonic
    for i in (0..m - 1).rev() {
        if adjusted_p_values[i] > adjusted_p_values[i + 1] {
            adjusted_p_values[i] = adjusted_p_values[i + 1];
        }
    }

    // Assign adjusted p-values back to DMRs
    for (dmr, adj_p_value) in dmrs.iter_mut().zip(adjusted_p_values.iter()) {
        dmr.p_value = *adj_p_value;
    }

    let mut writer = BufWriter::new(
        File::create(
            output
                .with_added_extension("segments")
                .with_added_extension("tsv"),
        )
        .expect("Failed to open output file"),
    );
    for dmr in dmrs {
        if dmr.p_value < filters.q_value {
            writeln!(
                &mut writer,
                "{}\t{}\t{}\t{}\t{}",
                dmr.chr.as_str(),
                dmr.start,
                dmr.end,
                dmr.p_value,
                dmr.meth_diff()
            )
            .unwrap();
        }
    }

    writer.flush().unwrap()
}
