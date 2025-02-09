use crate::{expand_wildcards, DmrContext};
use _lib::exports::anyhow::anyhow;
use _lib::exports::itertools::Itertools;
use _lib::exports::serde::Serialize;
use _lib::tools::dmr::dmr_fast::DmrConfig;
use _lib::tools::dmr::penalty_segment::PenaltySegmentModel;
use _lib::tools::dmr::DMRegion;
use _lib::utils::types::Context;
use clap::Args;
use console::style;
use dialoguer::Confirm;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::iter::repeat_n;
use std::path::PathBuf;

#[derive(Args, Debug, Clone)]
pub(crate) struct FilterArgs {
    #[clap(short, long, value_enum, default_value_t = DmrContext::CG, help_heading="FILTER ARGS", help = "Select cytosine methylation context.")]
    pub(crate) context: DmrContext,
    #[arg(
        short,
        long,
        default_value_t = 0,
        help_heading = "FILTER ARGS",
        help = "Set missing values threshold. Cytosines with no data in > n_samples will be discarded."
    )]
    pub(crate) n_missing: usize,
    #[arg(
        short = 'v',
        long,
        default_value_t = 5,
        help_heading = "FILTER ARGS",
        help = "Cytosines with coverage below threshold will be discarded."
    )]
    pub(crate) min_coverage: i16,
    #[arg(
        short,
        long,
        default_value_t = 10,
        help_heading = "FILTER ARGS",
        help = "DMRs with cytosine count below threshold will be discarded."
    )]
    pub(crate) min_cytosines: usize,
    #[arg(
        short,
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "DMRs with difference between samples below threshold will be discarded."
    )]
    pub(crate) diff_threshold: f64,
    #[arg(
        short = 'q',
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "Q-value for DMR identification."
    )]
    pub(crate) q_value: f64,
}

#[derive(Args, Debug, Clone)]
pub(crate) struct SegmentationArgs {
    #[arg(
        long,
        default_value_t = 100,
        help_heading = "SEGMENTATION ARGS",
        help = "Number of optimal lambda approximation iterations."
    )]
    nlambda: usize,
    #[arg(
        long,
        default_value_t = 1e-3,
        help_heading = "SEGMENTATION ARGS",
        help = "Initial lambda value. Affects segments length. Must be at least 10 times less than expected methylation density range."
    )]
    initial_l: f64,
    #[arg(
        long,
        default_value_t = 1e-6,
        help_heading = "SEGMENTATION ARGS",
        help = "Tolerance of segment identifying after denoising step (should be low)."
    )]
    tolerance: f64,
    #[arg(
        long,
        default_value_t = 1e5,
        help_heading = "SEGMENTATION ARGS",
        help = "Segment length penalty in optimal lambda approximation. Longer segments -> less penalty."
    )]
    length_penalty: f64,
    #[arg(
        long,
        default_value_t = 1e4,
        help_heading = "SEGMENTATION ARGS",
        help = "Segments count penalty in optimal lambda approximation. Fewer segments -> less penalty."
    )]
    count_penalty: f64,
    #[arg(
        long,
        default_value_t = 1e-2,
        help_heading = "SEGMENTATION ARGS",
        help = "Welch t-Test P-value for segments merging step. Smaller p-value -> more iterations, less falsely merged segments."
    )]
    merge_p: f64,
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
    let (sample_paths, sample_labels) =
        if let Ok(result) = init(group_a, group_b, output.clone(), force, threads) {
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
        segment_model: PenaltySegmentModel {
            min_seg_length: filters.min_cytosines,
            seg_count_weight: segmentation.count_penalty,
            penalty_weight: segmentation.length_penalty,
            segmentation_tol: segmentation.tolerance,
            merge_pvalue: segmentation.merge_p,
            num_lambdas: segmentation.nlambda,
            lambda_low: segmentation.initial_l,
            union_threshold: segmentation.union_threshold,
            max_dist: segmentation.max_dist,
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

pub(crate) fn init_pbar(total: usize) -> _lib::exports::anyhow::Result<ProgressBar> {
    let progress_bar = ProgressBar::new(total as u64);
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}, ETA: {eta}] [{bar:40.cyan/blue}] {pos:>5.green}/{len:5} {msg}")?
            .progress_chars("#>-"),
    );
    progress_bar.set_message("Processing...");
    Ok(progress_bar)
}

pub(crate) fn init(
    group_a: Vec<String>,
    group_b: Vec<String>,
    output: PathBuf,
    force: bool,
    threads: usize,
) -> _lib::exports::anyhow::Result<(Vec<String>, Vec<String>)> {
    _lib::exports::rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;
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
            return Err(anyhow!("User aborted the process."));
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
    Ok((sample_paths, sample_labels))
}

pub(crate) fn run_dmr(
    iterator: Box<dyn Iterator<Item = (usize, DMRegion)>>,
    mut progress_bar: Option<ProgressBar>,
    filters: FilterArgs,
    output: PathBuf,
) {
    let mut p_values = Vec::new();
    let mut last_batch_idx = 0;

    let all_segments_path = output
        .with_added_extension("segments")
        .with_added_extension("tsv");

    let mut csv_writer = csv::WriterBuilder::default()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(&all_segments_path)
        .expect("Failed to open output file");
    for (batch_idx, dmr) in iterator {
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

        if dmr.n_cytosines >= filters.min_cytosines {
            p_values.push(dmr.p_value);
            csv_writer
                .serialize(dmr)
                .expect("Failed to write DMR to file");
        }
    }

    csv_writer.flush().expect("Failed to write DMRs to file");

    let dmrs = csv::ReaderBuilder::default()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(&all_segments_path)
        .expect("Failed to open output file")
        .deserialize::<DMRegion>()
        .collect::<Result<Vec<_>, _>>()
        .expect("Failed to read DMRs from file");

    let adjusted_p_values = _lib::exports::adjustp::adjust(
        &*dmrs.iter().map(|dmr| dmr.p_value).collect_vec(),
        _lib::exports::adjustp::Procedure::BenjaminiHochberg,
    );

    let filtered_path = output
        .with_added_extension("filtered")
        .with_added_extension("tsv");

    let mut csv_writer = csv::WriterBuilder::default()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(&filtered_path)
        .expect("Failed to open output file");

    let mut dmr_count = 0;
    for (mut dmr, padj) in dmrs.into_iter().zip(adjusted_p_values.into_iter()) {
        dmr.p_value = padj;
        if dmr.meth_diff() >= filters.diff_threshold && dmr.p_value <= filters.q_value {
            dmr_count += 1;
            p_values.push(dmr.p_value);
            csv_writer
                .serialize(dmr)
                .expect("Failed to write DMR to file");
        }
    }

    progress_bar.map(|pb| pb.finish());
    println!(
        "{}",
        style(format!("Found {} DMRs.", dmr_count)).green().bold()
    );
}
