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
use std::fs::File;
use std::iter::repeat_n;
use std::path::PathBuf;

use bsxplorer2::exports::anyhow::anyhow;
use bsxplorer2::tools::dmr::DMRegion;
use clap::{Args, ValueEnum};
use console::style;
use dialoguer::Confirm;
use indicatif::ProgressBar;
use serde::Serialize;

use crate::utils::init_pbar;
use crate::{expand_wildcards,
            init_logger,
            init_rayon_threads,
            DmrContext,
            UtilsArgs};

#[derive(Args, Debug, Clone)]
pub(crate) struct DmrArgs {
    #[arg(
        value_parser,
        short = 'A',
        long,
        required = true,
        help = "Paths to BSX files of the first sample group."
    )]
    group_a: Vec<String>,
    #[arg(
        value_parser,
        short = 'B',
        long,
        required = true,
        help = "Paths to BSX files of the second sample group."
    )]
    group_b: Vec<String>,
    #[arg(
        short = 'o',
        long,
        required = true,
        help = "Prefix for the generated output files."
    )]
    output:  PathBuf,
    #[arg(
        short,
        long,
        required = false,
        default_value_t = false,
        help = "Automatically confirm selected paths."
    )]
    force:   bool,

    #[clap(
        short,
        long,
        value_enum,
        default_value_t = DmrContext::CG,
        help_heading = "FILTER ARGS",
        help = "Select cytosine methylation context. Only cytosines in this context will be used for DMR calling. CG/CHG/CHH."
    )]
    pub context: DmrContext,

    #[arg(
        short,
        long,
        default_value_t = 0,
        help_heading = "FILTER ARGS",
        help = "Set missing values threshold. Cytosines with no data_structs \
                in more than specified number of samples will be discarded."
    )]
    pub n_missing: usize,

    #[arg(
        short = 'v',
        long,
        default_value_t = 5,
        help_heading = "FILTER ARGS",
        help = "Set coverage threshold. Cytosines with coverage below this \
                value in any of the samples will be discarded."
    )]
    pub min_coverage: i16,

    #[arg(
        short,
        long,
        default_value_t = 10,
        help_heading = "FILTER ARGS",
        help = "Set minimum number of cytosines threshold. DMRs with cytosine \
                count below this value will be discarded."
    )]
    pub min_cytosines: usize,

    #[arg(
        short,
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "Set minimum difference threshold. DMRs with an absolute \
                difference in methylation proportion between the two groups \
                smaller than this value will be discarded."
    )]
    pub diff_threshold: f64,

    #[arg(
        short = 'D',
        long,
        default_value_t = 100,
        help_heading = "SEGMENTATION ARGS",
        help = "Maximum distance between adjacent cytosines in a segment.  \
                Cytosines further apart than this distance will be in \
                separate segments."
    )]
    max_dist: u32,

    #[arg(
        short = 'L',
        long,
        default_value_t = 2.0,
        help_heading = "SEGMENTATION ARGS",
        help = "Initial regularization parameter for the Condat algorithm.  \
                Larger values result in stronger smoothing."
    )]
    initial_l: f64,

    #[arg(
        short = 'l',
        long,
        default_value_t = 1e-3,
        help_heading = "SEGMENTATION ARGS",
        help = "Minimum value for the regularization parameter.  The \
                regularization parameter is decreased during segmentation \
                until it is smaller than this value."
    )]
    l_min: f64,

    #[arg(
        long = "coef",
        default_value_t = 1.5,
        help_heading = "SEGMENTATION ARGS",
        help = "Coefficient by which `initial_l` is divided in each iteration \
                of the segmentation algorithm. Smaller values perform more \
                segmentation iterations."
    )]
    l_coef: f64,

    #[arg(
        long,
        default_value_t = 1e-6,
        help_heading = "SEGMENTATION ARGS",
        help = "Tolerance for merging adjacent segments after the Total \
                Variation denoising step (Condat's algorithm).  Smaller \
                values result in more segments being merged. Should be very \
                small to avoid over-segmentation after denoising."
    )]
    tolerance: f64,

    #[arg(
        long,
        default_value_t = 1e-2,
        help_heading = "SEGMENTATION ARGS",
        help = "Mann-Whitney U-test P-value threshold for merging adjacent \
                segments during recursive segmentation. Smaller p-values \
                result in more iterations and fewer falsely merged segments."
    )]
    merge_p: f64,

    #[arg(
        short = 'p',
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "Adjusted P-value threshold for DMR identification using \
                2D-Kolmogorov-Smirnov test. Segments with a p-value smaller \
                than specified will be reported as DMRs."
    )]
    padj: f64,

    #[clap(
        long="pmethod",
        value_enum,
        default_value_t = PadjMethod::BH,
        help_heading = "FILTER ARGS",)]
    pmethod: PadjMethod,
}

#[derive(Debug, Clone, ValueEnum)]
pub(crate) enum PadjMethod {
    Bonf,
    BH,
    BY,
    None,
}

pub fn init(
    group_a: Vec<String>,
    group_b: Vec<String>,
    output: PathBuf,
    force: bool,
    threads: usize,
    logger: bool,
) -> bsxplorer2::exports::anyhow::Result<(Vec<String>, Vec<String>)> {
    init_rayon_threads(threads)?;
    init_logger(logger)?;

    let a_paths = expand_wildcards(group_a);
    let b_paths = expand_wildcards(group_b);

    if !force {
        let prompt = format!(
            "Do you want to proceed with the following paths?\n\nGroup A: \
             {:?}\nGroup B: {:?}\nOutput:{:?}",
            a_paths, b_paths, output
        );
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

pub fn run(
    args: DmrArgs,
    utils: UtilsArgs,
) {
    let (sample_paths, sample_labels) = if let Ok(result) = init(
        args.group_a,
        args.group_b,
        args.output.clone(),
        args.force,
        utils.threads,
        utils.verbose,
    ) {
        result
    }
    else {
        return;
    };
    let context = args.context.tobsxplorer2();

    let run_config = bsxplorer2::tools::dmr::config::DmrConfig {
        context,
        n_missing: args.n_missing,
        min_coverage: args.min_coverage,
        diff_threshold: args.diff_threshold,
        min_cpgs: args.min_cytosines,
        max_dist: args.max_dist as u64,
        initial_l: args.initial_l,
        l_min: args.l_min,
        l_coef: args.l_coef,
        seg_tolerance: args.tolerance,
        merge_pvalue: args.merge_p,
        seg_pvalue: args.padj,
    };

    let files = sample_paths.iter().map(|path| {
        File::open(path)
            .expect(format!("Could not open file {}", path).as_str())
    });
    let labeled = sample_labels
        .into_iter()
        .zip(files)
        .collect::<Vec<_>>();

    let dmr_iterator = run_config
        .try_finish(labeled)
        .expect("Failed to initialize DMR iterator");

    let progress_bar = if utils.progress {
        let progress_bar = init_pbar(dmr_iterator.blocks_total())
            .expect("Failed to initialize progress bar");
        progress_bar
    }
    else {
        ProgressBar::hidden()
    };

    let mut last_batch_idx = 0;

    let all_segments_path = args
        .output
        .with_added_extension("segments")
        .with_added_extension("tsv");

    let mut csv_writer = csv::WriterBuilder::default()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(&all_segments_path)
        .expect("Failed to open output file");

    let mut dmr_count = 0;
    for (batch_idx, dmr) in dmr_iterator {
        if batch_idx != last_batch_idx {
            progress_bar.inc(batch_idx as u64 - last_batch_idx as u64);
            last_batch_idx = batch_idx;
        }

        progress_bar.set_message(format!(
            "{}{}",
            style(format!("{}:", dmr.chr.as_str())).blue(),
            style(format!("{}-{}", dmr.start, dmr.end)).green(),
        ));

        if dmr.n_cytosines >= args.min_cytosines {
            dmr_count += 1;
            csv_writer
                .serialize(dmr)
                .expect("Failed to write DMR to file");
        }
    }

    csv_writer
        .flush()
        .expect("Failed to write DMRs to file");

    progress_bar.finish();

    let all_segments = csv::ReaderBuilder::default()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(all_segments_path)
        .unwrap()
        .deserialize::<DMRegion>()
        .collect::<Result<Vec<_>, _>>()
        .expect("Failed to deserialize DMR segments file");

    use bsxplorer2::exports::adjustp;
    let padj = if !matches!(args.pmethod, PadjMethod::None) {
        adjustp::adjust(
            &all_segments
                .iter()
                .map(|s| s.p_value)
                .collect::<Vec<_>>(),
            match args.pmethod {
                PadjMethod::BH => adjustp::Procedure::BenjaminiHochberg,
                PadjMethod::Bonf => adjustp::Procedure::Bonferroni,
                PadjMethod::BY => adjustp::Procedure::BenjaminiYekutieli,
                _ => unreachable!(),
            },
        )
    }
    else {
        all_segments
            .iter()
            .map(|s| s.p_value)
            .collect::<Vec<_>>()
    };

    let filtered = all_segments
        .into_iter()
        .zip(padj)
        .map(|(dmr, p)| DmrFilteredRow::from_dmr(dmr, p))
        .filter(|dmr| dmr.padj <= args.padj)
        .filter(|dmr| dmr.meth_diff.abs() >= args.diff_threshold)
        .filter(|dmr| dmr.n_cytosines >= args.min_cytosines)
        .collect::<Vec<_>>();
    let dmr_count = filtered.len();

    let filtered_path = args
        .output
        .with_added_extension("filtered")
        .with_added_extension("tsv");

    let mut csv_writer = csv::WriterBuilder::default()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(&filtered_path)
        .expect("Failed to open output file");

    filtered
        .iter()
        .for_each(|dmr| csv_writer.serialize(dmr).unwrap());

    csv_writer
        .flush()
        .expect("Failed to write DMRs to file");

    println!(
        "{}",
        style(format!("Found {} DMRs.", dmr_count))
            .green()
            .bold()
    );
}

#[derive(Debug, Serialize)]
struct DmrFilteredRow {
    chr:         String,
    start:       u32,
    end:         u32,
    n_cytosines: usize,
    padj:        f64,
    p_utest:     f64,
    group_a:     f64,
    group_b:     f64,
    meth_diff:   f64,
    meth_mean:   f64,
}

impl DmrFilteredRow {
    fn from_dmr(
        dmr: DMRegion,
        padj: f64,
    ) -> Self {
        Self {
            padj,
            chr: dmr.chr,
            start: dmr.start,
            end: dmr.end,
            n_cytosines: dmr.n_cytosines,
            p_utest: dmr.p_value,
            group_a: dmr.meth_left,
            group_b: dmr.meth_right,
            meth_diff: dmr.meth_diff,
            meth_mean: dmr.meth_mean,
        }
    }
}
