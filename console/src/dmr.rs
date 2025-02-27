use crate::utils::init_pbar;
use crate::{expand_wildcards, DmrContext};
use _lib::exports::anyhow::anyhow;
use clap::Args;
use console::style;
use dialoguer::Confirm;
use std::fs::File;
use std::iter::repeat_n;
use std::path::PathBuf;

#[derive(Args, Debug, Clone)]
pub(crate) struct MetileneArgs {
    #[clap(short, long, value_enum, default_value_t = DmrContext::CG, help_heading="FILTER ARGS", help = "Select cytosine methylation context.")]
    pub context: DmrContext,

    #[arg(
        short,
        long,
        default_value_t = 0,
        help_heading = "FILTER ARGS",
        help = "Set missing values threshold. Cytosines with no data in > n_samples will be discarded."
    )]
    pub n_missing: usize,

    #[arg(
        short = 'v',
        long,
        default_value_t = 5,
        help_heading = "FILTER ARGS",
        help = "Cytosines with coverage below threshold will be discarded."
    )]
    pub min_coverage: i16,

    #[arg(
        short,
        long,
        default_value_t = 10,
        help_heading = "FILTER ARGS",
        help = "DMRs with cytosine count below threshold will be discarded."
    )]
    pub min_cytosines: usize,

    #[arg(
        short,
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "DMRs with difference between samples below threshold will be discarded."
    )]
    pub diff_threshold: f64,

    #[arg(
        short = 'D',
        long,
        default_value_t = 100,
        help_heading = "SEGMENTATION ARGS",
        help = "Maximum distance between adjacent cytosines in a segment."
    )]
    max_dist: u32,

    #[arg(
        short = 'L',
        long,
        default_value_t = 2.0,
        help_heading = "SEGMENTATION ARGS"
    )]
    initial_l: f64,

    #[arg(
        short = 'l',
        long,
        default_value_t = 1e-3,
        help_heading = "SEGMENTATION ARGS"
    )]
    l_min: f64,

    #[arg(
        long = "coef",
        default_value_t = 1.5,
        help_heading = "SEGMENTATION ARGS"
    )]
    l_coef: f64,

    #[arg(
        long,
        default_value_t = 1e-6,
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
        short = 'p',
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "P-value for DMR identification."
    )]
    seg_pvalue: f64,
}

pub fn init(
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

pub fn run(
    args: MetileneArgs,
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
    let context = args.context.to_lib();

    let run_config = _lib::tools::dmr::config::DmrConfig {
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
        seg_pvalue: args.seg_pvalue,
    };

    let files = sample_paths
        .iter()
        .map(|path| File::open(path).expect(format!("Could not open file {}", path).as_str()));
    let labeled = sample_labels.into_iter().zip(files).collect::<Vec<_>>();

    let dmr_iterator = run_config
        .try_finish(labeled)
        .expect("Failed to initialize DMR iterator");

    let mut progress_bar = if progress {
        let progress_bar =
            init_pbar(dmr_iterator.blocks_total()).expect("Failed to initialize progress bar");

        Some(progress_bar)
    } else {
        None
    };

    let mut last_batch_idx = 0;

    let all_segments_path = output
        .with_added_extension("segments")
        .with_added_extension("tsv");

    let mut csv_writer = csv::WriterBuilder::default()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(&all_segments_path)
        .expect("Failed to open output file");

    let mut dmr_count = 0;
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

            if dmr.n_cytosines >= args.min_cytosines {
                dmr_count += 1;
                csv_writer
                    .serialize(dmr)
                    .expect("Failed to write DMR to file");
            }
        }
    }

    csv_writer.flush().expect("Failed to write DMRs to file");

    progress_bar.map(|pb| pb.finish());
    println!(
        "{}",
        style(format!("Found {} DMRs.", dmr_count)).green().bold()
    );
}
