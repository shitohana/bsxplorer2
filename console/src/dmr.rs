use crate::{expand_wildcards, DmrContext};
use _lib::tools::dmr::{MethyLassoConfig, MethylLassoRunConfig, RegionType, SegmentModel};
use _lib::utils::types::Context;
use clap::Args;
use console::style;
use dialoguer::Confirm;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::iter::repeat_n;
use std::path::PathBuf;

#[derive(Args, Debug, Clone)]
pub(crate) struct FilterArgs {
    #[clap(short, long, value_enum, default_value_t = DmrContext::CG, help_heading="FILTER ARGS", help = "Select cytosine methylation context.")]
    context: DmrContext,
    #[arg(
        short,
        long,
        default_value_t = 0,
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
        short,
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "P-value for Welch t-Test for DMR identification. Regions with p-value above threshold will be discarded."
    )]
    p_value: f64,
    #[arg(
        short,
        long,
        default_value_t = 10,
        help_heading = "FILTER ARGS",
        help = "DMRs with cytosine count below threshold will be discarded."
    )]
    min_cytosines: usize,
    #[arg(
        short,
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "DMRs with difference between samples below threshold will be discarded."
    )]
    diff_threshold: f64,
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
}

pub fn run(
    filters: FilterArgs,
    segmentation: SegmentationArgs,
    group_a: Vec<String>,
    group_b: Vec<String>,
    output: PathBuf,
    force: bool,
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

    let run_config = MethylLassoRunConfig {
        analysis_config: MethyLassoConfig {
            context,
            n_missing: filters.n_missing,
            min_coverage: filters.min_coverage,
            diff_threshold: filters.diff_threshold,
            p_value: filters.p_value,
            type_density: Default::default(),
            min_cpgs: filters.min_cytosines,
            segment_model: SegmentModel {
                min_seg_length: filters.min_cytosines,
                num_lambdas: segmentation.nlambda,
                lambda_low: segmentation.initial_l,
                penalty_weight: segmentation.length_penalty,
                seg_count_weight: segmentation.count_penalty,
                merge_pvalue: segmentation.merge_p,
                segmentation_tol: segmentation.tolerance,
            },
        },

        sample_paths,
        sample_labels,
        selected_regions: vec![RegionType::DMR],
        output: output.to_string_lossy().to_string(),
    };

    let files = run_config
        .sample_paths
        .iter()
        .map(|path| File::open(path).expect(format!("Could not open file {}", path).as_str()))
        .map(BufReader::new)
        .collect::<Vec<_>>();

    let iterator = run_config
        .analysis_config
        .clone()
        .finish(
            run_config
                .sample_labels
                .iter()
                .cloned()
                .zip(files.into_iter())
                .collect::<Vec<_>>(),
        )
        .expect("Failed to initialize sample reader");

    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
            .template("{spinner} {msg}")
            .expect("Failed to set spinner template"),
    );
    spinner.set_message("Processing...");

    let mut sink =
        BufWriter::new(File::create(&run_config.output).expect("Failed to open output file"));

    let mut dmr_count = 0;
    for region in iterator {
        if run_config.selected_regions.contains(&region.region_type) {
            spinner.set_message(format!("Found {} DMRs", style(dmr_count).green().bold()));
            dmr_count += 1;

            if dmr_count % 10 == 0 {
                spinner.inc(1);
            }

            let meth_diff = region.mean_right - region.mean_left;
            let overall = (region.mean_left + region.mean_right) / 2.0;
            let n_cytosines = region.pairwise_sites.len();
            let ref_id = region.pairwise_sites.left.chr.as_str();
            let start = region.pairwise_sites.left.positions.first().unwrap();
            let end = region.pairwise_sites.left.positions.last().unwrap();
            let region_type = region.region_type;
            let mean_a = region.mean_left;
            let mean_b = region.mean_right;
            let pvalue = region.p_value;

            let line = format!("{ref_id}\t{start}\t{end}\t{n_cytosines}\t{region_type:?}\t{pvalue:.3e}\t{mean_a:.5}\t{mean_b:.5}\t{meth_diff:.5}\t{overall}");
            writeln!(sink, "{line}").expect("Failed to write line");
        }
    }
}
