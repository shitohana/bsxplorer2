use std::fs::File;
use std::iter::repeat_n;
use std::path::PathBuf;

use anyhow::anyhow;
use bsxplorer2::data_structs::enums::Context;
use bsxplorer2::data_structs::typedef::{CountType, DensityType, PosType};
use bsxplorer2::tools::dmr::DMRegion;
use clap::{Args, ValueEnum};
use console::style;
use dialoguer::Confirm;
use indicatif::ProgressBar;
use serde::Serialize;

use crate::utils::{expand_wildcards, init_pbar, UtilsArgs};

#[derive(Args, Debug, Clone)]
pub(crate) struct DmrArgs {
    #[arg(
        value_parser,
        short = 'A',
        long,
        num_args=1..,
        required = true,
        help = "Paths to BSX files of the first sample group."
    )]
    group_a: Vec<String>,
    #[arg(
        value_parser,
        short = 'B',
        long,
        num_args=1..,
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
        default_value_t = Context::CG,
        help_heading = "FILTER ARGS",
        help = "Select cytosine methylation context. Only cytosines in this context will be used for DMR calling. CG/CHG/CHH."
    )]
    pub context: Context,

    #[arg(
        short,
        long,
        default_value_t = 0,
        help_heading = "FILTER ARGS",
        help = "Set missing values threshold. Cytosines with no data_structs in more \
                than specified number of samples will be discarded."
    )]
    pub n_missing: usize,

    #[arg(
        short = 'v',
        long,
        default_value_t = 5,
        help_heading = "FILTER ARGS",
        help = "Set coverage threshold. Cytosines with coverage below this value in \
                any of the samples will be discarded."
    )]
    pub min_coverage: CountType,

    #[arg(
        short,
        long,
        default_value_t = 10,
        help_heading = "FILTER ARGS",
        help = "Set minimum number of cytosines threshold. DMRs with cytosine count \
                below this value will be discarded."
    )]
    pub min_cytosines: usize,

    #[arg(
        short,
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "Set minimum difference threshold. DMRs with an absolute difference in \
                methylation proportion between the two groups smaller than this value \
                will be discarded."
    )]
    pub diff_threshold: DensityType,

    #[arg(
        short = 'D',
        long,
        default_value_t = 100,
        help_heading = "SEGMENTATION ARGS",
        help = "Maximum distance between adjacent cytosines in a segment.  Cytosines \
                further apart than this distance will be in separate segments."
    )]
    max_dist: PosType,

    #[arg(
        short = 'L',
        long,
        default_value_t = 2.0,
        help_heading = "SEGMENTATION ARGS",
        help = "Initial regularization parameter for the Condat algorithm.  Larger \
                values result in stronger smoothing."
    )]
    initial_l: f64,

    #[arg(
        short = 'l',
        long,
        default_value_t = 1e-3,
        help_heading = "SEGMENTATION ARGS",
        help = "Minimum value for the regularization parameter.  The regularization \
                parameter is decreased during segmentation until it is smaller than \
                this value."
    )]
    l_min: f64,

    #[arg(
        long = "coef",
        default_value_t = 1.5,
        help_heading = "SEGMENTATION ARGS",
        help = "Coefficient by which `initial_l` is divided in each iteration of the \
                segmentation algorithm. Smaller values perform more segmentation \
                iterations."
    )]
    l_coef: f64,

    #[arg(
        long,
        default_value_t = 1e-6,
        help_heading = "SEGMENTATION ARGS",
        help = "Tolerance for merging adjacent segments after the Total Variation \
                denoising step (Condat's algorithm).  Smaller values result in more \
                segments being merged. Should be very small to avoid \
                over-segmentation after denoising."
    )]
    tolerance: DensityType,

    #[arg(
        long,
        default_value_t = 1e-2,
        help_heading = "SEGMENTATION ARGS",
        help = "Mann-Whitney U-test P-value threshold for merging adjacent segments \
                during recursive segmentation. Smaller p-values result in more \
                iterations and fewer falsely merged segments."
    )]
    merge_p: f64,

    #[arg(
        short = 'p',
        long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = "Adjusted P-value threshold for DMR identification using \
                2D-Kolmogorov-Smirnov test. Segments with a p-value smaller than \
                specified will be reported as DMRs."
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

impl DmrArgs {
    pub fn run(
        &self,
        utils: &UtilsArgs,
    ) -> anyhow::Result<()> {
        let a_paths = expand_wildcards(self.group_a.clone());
        let b_paths = expand_wildcards(self.group_b.clone());

        if !self.force {
            let prompt = format!(
                "Do you want to proceed with the following paths?\n\nGroup A: \
                 {:?}\nGroup B: {:?}\nOutput:{:?}",
                a_paths, b_paths, self.output
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

        if a_paths.len() == 0 {
            return Err(anyhow::anyhow!("Group A files must not be empty"));
        }
        if b_paths.len() == 0 {
            return Err(anyhow::anyhow!("Group B files must not be empty"));
        }

        for path in a_paths.iter().chain(b_paths.iter()) {
            if !path.exists() {
                eprintln!("Path {} does not exist.", style(path.display()).red());
            }
            if !path.is_file() {
                eprintln!("Path {} is not a file.", style(path.display()).red());
            }
        }

        if self.output.is_dir() {
            eprintln!(
                "Output path {} is a directory.",
                style(self.output.display()).red()
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

        let run_config = bsxplorer2::tools::dmr::config::DmrConfig {
            context:        self.context,
            n_missing:      self.n_missing,
            min_coverage:   self.min_coverage,
            diff_threshold: self.diff_threshold,
            min_cpgs:       self.min_cytosines,
            max_dist:       self.max_dist,
            initial_l:      self.initial_l,
            l_min:          self.l_min,
            l_coef:         self.l_coef,
            seg_tolerance:  self.tolerance,
            merge_pvalue:   self.merge_p,
            seg_pvalue:     self.padj,
        };

        let files = sample_paths.iter().map(|path| {
            File::open(path).expect(format!("Could not open file {}", path).as_str())
        });
        let labeled = sample_labels.into_iter().zip(files).collect::<Vec<_>>();

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

        let all_segments_path = format!(
            "{}.segments.tsv",
            self.output
                .to_str()
                .unwrap_or_else(|| panic!("Path is empty!"))
        );

        let mut csv_writer = csv::WriterBuilder::default()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(&all_segments_path)
            .expect("Failed to open output file");

        let mut _dmr_count = 0;
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

            if dmr.n_cytosines >= self.min_cytosines {
                _dmr_count += 1;
                csv_writer
                    .serialize(dmr)
                    .expect("Failed to write DMR to file");
            }
        }

        csv_writer.flush().expect("Failed to write DMRs to file");

        progress_bar.finish();

        let all_segments = csv::ReaderBuilder::default()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(all_segments_path)
            .unwrap()
            .deserialize::<DMRegion>()
            .collect::<Result<Vec<_>, _>>()
            .expect("Failed to deserialize DMR segments file");

        let padj = if !matches!(self.pmethod, PadjMethod::None) {
            adjustp::adjust(
                &all_segments.iter().map(|s| s.p_value).collect::<Vec<_>>(),
                match self.pmethod {
                    PadjMethod::BH => adjustp::Procedure::BenjaminiHochberg,
                    PadjMethod::Bonf => adjustp::Procedure::Bonferroni,
                    PadjMethod::BY => adjustp::Procedure::BenjaminiYekutieli,
                    _ => unreachable!(),
                },
            )
        }
        else {
            all_segments.iter().map(|s| s.p_value).collect::<Vec<_>>()
        };

        let filtered = all_segments
            .into_iter()
            .zip(padj)
            .map(|(dmr, p)| DmrFilteredRow::from_dmr(dmr, p))
            .filter(|dmr| dmr.padj <= self.padj)
            .filter(|dmr| dmr.meth_diff.abs() >= self.diff_threshold)
            .filter(|dmr| dmr.n_cytosines >= self.min_cytosines)
            .collect::<Vec<_>>();
        let dmr_count = filtered.len();

        let filtered_path = format!(
            "{}.filtered.tsv",
            self.output
                .to_str()
                .unwrap_or_else(|| panic!("Path is empty!"))
        );
        let mut csv_writer = csv::WriterBuilder::default()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(&filtered_path)
            .expect("Failed to open output file");

        filtered
            .iter()
            .for_each(|dmr| csv_writer.serialize(dmr).unwrap());

        csv_writer.flush().expect("Failed to write DMRs to file");

        println!(
            "{}",
            style(format!("Found {} DMRs.", dmr_count)).green().bold()
        );

        Ok(())
    }
}

#[derive(Debug, Serialize)]
struct DmrFilteredRow {
    chr:         String,
    start:       u32,
    end:         u32,
    n_cytosines: usize,
    padj:        f64,
    p_utest:     f64,
    group_a:     f32,
    group_b:     f32,
    meth_diff:   f32,
    meth_mean:   f32,
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
