use std::fs::File;
use std::path::PathBuf;

use bsxplorer2::data_structs::typedef::{
    CountType,
    DensityType,
    PosType,
};
use bsxplorer2::data_structs::Context;
use bsxplorer2::tools::dmr::{
    DMRegion,
    DmrConfig,
    DmrReader,
};
use clap::{
    Args,
    ValueEnum,
};
use console::style;
use dialoguer::Confirm;
use indicatif::ProgressBar;
use serde::Serialize;

use crate::strings::dmr as strings;
use crate::utils::{
    expand_wildcards,
    init_pbar,
};
use crate::{
    assert_or_exit,
    exit_with_msg,
    init_progress,
};

#[derive(Args, Debug, Clone)]
pub(crate) struct DmrArgs {
    #[arg(
        value_parser,
        short = 'A', long,
        num_args=1..,
        required = true,
        help = strings::GROUP_A
    )]
    group_a: Vec<String>,
    #[arg(
        value_parser,
        short = 'B', long,
        num_args=1..,
        required = true,
        help = strings::GROUP_B
    )]
    group_b: Vec<String>,
    #[arg(
        short = 'o', long,
        required = true,
        help = strings::OUTPUT
    )]
    output:  PathBuf,
    #[arg(
        short, long,
        default_value_t = false,
        help = strings::FORCE
    )]
    force:   bool,

    #[clap(
        short, long,
        value_enum,
        default_value_t = Context::CG,
        help_heading = "FILTER ARGS",
        help = strings::CONTEXT
    )]
    context: Context,

    #[arg(
        short, long,
        default_value_t = 0,
        help_heading = "FILTER ARGS",
        help = strings::N_MISSING
    )]
    n_missing: usize,

    #[arg(
        short = 'v', long,
        default_value_t = 5,
        help_heading = "FILTER ARGS",
        help = strings::MIN_COVERAGE
    )]
    min_coverage: CountType,

    #[arg(
        short, long,
        default_value_t = 10,
        help_heading = "FILTER ARGS",
        help = strings::MIN_CYTOSINES
    )]
    min_cytosines: usize,

    #[arg(
        short, long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = strings::DIFF_THRESHOLD
    )]
    diff_threshold: DensityType,

    #[arg(
        short = 'D', long,
        default_value_t = 100,
        help_heading = "SEGMENTATION ARGS",
        help = strings::MAX_DIST
    )]
    max_dist: PosType,

    #[arg(
        short = 'L', long,
        default_value_t = 2.0,
        help_heading = "SEGMENTATION ARGS",
        help = strings::INITIAL_L
    )]
    initial_l: f64,

    #[arg(
        short = 'l', long,
        default_value_t = 1e-3,
        help_heading = "SEGMENTATION ARGS",
        help = strings::L_MIN
    )]
    l_min: f64,

    #[arg(
        long = "coef",
        default_value_t = 1.5,
        help_heading = "SEGMENTATION ARGS",
        help = strings::L_COEF
    )]
    l_coef: f64,

    #[arg(
        long,
        default_value_t = 1e-6,
        help_heading = "SEGMENTATION ARGS",
        help = strings::TOLERANCE
    )]
    tolerance: DensityType,

    #[arg(
        long,
        default_value_t = 1e-2,
        help_heading = "SEGMENTATION ARGS",
        help = strings::MERGE_P
    )]
    merge_p: f64,

    #[arg(
        short = 'p', long,
        default_value_t = 0.05,
        help_heading = "FILTER ARGS",
        help = strings::PADJ
    )]
    padj: f64,

    #[clap(
        long="pmethod",
        value_enum,
        default_value_t = PadjMethod::BH,
        help_heading = "FILTER ARGS"
    )]
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
    pub fn run(&self) -> anyhow::Result<()> {
        let a_paths = expand_wildcards(self.group_a.clone());
        let b_paths = expand_wildcards(self.group_b.clone());
        self.confirm(&a_paths, &b_paths);

        assert_or_exit!(
            !self.output.is_dir(),
            "Output {} is a directory.",
            self.output.display()
        );
        assert_or_exit!(
            self.output.to_str().is_some_and(|v| !v.is_empty()),
            "Output prefix is empty!"
        );

        assert_or_exit!(a_paths.len() != 0, "Group {} files must not be empty", "A");
        assert_or_exit!(b_paths.len() != 0, "Group {} files must not be empty", "B");
        for path in a_paths.iter().chain(b_paths.iter()) {
            assert_or_exit!(path.exists(), "Path {} does not exist.", path.display());
            assert_or_exit!(path.is_file(), "Path {} is not a file.", path.display());
        }

        let mut dmr_iterator = DmrReader::from_readers(
            a_paths
                .into_iter()
                .map(|p| {
                    File::open(&p)
                        .expect(format!("Could not open file {:?}", p).as_str())
                })
                .collect::<Vec<_>>(),
            b_paths
                .into_iter()
                .map(|p| {
                    File::open(&p)
                        .expect(format!("Could not open file {:?}", p).as_str())
                })
                .collect::<Vec<_>>(),
            DmrConfig::from(self),
        )
        .expect("Failed to initialize DMR iterator");

        let progress_bar = init_progress!(dmr_iterator.blocks_total())?;

        let mut last_batch_idx = 0;

        let all_segments_path = self.output.join(".segments.tsv");

        let mut csv_writer = csv::WriterBuilder::default()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(&all_segments_path)
            .expect("Failed to open output file");

        let mut _dmr_count = 0;
        let mut iter = dmr_iterator.iter();

        while let Some(dmr) = iter.next() {
            let dmr = dmr.expect("Error during DMR identification.");

            let batch_idx = iter.batch_num();
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

    fn confirm(
        &self,
        a_paths: &[PathBuf],
        b_paths: &[PathBuf],
    ) {
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
                exit_with_msg!("Process aborted by the user.");
            }
        }
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

impl From<&DmrArgs> for DmrConfig {
    fn from(value: &DmrArgs) -> Self {
        bsxplorer2::tools::dmr::DmrConfig {
            context:        value.context,
            n_missing:      value.n_missing,
            min_coverage:   value.min_coverage,
            diff_threshold: value.diff_threshold,
            min_cpgs:       value.min_cytosines,
            max_dist:       value.max_dist,
            initial_l:      value.initial_l,
            l_min:          value.l_min,
            l_coef:         value.l_coef,
            seg_tolerance:  value.tolerance,
            merge_pvalue:   value.merge_p,
            seg_pvalue:     value.padj,
        }
    }
}
