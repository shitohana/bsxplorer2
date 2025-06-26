use std::fs::File;
use std::path::{
    Path,
    PathBuf,
};

use anyhow::{
    anyhow,
    Context as AnyhowContext,
    Result,
};
use bsxplorer2::data_structs::typedef::{
    CountType,
    DensityType,
    PosType,
};
use bsxplorer2::data_structs::Context;
use bsxplorer2::prelude::MultiBsxFileReader;
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
use csv::WriterBuilder;
use dialoguer::Confirm;
use log::{
    debug,
    info,
};
use serde::Serialize;
use spipe::spipe;

use crate::strings::dmr as strings;
use crate::utils::{
    check_validate_paths,
    init_progress,
    init_readers,
};
use crate::{
    exit_with_msg,
    PipelineCommand,
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

#[derive(Debug, Clone, ValueEnum, Copy)]
pub(crate) enum PadjMethod {
    Bonf,
    BH,
    BY,
    None,
}

fn get_csv_writer(path: &Path) -> Result<csv::Writer<File>> {
    WriterBuilder::default()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .map_err(|e| anyhow!(e))
        .context(format!("{}", path.to_string_lossy()))
}

fn get_csv_reader(path: &Path) -> Result<csv::Reader<File>> {
    csv::ReaderBuilder::default()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .map_err(|e| anyhow!(e))
        .context(format!("{}", path.to_string_lossy()))
}

fn write_dmr_segments(
    left_paths: &[String],
    right_paths: &[String],
    config: DmrConfig,
    output: PathBuf,
) -> Result<PathBuf> {
    let left_reader = spipe!(
        left_paths
            =>  check_validate_paths()
            =>& init_readers
            =>? MultiBsxFileReader::from_iter
    );
    info!("Initialized readers for {}", left_paths.join(","));

    let right_reader = spipe!(
        right_paths
            =>  check_validate_paths()
            =>& init_readers
            =>? MultiBsxFileReader::from_iter
    );
    info!("Initialized readers for {}", left_paths.join(","));

    let mut dmr_iterator = DmrReader::new(config.clone(), (left_reader, right_reader));
    debug!("Initialized DMR iterator");

    let segments_path = output.join(".segments.tsv");
    info!(
        "Segments will be written to {}",
        segments_path.to_string_lossy()
    );

    let mut csv_writer = get_csv_writer(segments_path.as_path())?;
    info!("Initialized CSV writer for segments");

    let progress_bar = init_progress(Some(dmr_iterator.blocks_total()))?;

    let mut iter = dmr_iterator.iter();
    debug!("Starting DMR iteration");

    while let Some(dmr) = iter.next() {
        let batch_idx = iter.batch_num();
        let dmr = dmr?;
        debug!("Received DMR from batch {}: {}", batch_idx, dmr.as_contig());

        if batch_idx as u64 != progress_bar.position() {
            progress_bar.set_position(batch_idx as u64);
        }

        if dmr.n_cytosines >= config.min_coverage as usize {
            csv_writer
                .serialize(dmr)
                .context(format!("Batch idx: {batch_idx}"))?;
            debug!("DMR written to file");
        }
        else {
            debug!("Short DMR - skipping")
        }
    }
    csv_writer.flush()?;
    info!("CSV writer finished");

    progress_bar.finish_with_message(format!(
        "Segments have been calculated and written to: {}",
        output.to_string_lossy()
    ));

    info!("Done writing DMR segments");
    Ok(segments_path)
}

fn write_filtered_segments(
    filtered_rows: Vec<DmrFilteredRow>,
    output: PathBuf,
) -> Result<PathBuf> {
    let mut writer = get_csv_writer(output.as_path())?;
    info!("Initialized writing filtered DMRs to {}", output.display());

    for row in filtered_rows {
        writer.serialize(row)?;
    }

    writer.flush()?;
    info!("Finished writing filtered DMRs to {}", output.display());

    Ok(output.to_owned())
}

fn read_dmr_segments(path: &Path) -> Result<Vec<DMRegion>> {
    get_csv_reader(path)?
        .deserialize::<DMRegion>()
        .collect::<Result<Vec<_>, csv::Error>>()
        .map_err(|e| anyhow!(e))
}

fn filter_dmrs(
    dmrs: Vec<DMRegion>,
    pmethod: PadjMethod,
    padj_threshold: f64,
    diff_threshold: DensityType,
    min_cytosines: usize,
) -> Vec<DmrFilteredRow> {
    let padj = if !matches!(pmethod, PadjMethod::None) {
        adjustp::adjust(
            &dmrs.iter().map(|s| s.p_value).collect::<Vec<_>>(),
            match pmethod {
                PadjMethod::BH => adjustp::Procedure::BenjaminiHochberg,
                PadjMethod::Bonf => adjustp::Procedure::Bonferroni,
                PadjMethod::BY => adjustp::Procedure::BenjaminiYekutieli,
                _ => unreachable!(),
            },
        )
    }
    else {
        dmrs.iter().map(|s| s.p_value).collect::<Vec<_>>()
    };

    if !matches!(pmethod, PadjMethod::None) {
        debug!("Calculated adjusted p-values");
    }

    let filtered = dmrs
        .into_iter()
        .zip(padj)
        .map(|(dmr, p)| DmrFilteredRow::from_dmr(dmr, p))
        .filter(|dmr| dmr.padj <= padj_threshold)
        .filter(|dmr| dmr.meth_diff.abs() >= diff_threshold)
        .filter(|dmr| dmr.n_cytosines >= min_cytosines)
        .collect::<Vec<_>>();

    debug!("Done filtering DMRs");

    filtered
}

fn confirm_paths(
    left_paths: &[String],
    right_paths: &[String],
    output: &Path,
    force: bool,
) {
    if !force {
        let prompt = format!(
            "Do you want to proceed with the following paths?\n\nGroup A: {:?}\nGroup \
             B: {:?}\nOutput:{:?}",
            left_paths, right_paths, output
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

impl PipelineCommand for DmrArgs {
    fn run(&self) -> anyhow::Result<()> {
        confirm_paths(&self.group_a, &self.group_b, &self.output, self.force);

        let segments_path = write_dmr_segments(
            &self.group_a,
            &self.group_b,
            DmrConfig::from(self),
            self.output.clone(),
        )?;

        let dmrs = read_dmr_segments(&segments_path)?;

        let filtered = filter_dmrs(
            dmrs,
            self.pmethod,
            self.padj,
            self.diff_threshold,
            self.min_cytosines,
        );
        let dmr_count = filtered.len();

        let output_path = write_filtered_segments(filtered, self.output.clone())?;

        println!(
            "Found {} DMRs.\nWritten segments to: {}\nWritten filtered DMRs to: {}",
            style(dmr_count).green().bold(),
            segments_path.display(),
            output_path.display()
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
