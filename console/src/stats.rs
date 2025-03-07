use crate::utils::init_pbar;
use crate::{init_logger, init_rayon_threads, UtilsArgs};
use bsxplorer2::data_structs::region_data::RegionData;
use bsxplorer2::exports::itertools::Itertools;
use bsxplorer2::io::bsx::read::BsxFileReader;
use bsxplorer2::io::bsx::region_read::RegionReader;
use bsxplorer2::tools::meth_stats::{MethylationStatFlat, MethylationStats};
use bsxplorer2::utils::types::{IPCEncodedEnum, Strand};
use clap::{Args, ValueEnum};
use console::style;
use serde::{Serialize, Serializer};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

#[derive(Args, Debug, Clone)]
pub(crate) struct StatsArgs {
    #[arg(help = "Path of the input file.")]
    input: PathBuf,
    #[arg(
        short = 'o',
        long,
        required = true,
        help = "Path for the generated output file."
    )]
    output: PathBuf,
    #[clap(short, long, value_enum, default_value_t = StatsMode::Genomewide, help = "Stats mode.")]
    mode: StatsMode,
    #[clap(short, long, value_enum, default_value_t = AnnotationFormat::Gff, help = "Annotation format.")]
    format: AnnotationFormat,
    #[arg(
        short = 'a',
        long,
        required = false,
        help = "Path for the generated output file."
    )]
    annot_path: Option<PathBuf>,
    #[arg(
        long,
        required = false,
        help = "Feature type to filter.",
        default_value = "gene"
    )]
    feature_type: Option<String>,
    #[arg(
        long,
        required = false,
        help = "Number of threads to use.",
        default_value_t = 1
    )]
    threads: usize,
}

#[derive(Debug, Clone, ValueEnum)]
enum StatsMode {
    Genomewide,
    Regions,
}

#[derive(Debug, Clone, ValueEnum)]
enum AnnotationFormat {
    Gff,
    Bed,
}

pub(crate) fn run(args: StatsArgs, utils: UtilsArgs) {
    init(&args);
    init_rayon_threads(utils.threads).expect("Failed to create global thread pool");
    init_logger(utils.verbose).expect("Failed to set up logger");

    match args.mode {
        StatsMode::Genomewide => {
            let reader = BsxFileReader::new(
                File::open(&args.input).expect("Error: failed to open input file."),
            );

            let mut initial_stats = MethylationStats::new();
            let pbar =
                init_pbar(reader.blocks_total()).expect("Error: failed to create progress bar.");

            for batch in reader {
                let batch = batch.expect("Error: failed to read batch.");
                let batch_stats = batch
                    .get_methylation_stats()
                    .expect("Error: failed to get methylation stats.");
                initial_stats.merge(&batch_stats);
                pbar.inc(1);
            }
            initial_stats.finalize_methylation();
            let json = bsxplorer2::exports::serde_json::to_string_pretty(&initial_stats)
                .expect("Serialization failed");
            let mut file =
                File::create(&args.output).expect("Error: failed to create output file.");
            file.write_all(json.as_bytes())
                .expect("Error: failed to write to output file.");
            pbar.finish();
            pbar.set_message("Done.");
        }
        StatsMode::Regions => {
            use bsxplorer2::utils::types::Strand;

            let pbar = init_pbar(0).expect("Error: failed to create progress bar.");
            pbar.set_message("Reading annotation...");

            let annotation = match args.format {
                AnnotationFormat::Gff => {
                    use bsxplorer2::exports::bio::bio_types::strand::Strand as BioStrand;
                    use bsxplorer2::exports::bio::io::gff::*;

                    let mut reader = Reader::from_file(&args.annot_path.unwrap(), GffType::GFF3)
                        .expect("Error: failed to open annotation file.");

                    reader
                        .records()
                        .into_iter()
                        .flatten()
                        .filter(|r| {
                            if let Some(filter) = args.feature_type.as_ref() {
                                r.feature_type() == filter
                            } else {
                                true
                            }
                        })
                        .map(|r| {
                            let start = r.start();
                            let end = r.end();
                            let chr = r.seqname().to_string();
                            let strand = r
                                .strand()
                                .map(|s| match s {
                                    BioStrand::Forward => Strand::Forward,
                                    BioStrand::Reverse => Strand::Reverse,
                                    BioStrand::Unknown => Strand::None,
                                })
                                .unwrap_or(Strand::None);
                            let attributes = HashMap::from_iter(
                                r.attributes()
                                    .into_iter()
                                    .map(|(k, v)| (k.clone(), v.join(", "))),
                            );
                            RegionData::new(chr, *start, *end, strand, (), attributes)
                        })
                        .collect_vec()
                }
                AnnotationFormat::Bed => {
                    use bsxplorer2::exports::bio::bio_types::strand::Strand as BioStrand;
                    use bsxplorer2::exports::bio::io::bed::*;

                    let mut reader = Reader::from_file(&args.annot_path.unwrap())
                        .expect("Error: failed to open annotation file.");

                    reader
                        .records()
                        .into_iter()
                        .flatten()
                        .map(|r| {
                            let start = r.start();
                            let end = r.end();
                            let chr = r.chrom().to_string();
                            let strand = r
                                .strand()
                                .map(|s| match s {
                                    BioStrand::Forward => Strand::Forward,
                                    BioStrand::Reverse => Strand::Reverse,
                                    BioStrand::Unknown => Strand::None,
                                })
                                .unwrap_or(Strand::None);
                            RegionData::new(chr, start, end, strand, (), Default::default())
                        })
                        .collect_vec()
                }
            };

            pbar.set_message("Reading regions...");
            pbar.set_length(annotation.len() as u64);

            let region_reader = RegionReader::try_new(
                File::open(&args.input).expect("Error: failed to open input file."),
                &annotation,
            )
            .expect("Error: failed to create region reader.");

            let mut writer = csv::WriterBuilder::default()
                .delimiter(b'\t')
                .has_headers(true)
                .from_path(&args.output)
                .expect("Error: failed to create output file.");

            for region_data in region_reader {
                let mut stats = region_data
                    .data()
                    .get_methylation_stats()
                    .expect("Error: failed to get methylation stats.");
                stats.finalize_methylation();
                let stats_row: MethylationStatFlat = stats.into();

                let stats_row = RegionRow {
                    chr: region_data.chr().clone(),
                    start: region_data.start(),
                    end: region_data.end(),
                    strand: region_data.strand(),
                    mean_methylation: stats_row.mean_methylation,
                    methylation_var: stats_row.methylation_var,
                    mean_coverage: stats_row.mean_coverage,
                    cg_mean_methylation: stats_row.cg_mean_methylation,
                    cg_coverage: stats_row.cg_coverage,
                    chg_mean_methylation: stats_row.chg_mean_methylation,
                    chg_coverage: stats_row.chg_coverage,
                    chh_mean_methylation: stats_row.chh_mean_methylation,
                    chh_coverage: stats_row.chh_coverage,
                    fwd_mean_methylation: stats_row.fwd_mean_methylation,
                    fwd_coverage: stats_row.fwd_coverage,
                    rev_mean_methylation: stats_row.rev_mean_methylation,
                    rev_coverage: stats_row.rev_coverage,
                };
                writer
                    .serialize(stats_row)
                    .expect("Error: failed to write to output file.");

                pbar.inc(1);
            }
            pbar.finish();
            writer.flush().unwrap();
            pbar.set_message("Done.");
        }
    }
}

fn serialize_strand<S>(strand: &Strand, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(strand.to_string().as_str())
}

// Unfortunately, csv_serde does not support serde(flatten) attribute
#[derive(Serialize)]
struct RegionRow {
    chr: String,
    start: u64,
    end: u64,
    #[serde(serialize_with = "serialize_strand")]
    strand: Strand,
    mean_methylation: f64,
    methylation_var: f64,
    mean_coverage: f64,
    cg_mean_methylation: f64,
    cg_coverage: u32,
    chg_mean_methylation: f64,
    chg_coverage: u32,
    chh_mean_methylation: f64,
    chh_coverage: u32,
    fwd_mean_methylation: f64,
    fwd_coverage: u32,
    rev_mean_methylation: f64,
    rev_coverage: u32,
}

fn init(args: &StatsArgs) {
    bsxplorer2::exports::rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .expect("Error: failed to initialize thread pool.");
    if !args.input.exists() {
        eprintln!(
            "Error: input file {} not found.",
            style(args.input.to_str().unwrap()).red()
        );
        std::process::exit(1);
    }
    if !args.input.is_file() {
        eprintln!(
            "Error: input file {} is not a file.",
            style(args.input.to_str().unwrap()).red()
        );
        std::process::exit(1);
    }
    if args.output.is_dir() {
        eprintln!(
            "Error: output file {} is a directory.",
            style(args.output.to_str().unwrap()).red()
        );
        std::process::exit(1);
    }
    if matches!(args.mode, StatsMode::Regions) {
        if args.annot_path.is_none() {
            eprintln!("Error: missing required argument `annot_path` for `regions` mode.");
            std::process::exit(1);
        }
        if !args.annot_path.as_ref().unwrap().exists() {
            eprintln!(
                "Error: annotation file {} not found.",
                style(args.annot_path.as_ref().unwrap().to_str().unwrap()).red()
            );
            std::process::exit(1);
        }
    }
    if matches!(args.mode, StatsMode::Genomewide) {
        if args.annot_path.is_some() {
            println!(
                "{}",
                style("Warning: ignoring `annot_path` argument for `genomewide` mode.").red()
            );
        }
    }
}
