use std::collections::{
    BTreeMap,
    HashMap,
};
use std::fs::File;
use std::io::IsTerminal;
use std::path::PathBuf;
use std::sync::{
    Arc,
    RwLock,
};

use anyhow::{
    anyhow,
    Ok,
};
use bsxplorer2::data_structs::batch::{
    self,
    AggMethod,
    BsxBatch,
};
use bsxplorer2::data_structs::coords::Contig;
use bsxplorer2::data_structs::typedef::{
    BsxSmallStr,
    PosType,
};
use bsxplorer2::data_structs::{
    Context,
    Strand,
};
use bsxplorer2::io::bsx::{
    BsxFileReader,
    MultiBsxFileReader,
};
use bsxplorer2::tools::dimred::{
    pelt,
    MethDataBinom,
    SegmentationData,
};
use clap::Args;
use console::style;
use indicatif::ProgressBar;
use itertools::{
    izip,
    Itertools,
};
use std_semaphore::Semaphore;

use crate::utils::{
    expand_wildcards_single,
    init_pbar,
};

#[derive(Args, Debug, Clone)]
pub(crate) struct DimRedArgs {
    #[arg(
        value_parser,
        num_args=1..,
        required = true,
        help = "Paths to BSX files"
    )]
    files:        Vec<String>,
    #[arg(short, long, required = true, help = "Path to output file")]
    output:       PathBuf,
    #[arg(
        short, long,
        help = "Methylation context",
        default_value_t = Context::CG
    )]
    context:      Context,
    #[arg(short = 'C', long, help = "Minimal coverage", default_value_t = 5)]
    coverage:     u16,
    #[arg(
        short,
        long,
        help = "Beta parameter value for PELT algorithm. If none - BIC (ln(length)) \
                is selected."
    )]
    beta:         Option<f64>,
    #[arg(
        short,
        long,
        default_value_t = 20,
        help = "Minimal segment size (number of cytosines)"
    )]
    min_size:     usize,
    #[arg(long, default_value_t = 10000, help = "Chunk size for PELT algorithm")]
    chunk:        usize,
    #[arg(
        long,
        help = "Batch intersection size for parallel processing",
        default_value_t = 1000
    )]
    intersection: usize,
    #[arg(
        long,
        help = "Changepoint joint distance for parallel processing",
        default_value_t = 5
    )]
    joint:        usize,
}

macro_rules! read_lock {
    ($lock: expr) => {
        $lock
            .read()
            .map_err(|e| anyhow!(e.to_string()).context("Failed to acquire read lock"))
    };
}

macro_rules! write_lock {
    ($lock: expr) => {
        $lock
            .write()
            .map_err(|e| anyhow!(e.to_string()).context("Failed to acquire write lock"))
    };
}

struct DimRedRunner {
    cache:    Option<(String, BsxBatch)>,
    progress: ProgressBar,
    args:     DimRedArgs,
}

impl DimRedRunner {
    fn new(args: DimRedArgs) -> anyhow::Result<Self> {
        let progress = if std::io::stdin().is_terminal() {
            init_pbar(0)?
        }
        else {
            ProgressBar::hidden()
        };

        Ok(DimRedRunner {
            cache: None,
            args,
            progress,
        })
    }

    fn run(self) -> anyhow::Result<()> {
        // Expand file paths and validate input files
        let paths = self
            .args
            .files
            .iter()
            .flat_map(|path| expand_wildcards_single(path))
            .collect::<Vec<_>>();

        // Initialize file readers
        let mut reader = MultiBsxFileReader::from_iter(
            paths
                .clone()
                .into_iter()
                .map(|path| {
                    File::open(&path).map_err(|e| {
                        anyhow!("Failed to open file '{}': {}", path.display(), e)
                    })
                })
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .map(|file| {
                    BsxFileReader::try_new(file)
                        .map_err(|e| anyhow!("Failed to create BSX reader: {}", e))
                })
                .collect::<Result<Vec<_>, _>>()?,
        );

        self.progress.set_length(reader.blocks_total() as u64);
        self.progress
            .set_message(format!("Processing {} BSX files...", paths.len()));

        // Initialize shared state for parallel processing
        let args = Arc::new(self.args.clone());
        let self_ref = Arc::new(RwLock::new(self));
        let results = Arc::new(crossbeam::queue::SegQueue::new());
        let semaphore = Arc::new(Semaphore::new(
            bsxplorer2::utils::THREAD_POOL.current_num_threads() as isize,
        ));
        let total_cytosines_processed = Arc::new(RwLock::new(0));

        // Process batches in parallel
        bsxplorer2::utils::THREAD_POOL.scope(|s| -> anyhow::Result<()> {
            for (batch_idx, batch_res) in reader
                .iter_merged(batch::AggMethod::Sum, AggMethod::Mean)
                .enumerate()
            {
                let batch = batch_res?;
                let filtered_batch = read_lock!(self_ref)?.filter(batch)?;

                if filtered_batch.is_empty() {
                    continue;
                }

                read_lock!(self_ref)?.progress.inc(1);
                *write_lock!(total_cytosines_processed)? += filtered_batch.len();

                // Process batch if cache conditions are met
                if let Some((batch_to_process, is_final_segment)) =
                    write_lock!(self_ref)?.update_cache(filtered_batch)?
                {
                    let args_clone = Arc::clone(&args);
                    let results_clone = Arc::clone(&results);
                    let semaphore_clone = Arc::clone(&semaphore);

                    semaphore_clone.acquire();
                    s.spawn(move |_| {
                        let chromosome = unsafe {
                            batch_to_process.seqname().unwrap_unchecked().to_string()
                        };

                        match process_batch(
                            &batch_to_process,
                            &args_clone,
                            is_final_segment,
                        ) {
                            Result::Ok((positions, densities)) => {
                                results_clone.push((
                                    batch_idx, chromosome, positions, densities,
                                ));
                            },
                            Err(e) => {
                                eprintln!(
                                    "Error processing batch {}: {}",
                                    batch_idx, e
                                );
                            },
                        }

                        semaphore_clone.release();
                    });
                }
            }

            Ok(())
        })?;

        // Process any remaining cached data
        if let Some((_, remaining_batch)) = read_lock!(self_ref)?.cache.as_ref() {
            let chromosome =
                unsafe { remaining_batch.seqname().unwrap_unchecked().to_string() };

            let (positions, densities) =
                process_batch(&remaining_batch, &args, true)
                    .map_err(|e| anyhow!("Error processing final batch: {}", e))?;

            results.push((usize::MAX, chromosome, positions, densities));
        }

        // Collect and organize results
        let mut changepoints = ChangepointMap::new();
        while let Some((batch_idx, chromosome, positions, densities)) = results.pop() {
            changepoints.entry(chromosome.into()).or_default().extend(
                izip!(positions, densities)
                    .map(|(pos, density)| (pos, (batch_idx, density))),
            );
        }

        // Merge nearby changepoints across parallel segments
        let merged_changepoints =
            merge_nearby_changepoints(changepoints, read_lock!(self_ref)?.args.joint);

        let total_changepoints = merged_changepoints
            .values()
            .map(|changepoints| changepoints.len())
            .sum::<usize>();

        // Write results to BED format
        read_lock!(self_ref)?.progress.finish_and_clear();

        let output_path = read_lock!(self_ref)?.args.output.clone();
        let mut bed_writer =
            bio::io::bed::Writer::new(File::create(&output_path).map_err(|e| {
                anyhow!(
                    "Failed to create output file '{}': {}",
                    output_path.display(),
                    e
                )
            })?);

        for (chromosome, chromosome_changepoints) in merged_changepoints {
            let mut previous_position = 0;

            for (position, methylation_density) in chromosome_changepoints {
                let segment = Contig::new(
                    chromosome.clone(),
                    previous_position,
                    position,
                    Strand::None,
                );
                let mut bed_record: bio::io::bed::Record = segment.into();
                bed_record.set_score(&methylation_density.to_string());
                bed_writer.write(&bed_record)?;

                previous_position = position;
            }
        }

        // Print summary statistics
        let total_cytosines = *read_lock!(total_cytosines_processed)?;
        let compression_ratio = if total_cytosines > 0 {
            (total_changepoints as f64 / total_cytosines as f64) * 100.0
        }
        else {
            0.0
        };

        eprintln!();
        eprintln!("Dimensionality reduction completed successfully!");
        eprintln!("  Total cytosines processed: {}", total_cytosines);
        eprintln!("  Total changepoints identified: {}", total_changepoints);
        eprintln!(
            "  Compression ratio: {:.2}% of original data",
            compression_ratio
        );
        eprintln!("  Output file: {}", style(output_path.display()).green());

        Ok(())
    }

    fn filter(
        &self,
        batch: BsxBatch,
    ) -> anyhow::Result<BsxBatch> {
        let filtered = batch
            .lazy()
            .filter_context(self.args.context)
            .filter_coverage_gt(self.args.coverage.saturating_sub(1))
            .collect()?;

        if filtered.count_m().null_count() != 0
            || filtered.count_total().null_count() != 0
        {
            anyhow::bail!("Can not segment region. Missing counts data");
        }

        Ok(filtered)
    }

    fn update_cache(
        &mut self,
        batch: BsxBatch,
    ) -> anyhow::Result<Option<(BsxBatch, bool)>> {
        let current_chr = unsafe { batch.seqname().unwrap_unchecked().to_string() };

        match self.cache.as_mut() {
            Some((cached_chr, cached_batch)) if cached_chr == &current_chr => {
                // Same chromosome - combine batches
                unsafe { cached_batch.extend_unchecked(&batch)? };

                // Check if combined batch is long enough to process
                if cached_batch.len() >= self.args.chunk + self.args.intersection {
                    // Split the batch: process the first part, keep the remainder
                    let split_point = cached_batch.len() - self.args.intersection;
                    let to_process = cached_batch.slice(0, cached_batch.len());
                    let to_cache =
                        cached_batch.slice(split_point as i64, cached_batch.len());

                    // Update cache with remainder
                    *cached_batch = to_cache;
                    return Ok(Some((to_process, false)));
                }
                else {
                    return Ok(None);
                }
            },
            Some((_cached_chr, cached_batch)) => {
                // Different chromosome - process cached batch and start new cache
                let output = std::mem::replace(cached_batch, batch);
                Ok(Some((output, true)))
            },
            None => {
                // No cache - start caching
                self.cache = Some((current_chr, batch));
                Ok(None)
            },
        }
    }
}

fn process_batch(
    batch: &BsxBatch,
    args: &DimRedArgs,
    is_last: bool,
) -> anyhow::Result<(Vec<PosType>, Vec<f64>)> {
    let (positions, count_m, count_total) = unsafe {
        (
            batch
                .position()
                .to_vec_null_aware()
                .left()
                .unwrap_unchecked(),
            batch
                .count_m()
                .to_vec_null_aware()
                .left()
                .unwrap_unchecked(),
            batch
                .count_total()
                .to_vec_null_aware()
                .left()
                .unwrap_unchecked(),
        )
    };
    let meth_data = MethDataBinom::new(&count_m, &count_total);

    let (segment_boundaries, _score) = pelt(
        &meth_data,
        args.beta.unwrap_or((meth_data.len() as f64).ln()),
        args.min_size,
    );
    let mut selected_positions = segment_boundaries
        .iter()
        .map(|&idx| positions[idx + 1])
        .collect_vec();
    let mut densities =
        batch.partition_density(&segment_boundaries, AggMethod::Mean.get_fn())?;

    if is_last {
        selected_positions.push(batch.last_pos().unwrap());
    }
    else {
        densities.pop();
    }

    debug_assert_eq!(
        selected_positions.len(),
        densities.len(),
        "Length mismatch {selected_positions:?}, {densities:?}, {is_last}"
    );
    Ok((selected_positions, densities))
}

impl DimRedArgs {
    pub fn run(&self) -> anyhow::Result<()> {
        DimRedRunner::new(self.clone())?.run()
    }
}

/// Global genomic position -> (chunk_id, density)
type ChrChangepoints = BTreeMap<PosType, (usize, f64)>;
/// Chromosome name -> Chromosome changepoints
type ChangepointMap = HashMap<BsxSmallStr, ChrChangepoints>;

fn merge_nearby_changepoints(
    map: ChangepointMap,
    r: usize,
) -> HashMap<BsxSmallStr, Vec<(PosType, f64)>> {
    let mut result = HashMap::new();

    for (chr, changepoints) in map {
        let mut merged = Vec::new();
        let changepoints_vec: Vec<(PosType, (usize, f64))> = changepoints
            .iter()
            .map(|(&pos, &(chunk_id, density))| (pos, (chunk_id, density)))
            .collect();

        if changepoints_vec.is_empty() {
            result.insert(chr.clone(), merged);
            continue;
        }

        let mut i = 0;
        while i < changepoints_vec.len() {
            let (_start_pos, (_start_chunk_id, start_density)) = changepoints_vec[i];
            let mut end_idx = i;
            let mut densities = vec![start_density];

            // Look for consecutive elements with different chunk_ids that are within
            // distance r
            while end_idx + 1 < changepoints_vec.len() {
                let (curr_pos, (curr_chunk_id, _curr_density)) =
                    changepoints_vec[end_idx];
                let (next_pos, (next_chunk_id, next_density)) =
                    changepoints_vec[end_idx + 1];

                if curr_chunk_id != next_chunk_id
                    && (next_pos - curr_pos) <= r as PosType
                {
                    densities.push(next_density);
                    end_idx += 1;
                }
                else {
                    break;
                }
            }

            // Calculate mean density
            let mean_density = densities.iter().sum::<f64>() / densities.len() as f64;

            // Use the position of the last element in the group
            let final_pos = changepoints_vec[end_idx].0;
            merged.push((final_pos, mean_density));

            i = end_idx + 1;
        }

        result.insert(chr.clone(), merged);
    }

    result
}

#[cfg(test)]
mod tests {
    use std::collections::{
        BTreeMap,
        HashMap,
    };

    use super::*;

    #[test]
    fn test_merge_nearby_changepoints_empty_map() {
        let map = HashMap::new();
        let result = merge_nearby_changepoints(map, 10);
        assert!(result.is_empty());
    }

    #[test]
    fn test_merge_nearby_changepoints_empty_chromosome() {
        let mut map = HashMap::new();
        map.insert("chr1".into(), BTreeMap::new());

        let result = merge_nearby_changepoints(map, 10);
        assert_eq!(result.len(), 1);
        assert!(result.get("chr1").unwrap().is_empty());
    }

    #[test]
    fn test_merge_nearby_changepoints_single_point() {
        let mut map = HashMap::new();
        let mut changepoints = BTreeMap::new();
        changepoints.insert(100, (0, 0.5));
        map.insert("chr1".into(), changepoints);

        let result = merge_nearby_changepoints(map, 10);
        assert_eq!(result.len(), 1);
        let chr1_result = result.get("chr1").unwrap();
        assert_eq!(chr1_result.len(), 1);
        assert_eq!(chr1_result[0], (100, 0.5));
    }

    #[test]
    fn test_merge_nearby_changepoints_same_chunk() {
        let mut map = HashMap::new();
        let mut changepoints = BTreeMap::new();
        changepoints.insert(100, (0, 0.5));
        changepoints.insert(105, (0, 0.6));
        changepoints.insert(110, (0, 0.7));
        map.insert("chr1".into(), changepoints);

        let result = merge_nearby_changepoints(map, 10);
        assert_eq!(result.len(), 1);
        let chr1_result = result.get("chr1").unwrap();
        // Should not merge points from same chunk
        assert_eq!(chr1_result.len(), 3);
    }

    #[test]
    fn test_merge_nearby_changepoints_different_chunks_within_distance() {
        let mut map = HashMap::new();
        let mut changepoints = BTreeMap::new();
        changepoints.insert(100, (0, 0.5));
        changepoints.insert(105, (1, 0.7));
        map.insert("chr1".into(), changepoints);

        let result = merge_nearby_changepoints(map, 10);
        assert_eq!(result.len(), 1);
        let chr1_result = result.get("chr1").unwrap();
        assert_eq!(chr1_result.len(), 1);
        assert_eq!(chr1_result[0], (105, 0.6)); // Mean of 0.5 and 0.7
    }

    #[test]
    fn test_merge_nearby_changepoints_different_chunks_beyond_distance() {
        let mut map = HashMap::new();
        let mut changepoints = BTreeMap::new();
        changepoints.insert(100, (0, 0.5));
        changepoints.insert(120, (1, 0.7));
        map.insert("chr1".into(), changepoints);

        let result = merge_nearby_changepoints(map, 10);
        assert_eq!(result.len(), 1);
        let chr1_result = result.get("chr1").unwrap();
        assert_eq!(chr1_result.len(), 2);
        assert_eq!(chr1_result[0], (100, 0.5));
        assert_eq!(chr1_result[1], (120, 0.7));
    }

    #[test]
    fn test_merge_nearby_changepoints_multiple_consecutive_merges() {
        let mut map = HashMap::new();
        let mut changepoints = BTreeMap::new();
        changepoints.insert(100, (0, 0.2));
        changepoints.insert(103, (1, 0.4));
        changepoints.insert(106, (2, 0.6));
        changepoints.insert(109, (3, 0.8));
        map.insert("chr1".into(), changepoints);

        let result = merge_nearby_changepoints(map, 5);
        assert_eq!(result.len(), 1);
        let chr1_result = result.get("chr1").unwrap();
        assert_eq!(chr1_result.len(), 1);
        assert_eq!(chr1_result[0], (109, 0.5)); // Mean of 0.2, 0.4, 0.6, 0.8
    }

    #[test]
    fn test_merge_nearby_changepoints_mixed_scenarios() {
        let mut map = HashMap::new();
        let mut changepoints = BTreeMap::new();
        changepoints.insert(100, (0, 0.1));
        changepoints.insert(103, (1, 0.3));
        changepoints.insert(106, (2, 0.5));
        changepoints.insert(150, (3, 0.7));
        changepoints.insert(153, (4, 0.9));
        map.insert("chr1".into(), changepoints);

        let result = merge_nearby_changepoints(map, 5);
        assert_eq!(result.len(), 1);
        let chr1_result = result.get("chr1").unwrap();
        assert_eq!(chr1_result.len(), 2);
        assert_eq!(chr1_result[0], (106, 0.3)); // Mean of 0.1, 0.3, 0.5
        assert_eq!(chr1_result[1], (153, 0.8)); // Mean of 0.7, 0.9
    }

    #[test]
    fn test_merge_nearby_changepoints_multiple_chromosomes() {
        let mut map = HashMap::new();

        let mut chr1_changepoints = BTreeMap::new();
        chr1_changepoints.insert(100, (0, 0.5));
        chr1_changepoints.insert(105, (1, 0.7));
        map.insert("chr1".into(), chr1_changepoints);

        let mut chr2_changepoints = BTreeMap::new();
        chr2_changepoints.insert(200, (0, 0.3));
        chr2_changepoints.insert(210, (1, 0.9));
        map.insert("chr2".into(), chr2_changepoints);

        let result = merge_nearby_changepoints(map, 10);
        assert_eq!(result.len(), 2);

        let chr1_result = result.get("chr1").unwrap();
        assert_eq!(chr1_result.len(), 1);
        assert_eq!(chr1_result[0], (105, 0.6));

        let chr2_result = result.get("chr2").unwrap();
        assert_eq!(chr2_result.len(), 1);
        assert_eq!(chr2_result[0], (210, 0.6));
    }
}
