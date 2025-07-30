use std::fs::File;
use std::path::PathBuf;
use std::sync::{
    Arc,
    Mutex,
};

use bsxplorer2::data_structs::batch::{
    AggMethod,
    BsxBatch,
};
use bsxplorer2::data_structs::coords::ContigIntervalMap;
use bsxplorer2::data_structs::typedef::PosType;
use bsxplorer2::prelude::*;
use bsxplorer2::tools::dimred::merge::EqFloat;
use bsxplorer2::tools::dimred::{
    pelt,
    MethDataBinom,
};
use bsxplorer2::utils::{
    BoundThreadExecutor,
    THREAD_POOL,
};
use clap::Args;
use crossbeam::channel::Sender;
use hashbrown::HashMap;
use itertools::{
    izip,
    Itertools,
};
use log::{
    debug,
    error,
    info,
};
use rayon::prelude::*;
use uuid::Uuid;

use crate::dimred::write_imap;
use crate::strings::dimred as strings;
use crate::utils::{
    check_validate_paths,
    init_progress,
    init_readers,
};
use crate::PipelineCommand;

#[derive(Args, Debug, Clone)]
pub(crate) struct DimRedArgs {
    #[arg(value_parser, num_args=1.., required = true, help = strings::FILES)]
    files:        Vec<String>,
    #[arg(short, long, required = true, help = strings::OUTPUT)]
    output:       PathBuf,
    #[arg(short, long, help = strings::CONTEXT, default_value_t = Context::CG)]
    context:      Context,
    #[arg(short = 'C', long, help = strings::COVERAGE, default_value_t = 5)]
    coverage:     u16,
    #[arg(short, long, help = strings::BETA)]
    beta:         Option<f64>,
    #[arg(short, long, default_value_t = 20, help = strings::MIN_SIZE)]
    min_size:     usize,
    #[arg(long, default_value_t = 10000, help = strings::CHUNK)]
    chunk:        usize,
    #[arg(long, default_value_t = 1000, help = strings::INTERSECTION)]
    intersection: usize,
    #[arg(long, help = strings::JOINT, default_value_t = 5)]
    joint:        usize,
}

fn read_thread(
    mut reader: MultiBsxFileReader,
    sender: Sender<BsxBatch>,
    min_coverage: u16,
    context: Context,
) {
    let pbar = init_progress(Some(reader.blocks_total())).expect("Should not fail");

    THREAD_POOL.spawn(move || {
        pbar.wrap_iter(reader.iter_merged(AggMethod::Sum, AggMethod::Mean))
            .enumerate()
            .filter_map(|(idx, batch_res)| {
                match batch_res {
                    Result::Err(e) => {
                        eprintln!("Error processing batch {}: {}", idx, e);
                        None
                    },
                    Result::Ok(batch) => Some(batch),
                }
            })
            .map(|b| {
                b.lazy()
                    .filter_coverage_gt(min_coverage)
                    .filter_context(context)
                    .collect()
                    .expect("Should not fail")
            })
            .filter(BsxBatch::is_empty)
            .try_for_each(|b| sender.send(b))
            .expect("Failed to send batch to main thread! This is unexpected.");

        pbar.finish();
    });
}

struct DimRedRunner {
    cache: Option<(String, BsxBatch)>,
    args:  DimRedArgs,
}

impl DimRedRunner {
    fn new(args: DimRedArgs) -> anyhow::Result<Self> {
        Ok(DimRedRunner { cache: None, args })
    }

    fn run(mut self) -> anyhow::Result<()> {
        let reader = spipe::spipe!(
            &self.args.files
                =>  check_validate_paths()
                =>& init_readers
                =>? MultiBsxFileReader::from_iter
        );
        info!("Opened MultiBsxFileReader");

        let (sc, rc) = crossbeam::channel::bounded(THREAD_POOL.current_num_threads());
        let executor = BoundThreadExecutor::new(&THREAD_POOL);

        read_thread(reader, sc, self.args.coverage, self.args.context);
        info!("Started reader thread");

        let results_map = Arc::new(Mutex::new(HashMap::new()));
        let mut cyt_sum = 0;

        // Process batches in parallel
        while let Ok(batch) = rc.recv() {
            cyt_sum += batch.len();

            // Process batch if cache conditions are met
            if let Some((current_chunk, is_final_chunk)) = self.update_cache(batch) {
                let results_map = Arc::clone(&results_map);
                executor.install(move || {
                    process_and_store(
                        current_chunk,
                        results_map.lock().as_deref_mut().unwrap(),
                        self.args.beta,
                        self.args.min_size,
                        is_final_chunk,
                    );
                });
            }
        }

        executor.join();

        // Process any remaining cached data
        if let Some((_, remaining_batch)) = self.cache.take() {
            process_and_store(
                remaining_batch,
                results_map.lock().as_deref_mut().unwrap(),
                self.args.beta,
                self.args.min_size,
                true,
            );
        }
        info!("Finished processing batches");

        let contig_imap = merge_changepoints(
            Arc::try_unwrap(results_map)
                .expect("Should not fail")
                .into_inner()?,
            self.args.joint,
        );
        print_stats(cyt_sum, contig_imap.n_intervals());
        info!("Changepoints merged");

        write_imap(File::create(self.args.output.as_path())?, &contig_imap)?;
        info!("Result written to {}", self.args.output.display());

        Ok(())
    }

    fn update_cache(
        &mut self,
        batch: BsxBatch,
    ) -> Option<(BsxBatch, bool)> {
        let this_chr = unsafe { batch.seqname().unwrap_unchecked().to_string() };

        match self.cache.as_mut() {
            Some((cached_chr, cached_batch)) if cached_chr == &this_chr => {
                // Same chromosome - combine batches
                unsafe { cached_batch.extend_unchecked(&batch) };

                // Check if combined batch is long enough to process
                if cached_batch.len() >= self.args.chunk + self.args.intersection {
                    // Split the batch: process the first part, keep the remainder
                    let split_point = cached_batch.len() - self.args.intersection;
                    let to_process = cached_batch.slice(0, cached_batch.len());
                    let to_cache =
                        cached_batch.slice(split_point as i64, cached_batch.len());

                    // Update cache with remainder
                    *cached_batch = to_cache;
                    Some((to_process, false))
                }
                else {
                    None
                }
            },
            Some((_cached_chr, cached_batch)) => {
                // Different chromosome - process cached batch and start new cache
                let output = std::mem::replace(cached_batch, batch);
                Some((output, true))
            },
            None => {
                // No cache - start caching
                self.cache = Some((this_chr, batch));
                None
            },
        }
    }
}

fn print_stats(
    before_dimred: usize,
    after_dimred: usize,
) {
    let compression_ratio = if before_dimred > 0 {
        (after_dimred as f64 / before_dimred as f64) * 100.0
    }
    else {
        0.0
    };

    println!();
    println!("  Total cytosines processed: {}", before_dimred);
    println!("  Total changepoints identified: {}", after_dimred);
    println!(
        "  Compression ratio: {:.2}% of original data",
        compression_ratio
    );
}

fn detect_changepoints(
    batch: &BsxBatch,
    beta: Option<f64>,
    min_size: usize,
) -> Option<(Vec<PosType>, Vec<f64>)> {
    if batch.count_m().null_count() != 0 || batch.count_total().null_count() != 0 {
        return None;
    }

    let positions = batch.positions_vec();
    let meth_data = MethDataBinom::from(batch);
    let beta = beta.unwrap_or((batch.len() as f64).ln());

    let (seg_ends, _score) = pelt(&meth_data, beta, min_size);

    // We add +1 because pelt points to last pos
    // of a segment, and we need to get start of segment
    let seg_pos = seg_ends.iter().map(|&idx| positions[idx + 1]).collect_vec();
    let densities = batch
        .partition_density(&seg_ends, AggMethod::Mean.get_fn())
        .ok()?;
    Some((seg_pos, densities))
}

fn process_and_store(
    batch: BsxBatch,
    storage: &mut HashMap<String, Vec<(PosType, f64, Uuid)>>,
    beta: Option<f64>,
    min_size: usize,
    is_final: bool,
) {
    if let Some((mut positions, mut densities)) =
        detect_changepoints(&batch, beta, min_size)
    {
        if is_final {
            positions.push(batch.last_pos().unwrap());
        }
        else {
            densities.pop();
        }

        // We can unwrap as we know batch is not empty
        let chromosome = batch.seqname().unwrap().to_string();
        debug!("Processed final? {is_final} batch {batch:?}");

        storage.entry(chromosome).or_insert(Vec::new()).extend(
            izip!(positions, densities)
                .map(|(pos, density)| (pos, density, uuid::Uuid::new_v4())),
        )
    }
    else {
        error!("Could not process batch {}", batch)
    }
}


impl PipelineCommand for DimRedArgs {
    fn run(&self) -> anyhow::Result<()> {
        DimRedRunner::new(self.clone())?.run()
    }
}

fn merge_changepoints<I: Eq + Send + Copy>(
    map: HashMap<String, Vec<(PosType, f64, I)>>,
    r: usize,
) -> ContigIntervalMap<EqFloat> {
    let mut result = HashMap::new();

    for (chr, mut changepoints_vec) in map {
        let mut merged = Vec::new();
        changepoints_vec.par_sort_by_key(|(p, ..)| *p);

        let mut i = 0;
        while i < changepoints_vec.len() {
            let (_start_pos, start_density, _start_id) = changepoints_vec[i];
            let mut end_idx = i;
            let mut densities = vec![start_density];

            // Look for consecutive elements with different chunk_ids that are within
            // distance r
            while end_idx + 1 < changepoints_vec.len() {
                let (curr_pos, _curr_density, curr_id) = changepoints_vec[end_idx];
                let (next_pos, next_density, next_id) = changepoints_vec[end_idx + 1];

                if curr_id != next_id && (next_pos - curr_pos) <= r as PosType {
                    densities.push(next_density);
                    end_idx += 1;
                }
                else {
                    break;
                }
            }

            // Calculate mean density
            let mean_density = densities.iter().sum::<f64>() / densities.len() as f64;
            let mean_density = match EqFloat::try_from(mean_density) {
                Ok(v) => v,
                Err(e) => {
                    error!("Invalid density value {}: {}", mean_density, e);
                    Default::default()
                },
            };

            // Use the position of the last element in the group
            let final_pos = changepoints_vec[end_idx].0;
            merged.push((final_pos, mean_density));

            i = end_idx + 1;
        }

        result.insert(chr.clone(), merged);
    }

    ContigIntervalMap::from_breakpoints(result)
}
