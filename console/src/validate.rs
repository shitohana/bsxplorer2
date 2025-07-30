use std::collections::HashMap;
use std::fmt::{
    Debug,
    Display,
};
use std::fs::File;
use std::hash::Hash;
use std::process::exit;
use std::rc::Rc;

use bsxplorer2::data_structs::batch::BsxBatch;
use bsxplorer2::io::bsx::BsxFileReader;
use clap::Args;
use console::style;
use indicatif::ProgressBar;
use itertools::Itertools;
use polars::error::PolarsResult;

use crate::utils::{
    expand_wildcard,
    init_progress,
};
use crate::PipelineCommand;

#[derive(Args, Debug, Clone)]
pub(crate) struct ValidateArgs {
    #[arg(
        value_parser,
        num_args=2..,
        required = true,
        help = "Paths to BSX files"
    )]
    files: Vec<String>,

    #[arg(short, long, default_value_t = false, help = "Validate each row")]
    deep: bool,
}

impl ValidateArgs {
    fn create_file_readers(
        &self,
        paths: &[std::path::PathBuf],
    ) -> anyhow::Result<HashMap<String, BsxFileReader>> {
        let mut readers = HashMap::new();

        for path in paths {
            let file = File::open(path).map_err(|e| {
                anyhow::anyhow!("Error when opening file {}: {}", path.display(), e)
            })?;

            let reader = BsxFileReader::try_new(file).map_err(|e| {
                anyhow::anyhow!(
                    "Error when creating reader for file {}: {}",
                    path.display(),
                    e
                )
            })?;

            readers.insert(path.to_string_lossy().to_string(), reader);
        }

        Ok(readers)
    }

    fn validate_batch_numbers(
        &self,
        readers: &HashMap<String, BsxFileReader>,
    ) -> anyhow::Result<usize> {
        let batch_numbers = readers
            .iter()
            .map(|(path, reader)| (reader.blocks_total(), path.clone()))
            .into_group_map();

        validate_consistency(&batch_numbers, ValidationErrorType::BatchNumber);

        let common_batch_num = batch_numbers
            .keys()
            .next()
            .copied()
            .ok_or_else(|| anyhow::anyhow!("No files found"))?;

        Ok(common_batch_num)
    }

    fn validate_batches(
        &self,
        iterators: &mut HashMap<
            Rc<String>,
            impl Iterator<Item = PolarsResult<BsxBatch>>,
        >,
        progress_bar: &ProgressBar,
    ) -> anyhow::Result<()> {
        let mut batch_count = 0;

        while let Some(new_batches) = self.read_next_batches(iterators)? {
            self.validate_batch_sizes(&new_batches, batch_count)?;
            self.validate_batch_ranges(&new_batches, batch_count)?;

            if self.deep {
                self.validate_batch_positions(&new_batches, batch_count)?;
            }

            batch_count += 1;
            progress_bar.inc(1);
        }

        Ok(())
    }

    fn read_next_batches(
        &self,
        iterators: &mut HashMap<
            Rc<String>,
            impl Iterator<Item = PolarsResult<BsxBatch>>,
        >,
    ) -> anyhow::Result<Option<Vec<(Rc<String>, BsxBatch)>>> {
        let batches: Result<Vec<_>, _> = iterators
            .iter_mut()
            .map(|(path, iterator)| {
                iterator
                    .next()
                    .map(|batch_result| -> Result<(Rc<String>, BsxBatch), _> {
                        let batch = batch_result.map_err(|e| {
                            anyhow::anyhow!(
                                "Error when reading batch from {}: {}",
                                path,
                                e
                            )
                        })?;
                        Ok::<(Rc<String>, BsxBatch), anyhow::Error>((
                            path.clone(),
                            batch,
                        ))
                    })
                    .transpose()
            })
            .collect();

        match batches? {
            batches if batches.iter().all(|b| b.is_some()) => {
                Ok(Some(batches.into_iter().map(|b| b.unwrap()).collect()))
            },
            batches if batches.iter().all(|b| b.is_none()) => Ok(None),
            _ => Err(anyhow::anyhow!("Files have different number of batches")),
        }
    }

    fn validate_batch_sizes(
        &self,
        batches: &[(Rc<String>, BsxBatch)],
        batch_count: usize,
    ) -> anyhow::Result<()> {
        let batch_sizes = batches
            .iter()
            .map(|(path, batch)| (batch.len(), path))
            .into_group_map();

        validate_consistency(&batch_sizes, ValidationErrorType::BatchSize {
            batch_num: batch_count,
        });

        Ok(())
    }

    fn validate_batch_ranges(
        &self,
        batches: &[(Rc<String>, BsxBatch)],
        batch_count: usize,
    ) -> anyhow::Result<()> {
        let contigs = batches
            .iter()
            .map(|(path, batch)| (batch.as_contig().unwrap(), path))
            .into_group_map();

        validate_consistency(&contigs, ValidationErrorType::BatchRange {
            batch_num: batch_count,
        });

        Ok(())
    }

    fn validate_batch_positions(
        &self,
        batches: &[(Rc<String>, BsxBatch)],
        batch_count: usize,
    ) -> anyhow::Result<()> {
        let positions = batches
            .iter()
            .map(|(path, batch)| {
                (batch.position().to_vec_null_aware().unwrap_left(), path)
            })
            .into_group_map();

        validate_consistency(&positions, ValidationErrorType::BatchPositions {
            batch_num: batch_count,
        });

        Ok(())
    }
}

#[derive(Debug, Clone)]
enum ValidationErrorType {
    BatchNumber,
    BatchSize { batch_num: usize },
    BatchRange { batch_num: usize },
    BatchPositions { batch_num: usize },
}

impl ValidationErrorType {
    fn inconsistency_message(&self) -> String {
        match self {
            ValidationErrorType::BatchNumber => {
                "inconsistency in batch number".to_string()
            },
            ValidationErrorType::BatchSize { batch_num } => {
                format!("inconsistency in batch size. Batch №{}", batch_num)
            },
            ValidationErrorType::BatchRange { batch_num } => {
                format!(
                    "inconsistency in batch values range. Batch №{:?}",
                    batch_num
                )
            },
            ValidationErrorType::BatchPositions { batch_num } => {
                format!("inconsistency in batch positions. Batch №{:?}", batch_num)
            },
        }
    }

    fn most_files_message<T: Debug>(
        &self,
        value: &T,
    ) -> String {
        match self {
            ValidationErrorType::BatchNumber => {
                format!("Most files have {:?} batches", value)
            },
            ValidationErrorType::BatchSize { .. } => {
                format!("Most files have {:?} rows", value)
            },
            ValidationErrorType::BatchRange { .. } => {
                format!("Most files cover {:?} region", value)
            },
            ValidationErrorType::BatchPositions { .. } => {
                format!("")
            },
        }
    }

    fn following_files_message<T: Debug>(
        &self,
        value: &T,
    ) -> String {
        match self {
            ValidationErrorType::BatchNumber => {
                format!("Following files have {:?} batches:", value)
            },
            ValidationErrorType::BatchSize { .. } => {
                format!("Following files have {:?} rows:", value)
            },
            ValidationErrorType::BatchRange { .. } => {
                format!("Following files cover {:?} region:", value)
            },
            ValidationErrorType::BatchPositions { .. } => {
                format!("Following files have different positions")
            },
        }
    }
}

fn validate_consistency<K, V>(
    grouped_data: &HashMap<K, Vec<V>>,
    error_type: ValidationErrorType,
) where
    K: Clone + Debug + Hash + Eq,
    V: Display, {
    if grouped_data.len() > 1 {
        let most_common_key = grouped_data
            .keys()
            .max_by_key(|k| grouped_data.get(k).unwrap().len())
            .unwrap();

        eprintln!(
            "Detected {}",
            style(error_type.inconsistency_message()).red()
        );
        eprintln!(
            "{}",
            style(error_type.most_files_message(most_common_key)).green()
        );

        for (key, values) in grouped_data
            .iter()
            .sorted_by_key(|(_k, v)| v.len())
            .take(grouped_data.len() - 1)
        {
            eprintln!("{}", style(error_type.following_files_message(key)).red());
            for value in values {
                eprintln!("\t{}", value);
            }
        }
        exit(-1);
    }
}

impl PipelineCommand for ValidateArgs {
    fn run(&self) -> anyhow::Result<()> {
        // Expand file paths with wildcards
        let paths = self
            .files
            .iter()
            .flat_map(|path| expand_wildcard(path))
            .flatten()
            .collect::<Vec<_>>();

        // Create readers for all files
        let mut readers = self.create_file_readers(&paths)?;

        // Validate that all files have the same number of batches
        let common_batch_num = self.validate_batch_numbers(&readers)?;

        println!(
            "[{}] All files have {} batches",
            style("V").green(),
            style(common_batch_num).green()
        );

        // Create iterators for batch processing
        let mut iterators = readers
            .iter_mut()
            .map(|(path, reader)| (Rc::new(path.clone()), reader.iter()))
            .collect::<HashMap<_, _>>();

        // Initialize progress bar
        let progress_bar = init_progress(Some(common_batch_num))?;

        // Process and validate each batch
        self.validate_batches(&mut iterators, &progress_bar)?;

        progress_bar.finish_and_clear();

        println!("[{}] All files have equivalent batches", style("V").green());
        println!("{}", style("Files are valid").green());
        Ok(())
    }
}
