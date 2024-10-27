use crate::utils::types::{Context, Strand};
use polars::export::arrow::io::ipc::read::{read_file_metadata, FileReader};
use polars::io::ArrowReader;
use polars::prelude::PolarsError;
use polars_core::chunked_array::ChunkedArray;
use polars_core::datatypes::{ListType, UInt32Chunked};
use polars_core::frame::DataFrame;
use polars_core::prelude::{IntoSeries, Series};
use std::collections::Bound::Excluded;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;

/// IpcIndex stores information about batches in the IPC file. The key for HashMap is
/// tuple of (chromosome name, [Strand], [Context]), the value is a BTree with batch
/// start position as key and batch number as value.
pub type IpcIndex = HashMap<(String, Strand, Context), BTreeMap<u64, usize>>;

/// Read and create instance of [IpcIndex] from IPC file
pub fn get_index(path: &str) -> Result<IpcIndex, PolarsError> {
    // Open IPC file
    let mut file = File::open(path).expect("could not open ipc file");
    let metadata = read_file_metadata(&mut file).expect("could not read ipc file metadata");
    let mut reader = FileReader::new(file, metadata.clone(), None, None);

    // Init batch counter and index hashmap
    let mut num_batch: usize = 0;
    let mut index: HashMap<(String, Strand, Context), BTreeMap<u64, usize>> = HashMap::new();

    // Iterate over batches from file
    while let Some(batch) = reader.next_record_batch()? {
        let df = DataFrame::try_from((batch, metadata.schema.clone().as_ref()))
            .expect("could not create dataframe");

        // Read first position of the batch
        let pos = df
            .column("position")?
            .u64()?
            .get(0)
            .expect("could not get position");
        // Create key instance for hashmap
        let key: (String, Strand, Context) = {
            let chr_col = df.column("chr")?.categorical()?;
            (
                chr_col
                    .get_rev_map()
                    .get(chr_col.physical().get(0).unwrap())
                    .to_owned(),
                Strand::from_bool(df.column("strand")?.bool()?.get(0)),
                Context::from_bool(df.column("context")?.bool()?.get(0)),
            )
        };
        // If key exists, append batch position
        match index.get_mut(&key) {
            Some(btree) => {
                btree.insert(pos, num_batch);
            }
            None => {
                let mut btree: BTreeMap<u64, usize> = BTreeMap::new();
                btree.insert(pos, num_batch);
                index.insert(key, btree);
            }
        };
        num_batch += 1;
    }
    Ok(index)
}

/// Convert Series of struct (chr, strand, start, end) to column with batch indexes
/// from [IpcIndex]. This is a low-level logic method for [get_index].
pub fn series_to_index(
    series: Series,
    context: Context,
    index_map: &IpcIndex,
) -> Result<Option<Series>, PolarsError> {
    let s = series.struct_()?;
    // Bindings
    let chr_col = &s.field_by_name("chr")?;
    let strand_col = &s.field_by_name("strand")?;
    let start_col = &s.field_by_name("start")?;
    let end_col = &s.field_by_name("end")?;
    // Series
    let chr_arr = chr_col.str()?;
    let strand_arr = strand_col.str()?;
    let start_arr = start_col.u64()?;
    let end_arr = end_col.u64()?;

    let out: ChunkedArray<ListType> = itertools::izip!(
        chr_arr.into_iter(),
        strand_arr.into_iter(),
        start_arr.into_iter(),
        end_arr.into_iter()
    )
    .map(|(opt_chr, opt_strand, opt_start, opt_end)| {
        find_indices(opt_chr, opt_strand, opt_start, opt_end, context, &index_map)
    })
    .collect();
    Ok(Some(Series::from(out)))
}

fn find_indices(
    opt_chr: Option<&str>,
    opt_strand: Option<&str>,
    opt_start: Option<u64>,
    opt_end: Option<u64>,
    context: Context,
    index_map: &IpcIndex,
) -> Option<Series> {
    match (opt_chr, opt_strand, opt_start, opt_end) {
        (Some(chr), Some(strand), Some(start), Some(end)) => {
            let btree = index_map.get(&(chr.to_string(), Strand::from_str(strand), context));
            if let Some(btree) = btree {
                let mut root = btree.lower_bound(Excluded(&start));
                let mut batches = vec![root.peek_prev()];
                while let Some(new_batch) = root.next() {
                    if *new_batch.0 <= end {
                        batches.push(Some(new_batch));
                    }
                    else {
                        break;
                    }
                }
                let idxs = batches
                    .iter()
                    .flatten()
                    .map(|x| *x.1 as u32)
                    .collect::<Vec<u32>>();

                let arr = UInt32Chunked::from_vec("index".into(), idxs);
                Some(arr.into_series())
            }
            else {
                None
            }
        }
        _ => None,
    }
}
