use crate::utils::ipc_reader::IpcFileReader as IPCFileReader;
use polars::prelude::*;
use polars_core::prelude::{ChunkedArray, DataFrame, ListType, Series};
use std::fs::File;
use std::ops::{BitOr, Not};

pub fn series_to_offsets(series: Series, ipc_path: &str) -> Result<Option<Series>, PolarsError> {
    let s = series.struct_()?;
    let indices_s = s.field_by_name("index")?;
    let starts_s = s.field_by_name("start")?;
    let ends_s = s.field_by_name("end")?;

    let indices_arr = indices_s.u32()?;
    let starts_arr = starts_s.u64()?;
    let ends_arr = ends_s.u64()?;

    let file = File::open(ipc_path).expect("could not open ipc path");
    let mut ipc_reader = IPCFileReader::new(file, None, None);
    let mut cur_batch_num = u32::MAX;
    let mut cur_batch: DataFrame = DataFrame::empty();

    let offsets: ChunkedArray<ListType> = itertools::izip!(
        starts_arr.into_iter(),
        ends_arr.into_iter(),
        indices_arr.into_iter()
    )
    .map(
        |(opt_start, opt_end, opt_index)| match (opt_start, opt_end, opt_index) {
            (Some(start), Some(end), Some(index)) => {
                if index != cur_batch_num {
                    cur_batch_num = index;
                    cur_batch = ipc_reader
                        .read_df_at(index as usize)
                        .expect("could not read df at index");
                }
                let pos = cur_batch.column("position").unwrap();
                let bin_range = pos
                    .lt(start)
                    .unwrap()
                    .bitor(pos.gt(end).unwrap())
                    .not();
                let (mut start_pos, mut end_pos) = (None, None);
                for (idx, val) in bin_range.into_iter().enumerate() {
                    match val {
                        Some(true) => {
                            if start_pos.is_none() {
                                start_pos = Some(idx as i32);
                            }
                        }
                        Some(false) => {
                            if start_pos.is_some() && end_pos.is_none() {
                                end_pos = Some(idx as i32);
                                break;
                            }
                        }
                        _ => continue,
                    }
                }

                if start_pos.is_none() && end_pos.is_none() {
                    return None;
                }

                let offsets = vec![
                    start_pos.unwrap_or(0),
                    end_pos.unwrap_or(bin_range.len() as i32),
                ];
                let arr = Int32Chunked::from_vec("offset".into(), offsets);
                Some(arr.into_series())
            }
            _ => None,
        },
    )
    .collect();
    Ok(Some(Series::from(offsets)))
}
