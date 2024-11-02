use crate::utils::ipc_reader::IpcFileReader;
use crate::utils::types::{Context, Strand};
use itertools::Itertools;
use polars::prelude::*;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use std::collections::Bound::Excluded;
use std::collections::{BTreeMap, HashMap};
use std::io::{Read, Seek};

#[derive(Eq, Hash, PartialEq, Clone)]
pub struct BSXKey(String, Strand, Context);

pub type BSXIndex = HashMap<BSXKey, BTreeMap<u64, usize>>;

pub struct BSXFile<R: Read + Seek> {
    reader: IpcFileReader<R>,
    index: BSXIndex,
}

impl<R: Read + Seek> BSXFile<R> {
    pub fn new(handle: R, projection: Option<Vec<usize>>) -> Result<Self, PolarsError> {
        let mut reader = IpcFileReader::new(handle, projection, None);
        let mut index: HashMap<BSXKey, BTreeMap<u64, usize>> = HashMap::new();
        for batch_num in 0..reader.metadata().blocks.len() {
            match reader.read_df_at(batch_num) {
                Ok(df) => {
                    // Read first position of the batch
                    let pos = df
                        .column("position")?
                        .u64()?
                        .get(0)
                        .expect("could not get position");
                    // Create key instance for hashmap
                    let key: BSXKey = {
                        let chr_col = df.column("chr")?.categorical()?;
                        BSXKey(
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
                            btree.insert(pos, batch_num);
                        }
                        None => {
                            let mut btree: BTreeMap<u64, usize> = BTreeMap::new();
                            btree.insert(pos, batch_num);
                            index.insert(key, btree);
                        }
                    };
                }
                Err(e) => return Err(e),
            };
        }
        Ok(Self { reader, index })
    }

    pub fn find_index(&self, key: BSXKey, range: (Option<u64>, Option<u64>)) -> Option<Vec<u32>> {
        if let Some(btree) = self.index.get(&key) {
            let (start, end) = (
                range.0.unwrap_or(0),
                range
                    .1
                    .unwrap_or(*btree.iter().max().unwrap().0),
            );

            let mut root = btree.lower_bound(Excluded(&start));
            let mut out: Vec<u32> = Vec::new();
            if let Some(prev) = root.peek_prev() {
                out.push(*prev.1 as u32);
            }
            while let Some(new_batch) = root.next() {
                if *new_batch.0 <= end {
                    out.push(*new_batch.1 as u32);
                }
                else {
                    break;
                };
            }
            Some(out)
        }
        else {
            None
        }
    }

    pub fn add_index(
        &self,
        annotation: &mut DataFrame,
        context: Option<Context>,
    ) -> Result<DataFrame, PolarsError> {
        let mut annotation_colnames = annotation.schema().iter_names();
        for colname in ["chr", "strand", "start", "end"] {
            if !annotation_colnames.contains(colname) {
                return Err(PolarsError::ColumnNotFound(colname.into()));
            }
        }

        let (chr, strand, start, end) = {
            let chr_s = annotation.column("chr")?;
            let strand_s = annotation.column("strand")?;
            let start_s = annotation.column("start")?;
            let end_s = annotation.column("end")?;

            assert_eq!(
                [chr_s, strand_s, start_s, end_s]
                    .iter()
                    .map(|x| { x.null_count() })
                    .all_equal_value()
                    .unwrap(),
                0,
                "Null values in annotation!"
            );
            (
                chr_s
                    .str()?
                    .into_iter()
                    .map(|x| x.unwrap().to_owned())
                    .collect::<Vec<String>>(),
                strand_s
                    .str()?
                    .into_iter()
                    .map(|x| x.unwrap().to_owned())
                    .collect::<Vec<String>>(),
                start_s
                    .u64()?
                    .into_iter()
                    .flatten()
                    .collect::<Vec<u64>>(),
                end_s
                    .u64()?
                    .into_iter()
                    .flatten()
                    .collect::<Vec<u64>>(),
            )
        };
        assert!(
            [chr.len(), strand.len(), start.len(), end.len()]
                .iter()
                .all_equal(),
            "Number of elements in non-null arrays differ!"
        );

        let result: ChunkedArray<ListType> = (&chr, &strand, &start, &end)
            .into_par_iter()
            .map(|(chr_val, strand_val, start_val, end_val)| {
                let key = BSXKey(
                    chr_val.to_owned(),
                    Strand::from_str(strand_val),
                    context.unwrap_or(Context::CG),
                );
                let indices = self._find_index(key, (Some(*start_val), Some(*end_val)));
                indices
                    .map(|indices| UInt32Chunked::from_vec("index".into(), indices).into_series())
            })
            .collect();
        let df = annotation.with_column(result.with_name("index"))?;
        Ok(df)
    }

    pub fn into_iter(self) -> BSXIterator<R> {
        BSXIterator::new(self.reader, self.index)
    }

    pub fn get_batch_num(&self) -> usize {
        self.reader.metadata().blocks.len()
    }
    pub fn get_arrow_schema(&self) -> ArrowSchema {
        self.reader.schema().clone()
    }
    pub fn get_polars_schema(&self) -> Schema {
        let schema = self.get_arrow_schema();
        Schema::from_arrow_schema(&schema)
    }
}

#[derive(Clone)]
pub struct BSXIterator<R: Read + Seek> {
    reader: IpcFileReader<R>,
    index: BSXIndex,
    last_batch: DataFrame,
    last_batch_num: Option<u32>,
}

impl<R: Read + Seek> BSXIterator<R> {
    pub fn new(reader: IpcFileReader<R>, index: BSXIndex) -> BSXIterator<R> {
        Self {
            reader,
            index,
            last_batch: DataFrame::default(),
            last_batch_num: None,
        }
    }
    fn get_last_batch(&mut self) -> DataFrame {
        self.last_batch.clone()
    }
    fn get_last_batch_num(&mut self) -> Option<u32> {
        self.last_batch_num
    }
    pub fn get_batch(&mut self, index: u32) -> Result<DataFrame, PolarsError> {
        if let Some(last_batch_num) = self.get_last_batch_num() {
            if last_batch_num == index {
                return Ok(self.get_last_batch());
            }
        }
        self.last_batch_num = Some(index);
        match self.reader.read_df_at(index as usize) {
            Ok(df) => {
                self.last_batch = df;
                self.get_batch(index)
            }
            Err(e) => Err(e),
        }
    }
    pub fn get_region(
        &mut self,
        indices: &[u32],
        start: u64,
        end: u64,
    ) -> Result<DataFrame, PolarsError> {
        let dfs: Vec<LazyFrame> = indices
            .iter()
            .map(|i| self.get_batch(i.clone()).unwrap())
            .map(|df| {
                df.lazy().filter(
                    col("position")
                        .gt_eq(lit(start))
                        .and(col("position").lt_eq(lit(end))),
                )
            })
            .collect();
        concat(dfs, UnionArgs::default())?.collect()
    }
    pub fn iter_regions(
        &mut self,
        indexed_annotation: &DataFrame,
    ) -> impl Iterator<Item = Result<DataFrame, PolarsError>> {
        let mut schema_cols = indexed_annotation.schema().iter_names();
        for col in ["start", "end", "index"] {
            if !schema_cols.contains(&col.into()) {
                return Err(PolarsError::ColumnNotFound(col.into()));
            }
        }
        let annotation = indexed_annotation.clone();
        let index = annotation
            .column("index")?
            .list()?
            .into_iter()
            .filter_map(|s| match s {
                Some(series) => {
                    let arr: ChunkedArray<UInt32Type> = series.u32()?.to_owned();
                    Some(
                        arr.into_iter()
                            .flatten()
                            .collect::<Vec<u32>>(),
                    )
                }
                _ => None,
            })
            .collect::<Vec<Vec<u32>>>();
        let start = annotation
            .column("start")?
            .u64()?
            .into_iter()
            .flatten()
            .collect::<Vec<u64>>();
        let end = annotation
            .column("end")?
            .u64()?
            .into_iter()
            .flatten()
            .collect::<Vec<u64>>();

        itertools::izip!(&index, &start, &end)
            .map(|(index_val, start_val, end_val)| self.get_region(index_val, start_val, end_val))
    }
}
