// Copyright (c) 2020 Ritchie Vink
// Some portions Copyright (c) 2024 NVIDIA CORPORATION & AFFILIATES. All rights
// reserved.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// This is a modified copy of
// polars_core::utils::arrow::io::ipc::read::FileReader with method to select
// specified RecordBatch of an IPC file. GitHub: shitohana

use std::default::Default;
use std::io::{Read, Seek};

use polars::export::arrow::array::Array;
use polars::export::arrow::io::ipc::read::{
    read_batch, read_file_dictionaries, read_file_metadata, Dictionaries,
    FileMetadata,
};
use polars::export::arrow::record_batch::RecordBatchT;
use polars::prelude::{
    ArrowSchema, DataFrame, PlHashMap, PolarsError, PolarsResult,
};

fn apply_projection(
    chunk: RecordBatchT<Box<dyn Array>>,
    map: &PlHashMap<usize, usize>,
) -> RecordBatchT<Box<dyn Array>> {
    // re-order according to projection
    let arrays = chunk.into_arrays();
    let mut new_arrays = arrays.clone();

    map.iter()
        .for_each(|(old, new)| new_arrays[*new] = arrays[*old].clone());

    RecordBatchT::new(Default::default(), new_arrays)
}

pub fn prepare_projection(
    schema: &ArrowSchema,
    mut projection: Vec<usize>,
) -> (Vec<usize>, PlHashMap<usize, usize>, ArrowSchema) {
    let schema = projection
        .iter()
        .map(|x| {
            let (k, v) = schema.get_at_index(*x).unwrap();
            (k.clone(), v.clone())
        })
        .collect();

    let mut indices = (0..projection.len()).collect::<Vec<_>>();
    indices.sort_unstable_by_key(|&i| &projection[i]);
    let map = indices
        .iter()
        .copied()
        .enumerate()
        .fold(PlHashMap::default(), |mut acc, (index, new_index)| {
            acc.insert(index, new_index);
            acc
        });
    projection.sort_unstable();

    // check unique
    if !projection.is_empty() {
        let mut previous = projection[0];

        for &i in &projection[1..] {
            assert!(
                previous < i,
                "The projection on IPC must not contain duplicates"
            );
            previous = i;
        }
    }

    (projection, map, schema)
}

/// An iterator of [`RecordBatchT`]s from an Arrow IPC file.
pub struct IpcFileReader<R: Read + Seek> {
    reader: R,
    metadata: FileMetadata,
    // the dictionaries are going to be read
    dictionaries: Option<Dictionaries>,
    current_block: usize,
    blocks_total: usize,
    projection: Option<(Vec<usize>, PlHashMap<usize, usize>, ArrowSchema)>,
    remaining: usize,
    data_scratch: Vec<u8>,
    message_scratch: Vec<u8>,
}

impl<R: Read + Seek> IpcFileReader<R> {
    /// Creates a new [`IpcFileReader`]. Use `projection` to only take certain
    /// columns. # Panic
    /// Panics iff the projection is not in increasing order (e.g. `[1, 0]` nor
    /// `[0, 1, 1]` are valid)
    pub fn new(
        mut reader: R,
        projection: Option<Vec<usize>>,
        limit: Option<usize>,
    ) -> Self {
        let metadata = read_file_metadata(&mut reader)
            .expect("could not read ipc metadata");
        let projection = projection.map(|projection| {
            let (p, h, schema) =
                prepare_projection(&metadata.schema, projection);
            (p, h, schema)
        });
        let blocks_total = metadata.blocks.len();
        Self {
            reader,
            metadata,
            dictionaries: Default::default(),
            projection,
            remaining: limit.unwrap_or(usize::MAX),
            current_block: 0,
            blocks_total,
            data_scratch: Default::default(),
            message_scratch: Default::default(),
        }
    }

    /// Return the schema of the file
    pub fn schema(&self) -> &ArrowSchema {
        self.projection
            .as_ref()
            .map(|x| &x.2)
            .unwrap_or(&self.metadata.schema)
    }

    pub fn blocks_total(&self) -> usize {
        self.blocks_total
    }

    /// Returns the [`FileMetadata`]
    pub fn metadata(&self) -> &FileMetadata {
        &self.metadata
    }

    /// Consumes this FileReader, returning the underlying reader
    pub fn into_inner(self) -> R {
        self.reader
    }

    /// Get the inner memory scratches so they can be reused in a new writer.
    /// This can be utilized to save memory allocations for performance reasons.
    pub fn get_scratches(&mut self) -> (Vec<u8>, Vec<u8>) {
        (
            std::mem::take(&mut self.data_scratch),
            std::mem::take(&mut self.message_scratch),
        )
    }

    /// Set the inner memory scratches so they can be reused in a new writer.
    /// This can be utilized to save memory allocations for performance reasons.
    pub fn set_scratches(
        &mut self,
        scratches: (Vec<u8>, Vec<u8>),
    ) {
        (self.data_scratch, self.message_scratch) = scratches;
    }

    fn read_dictionaries(&mut self) -> PolarsResult<()> {
        if self.dictionaries.is_none() {
            self.dictionaries = Some(read_file_dictionaries(
                &mut self.reader,
                &self.metadata,
                &mut self.data_scratch,
            )?);
        };
        Ok(())
    }

    pub fn read_at(
        &mut self,
        num: usize,
    ) -> Option<PolarsResult<RecordBatchT<Box<dyn Array>>>> {
        if num >= self.metadata.blocks.len() {
            return None;
        }
        match self.read_dictionaries() {
            Ok(_) => {},
            Err(e) => return Some(Err(e)),
        };
        let chunk = read_batch(
            &mut self.reader,
            self.dictionaries.as_ref().unwrap(),
            &self.metadata,
            self.projection
                .as_ref()
                .map(|x| x.0.as_ref()),
            Some(self.remaining),
            num,
            &mut self.message_scratch,
            &mut self.data_scratch,
        );
        self.remaining -= chunk
            .as_ref()
            .map(|x| x.len())
            .unwrap_or_default();
        let chunk = if let Some((_, map, _)) = &self.projection {
            // re-order according to projection
            chunk.map(|chunk| apply_projection(chunk, map))
        } else {
            chunk
        };
        Some(chunk)
    }

    pub fn read_df_at(
        &mut self,
        num: usize,
    ) -> Result<DataFrame, PolarsError> {
        let batch = self
            .read_at(num)
            .expect("Failed to read");
        let result =
            DataFrame::try_from((batch?, self.metadata.schema.as_ref()));
        result
    }
}

impl<R: Read + Seek> Iterator for IpcFileReader<R> {
    type Item = PolarsResult<RecordBatchT<Box<dyn Array>>>;

    fn next(&mut self) -> Option<Self::Item> {
        // get current block
        if self.current_block == self.metadata.blocks.len() {
            return None;
        }

        match self.read_dictionaries() {
            Ok(_) => {},
            Err(e) => return Some(Err(e)),
        };

        let block = self.current_block;
        self.current_block += 1;

        let chunk = read_batch(
            &mut self.reader,
            self.dictionaries.as_ref().unwrap(),
            &self.metadata,
            self.projection
                .as_ref()
                .map(|x| x.0.as_ref()),
            Some(self.remaining),
            block,
            &mut self.message_scratch,
            &mut self.data_scratch,
        );
        self.remaining -= chunk
            .as_ref()
            .map(|x| x.len())
            .unwrap_or_default();

        let chunk = if let Some((_, map, _)) = &self.projection {
            // re-order according to projection
            chunk.map(|chunk| apply_projection(chunk, map))
        } else {
            chunk
        };
        Some(chunk)
    }
}
