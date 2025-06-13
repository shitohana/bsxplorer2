use std::collections::VecDeque;
use std::io::Cursor;
use std::os::fd::AsRawFd;
use std::sync::Arc;

use itertools::Itertools;
use memmap2::{
    Mmap,
    MmapOptions,
};
use polars::prelude::*;
use polars_arrow::io::ipc::read::{
    read_batch,
    read_file_dictionaries,
    read_file_metadata,
    Dictionaries,
    FileMetadata,
};
use rayon::prelude::*;

use crate::data_structs::batch::BsxBatch;
use crate::getter_fn;
use crate::utils::THREAD_POOL;

// ThreadLocalHandle no longer needs a lifetime parameter as it holds its own
// Arc<Mmap>
struct ThreadLocalHandle {
    _mmap:           Arc<Mmap>, // Hold a clone of the shared Mmap Arc
    handle:          Cursor<&'static [u8]>, /* Cursor borrowing from the Mmap inside
                                 * the Arc */
    data_scratch:    Vec<u8>,
    message_scratch: Vec<u8>,
}

impl ThreadLocalHandle {
    // Create a new handle from a clone of the shared mmap Arc
    pub fn new(mmap: Arc<Mmap>) -> Self {
        let slice: &[u8] = mmap.as_ref();
        let static_slice: &'static [u8] = unsafe { std::mem::transmute(slice) };
        let handle = Cursor::new(static_slice);

        Self {
            _mmap: mmap, // Store the Arc clone, ensuring the Mmap lives long enough
            handle,      // Store the cursor borrowing from the 'static transmuted slice
            data_scratch: Vec::new(),
            message_scratch: Vec::new(),
        }
    }

    pub fn read_batch(
        &mut self,
        index: usize,
        dictionaries: &Dictionaries,
        metadata: &FileMetadata,
    ) -> Option<PolarsResult<DataFrame>> {
        if index >= metadata.blocks.len() {
            return None;
        }
        let chunk = read_batch(
            &mut self.handle,
            dictionaries,
            metadata,
            None,
            None,
            index,
            &mut self.message_scratch,
            &mut self.data_scratch,
        );
        if let Err(err) = chunk {
            Some(Err(err))
        }
        else {
            let chunk = chunk.unwrap();
            let result = DataFrame::try_from((chunk, metadata.schema.as_ref()));
            Some(result)
        }
    }
}

/// A reader for .bsx files based on Apache Arrow IPC format, optimized for
/// reading batches in parallel.
pub struct BsxFileReader {
    thread_local_handles: Vec<ThreadLocalHandle>, // Holds the new struct
    cache:                VecDeque<BsxBatch>,
    metadata:             Arc<FileMetadata>,
    dictionaries:         Arc<Dictionaries>,
    blocks_total:         usize,
    _mmap:                Arc<Mmap>, /* Still need to hold the original Arc to keep
                                      * the Mmap alive */
}

impl Clone for BsxFileReader {
    fn clone(&self) -> Self {
        let thread_local_handles: Vec<ThreadLocalHandle> = (0..self.n_threads())
            .map(|_| ThreadLocalHandle::new(self._mmap.clone())) // Pass a clone of the Arc to each handle
            .collect_vec();

        Self {
            thread_local_handles,
            cache: self.cache.clone(),
            metadata: self.metadata.clone(),
            dictionaries: self.dictionaries.clone(),
            blocks_total: self.blocks_total,
            _mmap: self._mmap.clone(),
        }
    }
}

impl std::fmt::Debug for BsxFileReader {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.debug_struct("BsxFileReader")
            .field("cache", &self.cache)
            .field("metadata", &self.metadata)
            .field("dictionaries", &self.dictionaries)
            .field("blocks_total", &self.blocks_total)
            .finish()
    }
}

impl BsxFileReader {
    getter_fn!(cache, mut VecDeque<BsxBatch>);

    /// Creates a new `BsxFileReader` from a file handle.
    ///
    /// This will memory map the file and read its metadata and dictionaries.
    /// It also initializes a thread pool and thread-local handles for parallel
    /// reading.
    pub fn try_new<R: AsRawFd + 'static>(handle: R) -> anyhow::Result<Self> {
        let n_threads = THREAD_POOL.current_num_threads();
        Self::with_n_threads(handle, n_threads)
    }

    pub fn with_n_threads<R: AsRawFd + 'static>(
        handle: R,
        n_threads: usize,
    ) -> anyhow::Result<Self> {
        // 1. Create the shared Arc<Mmap>
        let mmap = Arc::new(unsafe { MmapOptions::new().map(&handle)? });
        // Create a temporary reader to read metadata and dictionaries
        // This reader borrows directly from the shared mmap Arc
        let mut reader = Cursor::new(mmap.as_ref());
        let metadata = read_file_metadata(&mut reader)?;
        let dictionaries = read_file_dictionaries(
            &mut reader,
            &metadata,
            &mut Vec::new(), // Use a new scratch vector for this part
        )?;

        let blocks_total = metadata.blocks.len();

        // 2. Create thread-local handles, each with a clone of the Arc<Mmap>
        let thread_local_handles: Vec<ThreadLocalHandle> = (0..n_threads)
            .map(|_| ThreadLocalHandle::new(mmap.clone())) // Pass a clone of the Arc to each handle
            .collect_vec();


        Ok(Self {
            thread_local_handles,
            cache: VecDeque::new(),
            metadata: Arc::new(metadata),
            dictionaries: Arc::new(dictionaries),
            blocks_total,
            _mmap: mmap, // Store the original Arc in the reader
        })
    }

    /// Reads a single batch from the file at the given index.
    pub fn get_batch(
        &mut self,
        batch_idx: usize,
    ) -> Option<PolarsResult<BsxBatch>> {
        let handle = &mut self.thread_local_handles[0];
        let df = handle.read_batch(
            batch_idx,
            self.dictionaries.as_ref(),
            self.metadata.as_ref(),
        );
        if let Some(Ok(df)) = df {
            Some(Ok(unsafe { BsxBatch::new_unchecked(df) }))
        }
        else {
            df.map(|e| Err(e.unwrap_err()))
        }
    }

    /// Reads multiple batches from the file at the given indices in parallel.
    pub fn get_batches(
        &mut self,
        batch_indices: &[usize],
    ) -> Vec<Option<PolarsResult<BsxBatch>>> {
        self.iter_batches(batch_indices.iter().cloned())
            .collect_vec()
    }

    /// Returns an iterator that reads batches from the file at the given
    /// indices in parallel.
    pub fn iter_batches<'a, I: IntoIterator<Item = usize> + 'a>(
        &'a mut self,
        batch_indices: I,
    ) -> impl Iterator<Item = Option<PolarsResult<BsxBatch>>> + 'a {
        let chunks_to_read = batch_indices
            .into_iter()
            .chunks(self.n_threads())
            .into_iter()
            .map(|c| c.collect_vec())
            .collect_vec();

        let iter = THREAD_POOL.install(|| {
            chunks_to_read
                .into_iter()
                .flat_map(|tasks_chunk| {
                    // Number of tasks in chunk is <= number of threads
                    self.thread_local_handles
                        .par_iter_mut()
                        .zip(tasks_chunk.par_iter())
                        .map(|(handle, idx)| {
                            handle.read_batch(*idx, &self.dictionaries, &self.metadata)
                        })
                        .collect::<Vec<_>>()
                })
                .map(|df_res| {
                    df_res.map(|df_res| {
                        df_res.map(|df| unsafe { BsxBatch::new_unchecked(df) })
                    })
                })
        });
        iter
    }

    /// Returns an iterator that reads all batches from the file sequentially.
    pub fn iter(&mut self) -> impl Iterator<Item = PolarsResult<BsxBatch>> + '_ {
        self.iter_batches(0..self.blocks_total())
            .map(|i| i.unwrap())
    }

    /// Returns the total number of data batches (blocks) in the file.
    pub fn blocks_total(&self) -> usize {
        self.blocks_total
    }

    /// Reads multiple batches from the file at the given indices in parallel
    /// and adds them to the internal cache.
    pub fn cache_batches(
        &mut self,
        batch_indices: &[usize],
    ) -> PolarsResult<()> {
        let mut new_batches = self
            .get_batches(batch_indices)
            .into_iter()
            .flatten()
            .collect::<PolarsResult<VecDeque<_>>>()?;
        self.cache.append(&mut new_batches);
        Ok(())
    }

    /// Returns the number of threads used by the internal thread pool.
    pub fn n_threads(&self) -> usize {
        self.thread_local_handles.len()
    }
}

impl IntoIterator for BsxFileReader {
    type IntoIter = BsxFileIterator;
    type Item = PolarsResult<BsxBatch>;

    fn into_iter(self) -> Self::IntoIter {
        BsxFileIterator {
            reader:        self,
            current_batch: 0,
        }
    }
}

/// An iterator over the batches in a `BsxFileReader`.
pub struct BsxFileIterator {
    reader:        BsxFileReader,
    current_batch: usize,
}

impl Iterator for BsxFileIterator {
    type Item = PolarsResult<BsxBatch>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(batch) = self.reader.cache_mut().pop_front() {
            self.current_batch += 1;
            Some(Ok(batch))
        }
        else if self.current_batch < self.reader.blocks_total() {
            let to_read = (self.current_batch
                ..(self.current_batch + self.reader.n_threads()))
                .collect_vec();
            let cache_res = self.reader.cache_batches(&to_read);

            if let Err(err) = cache_res {
                Some(Err(err))
            }
            else {
                self.next()
            }
        }
        else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let blocks = self.reader.blocks_total();
        (blocks, Some(blocks))
    }

    fn count(self) -> usize
    where
        Self: Sized, {
        self.reader.blocks_total() - self.current_batch
    }

    fn last(self) -> Option<Self::Item>
    where
        Self: Sized, {
        let mut reader = self.reader;
        reader.get_batch(reader.blocks_total() - 1)
    }

    fn nth(
        &mut self,
        n: usize,
    ) -> Option<Self::Item> {
        self.current_batch = n;
        self.reader.cache_mut().clear();
        self.reader.get_batch(n)
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::path::PathBuf;

    use rstest::{
        fixture,
        rstest,
    };

    use super::*;

    #[fixture]
    fn reader() -> BsxFileReader {
        BsxFileReader::try_new(
            File::open(
                PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/report.bsx"),
            )
            .expect("Error opening test report file"),
        )
        .unwrap()
    }

    #[rstest]
    fn create_reader(_reader: BsxFileReader) {}

    #[rstest]
    fn test_reading(mut _reader: BsxFileReader) -> anyhow::Result<()> {
        _reader.cache_batches(&[1, 2, 3, 4, 5, 6])?;
        assert!(!_reader.cache_mut().pop_front().unwrap().is_empty());
        Ok(())
    }

    #[rstest]
    fn test_iter(reader: BsxFileReader) {
        let mut batch_count = 0;
        let blocks_total = reader.blocks_total();
        for batch in reader.into_iter() {
            assert!(batch.is_ok());
            batch_count += 1;
        }
        assert_eq!(batch_count, blocks_total)
    }
}
