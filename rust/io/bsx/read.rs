use crate::data_structs::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use crate::io::bsx::ipc::IpcFileReader;
use log::{debug, trace, warn};
use polars::error::PolarsResult;
use polars::export::arrow::array::Array;
use polars::export::arrow::record_batch::RecordBatchT;
use polars::frame::DataFrame;
#[cfg(feature = "python")]
use pyo3::{pyclass, pymethods, PyRef, PyResult};
use std::collections::{BTreeMap, HashMap};
use std::io::{Read, Seek};

/// Chromosome indexing structure for BSX files
///
/// Maps chromosome names to a tree of positions, where each position points to
/// the corresponding batch index in the BSX file.
/// Structure: chromosome -> (start_position -> batch_index)
pub type BSXIndex = HashMap<String, BTreeMap<u64, usize>>;

/// Reader for BSX files
///
/// Provides functionality to read and process BSX files containing genomic data.
/// The reader implements the Iterator trait for sequential batch processing.
pub struct BsxFileReader<R: Read + Seek> {
    /// The underlying IPC file reader that handles the low-level Arrow format
    reader: IpcFileReader<R>,
}

impl<R: Read + Seek> BsxFileReader<R> {
    /// Creates a new BSX file reader from a seekable read source
    ///
    /// # Arguments
    ///
    /// * `handle` - Any type that implements the Read and Seek traits
    ///
    /// # Returns
    ///
    /// A new BsxFileReader instance
    pub fn new(handle: R) -> Self {
        debug!("Creating new BsxFileReader");
        Self {
            reader: IpcFileReader::new(handle, None, None),
        }
    }

    /// Processes a record batch into an EncodedBsxBatch
    ///
    /// # Arguments
    ///
    /// * `batch` - The raw record batch from the IPC reader
    ///
    /// # Returns
    ///
    /// A processed EncodedBsxBatch wrapped in a PolarsResult
    fn process_record_batch(
        &self,
        batch: PolarsResult<RecordBatchT<Box<dyn Array>>>,
    ) -> PolarsResult<EncodedBsxBatch> {
        trace!("Processing record batch");
        batch
            .map(|batch| {
                DataFrame::try_from((batch, self.reader.metadata().schema.as_ref()))
                    .expect("Failed to create DataFrame from batch - schema mismatch")
            })
            .map(|df| {
                trace!(
                    "Creating EncodedBsxBatch from DataFrame with {} rows",
                    df.height()
                );
                unsafe { EncodedBsxBatch::new_unchecked(df) }
            })
    }

    /// Retrieves a specific batch by index
    ///
    /// # Arguments
    ///
    /// * `batch_idx` - The index of the batch to retrieve
    ///
    /// # Returns
    ///
    /// The processed batch if it exists, None otherwise
    pub fn get_batch(&mut self, batch_idx: usize) -> Option<PolarsResult<EncodedBsxBatch>> {
        debug!("Getting batch at index {}", batch_idx);
        self.reader
            .read_at(batch_idx)
            .map(|res| self.process_record_batch(res))
    }

    /// Returns the total number of blocks in the BSX file
    ///
    /// # Returns
    ///
    /// The count of blocks (batches) in the file
    pub fn blocks_total(&self) -> usize {
        let count = self.reader.blocks_total();
        debug!("Total blocks in BSX file: {}", count);
        count
    }
}

impl<R: Read + Seek> Iterator for BsxFileReader<R> {
    type Item = PolarsResult<EncodedBsxBatch>;

    /// Advances the iterator and returns the next batch
    ///
    /// Returns None when there are no more batches or when an empty batch is encountered.
    fn next(&mut self) -> Option<Self::Item> {
        trace!("Getting next batch from BsxFileReader");
        let next = self.reader.next();
        let res = next.map(|res| self.process_record_batch(res));

        if let Some(Ok(data)) = res.as_ref() {
            if data.height() == 0 {
                debug!("Reached empty batch, stopping iteration");
                return None;
            }
        }

        if res.is_none() {
            debug!("No more batches available in BSX file");
        }

        res
    }
}

#[cfg(feature = "python")]
#[pyclass(name = "BsxFileReader")]
/// Python binding for the BsxFileReader
pub struct PyBsxFileReader {
    /// The wrapped Rust BsxFileReader implementation
    inner: BsxFileReader<File>,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyBsxFileReader {
    /// Creates a new PyBsxFileReader from a file path
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the BSX file
    ///
    /// # Returns
    ///
    /// A new PyBsxFileReader instance wrapped in a PyResult
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        debug!("Creating new PyBsxFileReader for file: {}", path);
        let file = File::open(&path)
            .with_context(|| format!("Failed to open BSX file at path: {}", path))
            .map_err(|e| {
                warn!("Failed to open BSX file: {}", e);
                pyo3::exceptions::PyIOError::new_err(e.to_string())
            })?;

        Ok(PyBsxFileReader {
            inner: BsxFileReader::new(file),
        })
    }

    /// Iterator protocol support for Python
    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        trace!("PyBsxFileReader.__iter__ called");
        slf
    }

    /// Advances the iterator and returns the next batch
    pub fn __next__(&mut self) -> Option<PyResult<EncodedBsxBatch>> {
        trace!("PyBsxFileReader.__next__ called");
        self.inner.next().map(|res| wrap_polars_result!(res))
    }
}
