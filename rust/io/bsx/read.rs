use crate::data_structs::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use crate::io::bsx::ipc::IpcFileReader;
use polars::error::PolarsResult;
use polars::export::arrow::array::Array;
use polars::export::arrow::record_batch::RecordBatchT;
use polars::frame::DataFrame;
#[cfg(feature = "python")]
use pyo3::{pyclass, pymethods, PyRef, PyResult};
use std::collections::{BTreeMap, HashMap};
use std::io::{Read, Seek};

/// chromosome -> (start_position -> batch_index)
pub type BSXIndex = HashMap<String, BTreeMap<u64, usize>>;

pub struct BsxFileReader<R: Read + Seek> {
    reader: IpcFileReader<R>,
}

impl<R: Read + Seek> BsxFileReader<R> {
    pub fn new(handle: R) -> Self {
        Self {
            reader: IpcFileReader::new(handle, None, None),
        }
    }

    fn process_record_batch(
        &self,
        batch: PolarsResult<RecordBatchT<Box<dyn Array>>>,
    ) -> PolarsResult<EncodedBsxBatch> {
        batch
            .map(|batch| {
                DataFrame::try_from((batch, self.reader.metadata().schema.as_ref()))
                    .expect("Failed to create dataFrame from batch")
            })
            .map(|df| unsafe { EncodedBsxBatch::new_unchecked(df) })
    }

    pub fn get_batch(&mut self, batch_idx: usize) -> Option<PolarsResult<EncodedBsxBatch>> {
        self.reader
            .read_at(batch_idx)
            .map(|res| self.process_record_batch(res))
    }

    pub fn blocks_total(&self) -> usize {
        self.reader.blocks_total()
    }
}

impl<R: Read + Seek> Iterator for BsxFileReader<R> {
    type Item = PolarsResult<EncodedBsxBatch>;
    fn next(&mut self) -> Option<Self::Item> {
        let next = self.reader.next();
        let res = next.map(|res| self.process_record_batch(res));
        if let Some(Ok(data)) = res.as_ref() {
            if data.height() == 0 {
                return None;
            }
        }
        res
    }
}

#[cfg(feature = "python")]
#[pyclass(name = "BsxFileReader")]
pub struct PyBsxFileReader {
    inner: BsxFileReader<File>,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyBsxFileReader {
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let file = File::open(path)?;
        Ok(PyBsxFileReader {
            inner: BsxFileReader::new(file),
        })
    }
    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    pub fn __next__(&mut self) -> Option<PyResult<EncodedBsxBatch>> {
        self.inner.next().map(|res| wrap_polars_result!(res))
    }
}
