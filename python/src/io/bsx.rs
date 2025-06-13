use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use bsxplorer2::io::bsx::{
    BsxFileReader as RsBsxFileReader,
    BsxFileWriter as RsBsxIpcWriter,
};
use polars::prelude::IpcCompression;
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use pyo3_polars::error::PyPolarsErr;

use crate::types::batch::PyBsxBatch;
use crate::utils::{
    FileOrFileLike,
    SinkHandle,
};


#[pyclass(name = "BsxFileReader", unsendable)]
#[derive(Debug, Clone)]
pub struct PyBsxFileReader {
    reader:            RsBsxFileReader,
    current_batch_idx: usize,
}

impl PyBsxFileReader {
    pub(crate) fn into_inner(self) -> RsBsxFileReader {
        self.reader
    }

    pub(crate) fn inner_mut(&mut self) -> &mut RsBsxFileReader {
        &mut self.reader
    }
}

#[pymethods]
impl PyBsxFileReader {
    #[new]
    pub fn new(file: FileOrFileLike) -> PyResult<Self> {
        let reader = match file {
            FileOrFileLike::File(path) => RsBsxFileReader::try_new(File::open(path)?)?,
            FileOrFileLike::ROnlyFileLike(handle) => RsBsxFileReader::try_new(handle)?,
            FileOrFileLike::RWFileLike(handle) => RsBsxFileReader::try_new(handle)?,
        };
        Ok(Self {
            reader,
            current_batch_idx: 0,
        })
    }

    pub fn get_batch(
        &mut self,
        batch_idx: usize,
    ) -> PyResult<Option<PyBsxBatch>> {
        match self.reader.get_batch(batch_idx) {
            Some(Ok(batch)) => Ok(Some(batch.into())),
            Some(Err(e)) => Err(PyPolarsErr::Polars(e).into()),
            None => Ok(None),
        }
    }

    pub fn get_batches(
        &mut self,
        batch_indices: Vec<usize>,
    ) -> PyResult<Vec<Option<PyBsxBatch>>> {
        let results = self.reader.get_batches(&batch_indices);
        let mut py_results = Vec::new();
        for result in results {
            match result {
                Some(Ok(batch)) => py_results.push(Some(batch.into())),
                Some(Err(e)) => return Err(PyPolarsErr::Polars(e).into()),
                None => py_results.push(None),
            }
        }
        Ok(py_results)
    }

    pub fn cache_batches(
        &mut self,
        batch_indices: Vec<usize>,
    ) -> PyResult<()> {
        self.reader
            .cache_batches(&batch_indices)
            .map_err(|e| PyPolarsErr::Polars(e).into())
    }

    pub fn n_threads(&self) -> usize {
        self.reader.n_threads()
    }

    pub fn blocks_total(&self) -> usize {
        self.reader.blocks_total()
    }

    pub fn __iter__(mut slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf.current_batch_idx = 0;
        slf
    }

    pub fn __next__(&mut self) -> Option<PyResult<PyBsxBatch>> {
        if let Some(batch) = self.reader.cache_mut().pop_front() {
            self.current_batch_idx += 1;
            Some(Ok(batch.into()))
        }
        else if self.current_batch_idx < self.reader.blocks_total() {
            let to_read = (self.current_batch_idx
                ..(self.current_batch_idx + self.reader.n_threads()))
                .collect::<Vec<_>>();
            let cache_res = self.reader.cache_batches(&to_read);
            if cache_res.is_ok() {
                self.__next__()
            }
            else {
                Some(Err(PyPolarsErr::Polars(cache_res.unwrap_err()).into()))
            }
        }
        else {
            None
        }
    }
}

#[pyclass(name = "IpcCompression", eq, eq_int)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum PyIpcCompression {
    LZ4,
    #[default]
    ZSTD,
}

impl From<PyIpcCompression> for IpcCompression {
    fn from(value: PyIpcCompression) -> Self {
        match value {
            PyIpcCompression::LZ4 => IpcCompression::LZ4,
            PyIpcCompression::ZSTD => IpcCompression::ZSTD,
        }
    }
}

impl From<IpcCompression> for PyIpcCompression {
    fn from(value: IpcCompression) -> Self {
        match value {
            IpcCompression::LZ4 => PyIpcCompression::LZ4,
            IpcCompression::ZSTD => PyIpcCompression::ZSTD,
        }
    }
}

#[pyclass(name = "BsxFileWriter", unsendable)]
pub struct PyBsxFileWriter {
    writer: Option<RsBsxIpcWriter<BufWriter<Box<dyn SinkHandle>>>>,
}

#[pymethods]
impl PyBsxFileWriter {
    #[new]
    #[pyo3(signature = (sink, chr_names, compression=None))]
    pub fn new(
        sink: FileOrFileLike,
        chr_names: Vec<String>,
        compression: Option<PyIpcCompression>,
    ) -> PyResult<Self> {
        let file: Box<dyn SinkHandle> = sink.get_writer()?;
        let writer = RsBsxIpcWriter::try_new(
            BufWriter::new(file),
            chr_names,
            compression.map(|x| x.into()),
        )
        .map_err(|e| PyIOError::new_err(e.to_string()))?;

        Ok(Self {
            writer: Some(writer),
        })
    }

    #[staticmethod]
    #[pyo3(signature = (sink, fai_path, compression=None))]
    pub fn from_sink_and_fai(
        sink: FileOrFileLike,
        fai_path: PathBuf,
        compression: Option<PyIpcCompression>,
    ) -> PyResult<Self> {
        let file: Box<dyn SinkHandle> = sink.get_writer()?;
        let writer = RsBsxIpcWriter::try_from_sink_and_fai(
            BufWriter::new(file),
            fai_path,
            compression.map(|x| x.into()),
        )
        .map_err(|e| PyIOError::new_err(e.to_string()))?;

        Ok(Self {
            writer: Some(writer),
        })
    }

    #[staticmethod]
    #[pyo3(signature = (sink, fasta_path, compression=None))]
    pub fn from_sink_and_fasta(
        sink: FileOrFileLike,
        fasta_path: PathBuf,
        compression: Option<PyIpcCompression>,
    ) -> PyResult<Self> {
        let file: Box<dyn SinkHandle> = sink.get_writer()?;
        let writer = RsBsxIpcWriter::try_from_sink_and_fasta(
            BufWriter::new(file),
            fasta_path,
            compression.map(|x| x.into()),
        )
        .map_err(|e| PyIOError::new_err(e.to_string()))?;

        Ok(Self {
            writer: Some(writer),
        })
    }

    #[allow(unsafe_code)]
    pub fn write_batch(
        &mut self,
        batch: PyBsxBatch,
    ) -> PyResult<()> {
        self.writer
            .as_mut()
            .ok_or_else(|| PyIOError::new_err("Writer is closed"))?
            .write_batch(batch.into())
            .map_err(|e| PyPolarsErr::Polars(e).into())
    }

    pub fn close(&mut self) -> PyResult<()> {
        if let Some(mut writer) = self.writer.take() {
            writer.close().map_err(|e| PyPolarsErr::Polars(e))?;
        }
        Ok(())
    }

    pub fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    pub fn __exit__(
        &mut self,
        _exc_type: PyObject,
        _exc_value: PyObject,
        _traceback: PyObject,
    ) -> PyResult<()> {
        self.close()
    }
}
