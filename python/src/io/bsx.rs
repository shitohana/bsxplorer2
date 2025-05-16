use std::io::BufWriter;
use std::path::PathBuf;

use bsxplorer2::data_structs::batch::BsxBatch;
use bsxplorer2::io::bsx::{BsxFileReader as RsBsxFileReader,
                          BsxFileWriter as RsBsxIpcWriter};
use polars::prelude::IpcCompression;
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use pyo3_polars::error::PyPolarsErr;

use crate::types::batch::PyBsxBatch;
use crate::types::coords::PyContig;
use crate::types::index::PyBatchIndex;
use crate::utils::{FileOrFileLike, ReadHandle, SinkHandle};


#[pyclass(name = "BsxFileReader", unsendable)]
pub struct PyBsxFileReader {
    reader: RsBsxFileReader<Box<dyn ReadHandle>>,
}

impl PyBsxFileReader {
    pub(crate) fn into_inner(self) -> RsBsxFileReader<Box<dyn ReadHandle>> {
        self.reader
    }
}

#[pymethods]
impl PyBsxFileReader {
    #[new]
    pub fn new(file: FileOrFileLike) -> PyResult<Self> {
        let file = file.get_reader()?;
        let reader = RsBsxFileReader::new(file);
        Ok(Self { reader })
    }

    #[staticmethod]
    pub fn from_file_and_index(
        file: FileOrFileLike,
        index: FileOrFileLike,
    ) -> PyResult<Self> {
        let file = file.get_reader()?;
        let mut index = index.get_reader()?;
        let inner = RsBsxFileReader::from_file_and_index(file, &mut index)?;
        Ok(Self { reader: inner })
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

    pub fn index(&mut self) -> PyResult<PyBatchIndex> {
        self.reader
            .index()
            .cloned()
            .map(PyBatchIndex::from)
            .map_err(|e| e.into())
    }

    pub fn set_index(
        &mut self,
        index: PyBatchIndex,
    ) -> () {
        self.reader.set_index(Some(index.into()))
    }

    pub fn query(
        &mut self,
        query: PyContig,
    ) -> PyResult<Option<PyBsxBatch>> {
        Ok(self.reader.query(&query.into())?.map(PyBsxBatch::from))
    }

    pub fn blocks_total(&self) -> usize {
        self.reader.blocks_total()
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyResult<PyBsxBatch>> {
        slf.reader.next().map(|res| {
            res.map(PyBsxBatch::from)
                .map_err(|e| PyPolarsErr::Polars(e).into())
        })
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
            None,
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
            None,
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
            None,
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
