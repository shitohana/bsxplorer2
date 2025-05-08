use std::io::{BufReader, BufWriter};
use std::path::PathBuf;

use bsxplorer2::data_structs::batch::{BsxBatch,
                                      BsxBatchMethods,
                                      EncodedBsxBatch};
use bsxplorer2::data_structs::coords::Contig;
use bsxplorer2::io::bsx::{BatchIndex,
                          BsxFileReader as RsBsxFileReader,
                          BsxIpcWriter as RsBsxIpcWriter};
use polars::prelude::IpcCompression;
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use pyo3_file::PyFileLikeObject;
use pyo3_polars::error::PyPolarsErr;
use pyo3_polars::PyDataFrame;

use crate::types::coords::PyContig;

#[pyclass]
pub struct PyBatchIndex {
    inner: BatchIndex<String, u32>,
}

#[pymethods]
impl PyBatchIndex {
    #[new]
    pub fn new() -> Self {
        Self {
            inner: BatchIndex::new(),
        }
    }

    pub fn insert(
        &mut self,
        contig: PyContig, // Using &PyAny to allow extraction
        batch_idx: usize,
    ) -> PyResult<()> {
        self.inner
            .insert((&contig).into(), batch_idx);
        Ok(())
    }

    pub fn sort(
        &self,
        contigs: Vec<PyContig>,
    ) -> Vec<PyContig> {
        let contigs_inner: Vec<Contig<String, u32>> = contigs
            .iter()
            .map(Contig::from)
            .collect();
        let sorted_inner = self
            .inner
            .sort(contigs_inner.into_iter());
        sorted_inner
            .into_iter()
            .map(PyContig::from)
            .collect()
    }

    pub fn find(
        &self,
        contig: PyContig,
    ) -> Option<Vec<usize>> {
        self.inner.find(&Contig::from(&contig))
    }

    pub fn chr_order(&self) -> Vec<String> {
        self.inner
            .get_chr_order()
            .iter()
            .cloned()
            .collect()
    }
}

#[pyclass(name = "BsxFileReader")]
pub struct PyBsxFileReader {
    reader: RsBsxFileReader<BufReader<PyFileLikeObject>>,
}

impl PyBsxFileReader {
    pub(crate) fn into_inner(
        self
    ) -> RsBsxFileReader<BufReader<PyFileLikeObject>> {
        self.reader
    }
}

#[pymethods]
impl PyBsxFileReader {
    #[new]
    pub fn new(file: PyObject) -> PyResult<Self> {
        let file_like = PyFileLikeObject::with_requirements(
            file, true, true, false, false,
        )?;
        let reader = RsBsxFileReader::new(BufReader::new(file_like));
        Ok(Self { reader })
    }

    pub fn get_batch(
        &mut self,
        batch_idx: usize,
    ) -> PyResult<Option<PyDataFrame>> {
        match self.reader.get_batch(batch_idx) {
            Some(Ok(batch)) => Ok(Some(PyDataFrame(batch.take()))),
            Some(Err(e)) => Err(PyPolarsErr::Polars(e).into()),
            None => Ok(None),
        }
    }

    pub fn blocks_total(&self) -> usize { self.reader.blocks_total() }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> { slf }

    pub fn __next__(
        mut slf: PyRefMut<'_, Self>
    ) -> Option<PyResult<PyDataFrame>> {
        slf.reader.next().map(|res| {
            res.map(|batch| PyDataFrame(batch.take()))
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

#[pyclass(name = "BsxIpcWriter")]
pub struct PyBsxFileWriter {
    writer: Option<RsBsxIpcWriter<BufWriter<PyFileLikeObject>>>,
}

#[pymethods]
impl PyBsxFileWriter {
    #[new]
    #[pyo3(signature = (sink, chr_names, compression=None))]
    pub fn new(
        sink: PyObject,
        chr_names: Vec<String>,
        compression: Option<PyIpcCompression>,
    ) -> PyResult<Self> {
        let file_like = PyFileLikeObject::with_requirements(
            sink, false, false, true, false,
        )?;

        let writer = RsBsxIpcWriter::try_new(
            BufWriter::new(file_like),
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
        sink: PyObject,
        fai_path: PathBuf,
        compression: Option<PyIpcCompression>,
    ) -> PyResult<Self> {
        let file_like = PyFileLikeObject::with_requirements(
            sink, false, false, true, false,
        )?;

        let writer = RsBsxIpcWriter::try_from_sink_and_fai(
            BufWriter::new(file_like),
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
        sink: PyObject,
        fasta_path: PathBuf,
        compression: Option<PyIpcCompression>,
    ) -> PyResult<Self> {
        let file_like = PyFileLikeObject::with_requirements(
            sink, false, false, true, false,
        )?;

        let writer = RsBsxIpcWriter::try_from_sink_and_fasta(
            BufWriter::new(file_like),
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
    pub fn write_encoded_batch(
        &mut self,
        batch: PyDataFrame,
    ) -> PyResult<()> {
        let encoded_batch = unsafe { EncodedBsxBatch::new_unchecked(batch.0) };
        self.writer
            .as_mut()
            .ok_or_else(|| PyIOError::new_err("Writer is closed"))?
            .write_encoded_batch(encoded_batch)
            .map_err(|e| PyPolarsErr::Polars(e).into())
    }

    #[allow(unsafe_code)]
    pub fn write_batch(
        &mut self,
        batch: PyDataFrame,
    ) -> PyResult<()> {
        let bsx_batch = unsafe { BsxBatch::new_unchecked(batch.0) };
        self.writer
            .as_mut()
            .ok_or_else(|| PyIOError::new_err("Writer is closed"))?
            .write_batch(bsx_batch)
            .map_err(|e| PyPolarsErr::Polars(e).into())
    }

    pub fn close(&mut self) -> PyResult<()> {
        if let Some(mut writer) = self.writer.take() {
            writer
                .close()
                .map_err(|e| PyPolarsErr::Polars(e))?;
        }
        Ok(())
    }

    pub fn __enter__(slf: Py<Self>) -> Py<Self> { slf }

    pub fn __exit__(
        &mut self,
        _exc_type: PyObject,
        _exc_value: PyObject,
        _traceback: PyObject,
    ) -> PyResult<()> {
        self.close()
    }
}
