use std::io::{BufReader, BufWriter};
use std::path::PathBuf;

use bsxplorer2::data_structs::batch::{BsxBatch, BsxBatchMethods, EncodedBsxBatch};
use bsxplorer2::io::bsx::{BsxFileReader as RsBsxFileReader, BsxIpcWriter as RsBsxIpcWriter};
use bsxplorer2::exports::polars::prelude::{IpcCompression};
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use pyo3_file::PyFileLikeObject;
use pyo3_polars::error::PyPolarsErr;
use pyo3_polars::PyDataFrame;

/// Reader for BSX files.
///
/// Parameters
/// ----------
/// file : str or file-like object
///     Path to the BSX file or a file-like object.
#[pyclass(name="BsxFileReader")]
pub struct PyBsxFileReader {
    reader: RsBsxFileReader<BufReader<PyFileLikeObject>>,
}

#[pymethods]
impl PyBsxFileReader {
    #[new]
    pub fn new(file: PyObject) -> PyResult<Self> {
        let file_like = PyFileLikeObject::with_requirements(file, true, true, false, false)?;
        let reader = RsBsxFileReader::new(BufReader::new(file_like));
        Ok(Self { reader })
    }

    /// Get a specific batch by index.
    ///
    /// Parameters
    /// ----------
    /// batch_idx : int
    ///     The index of the batch to retrieve.
    ///
    /// Returns
    /// -------
    /// PyDataFrame or None
    ///     The requested batch as a Polars DataFrame, or None if the index is out of bounds.
    pub fn get_batch(&mut self, batch_idx: usize) -> PyResult<Option<PyDataFrame>> {
        match self.reader.get_batch(batch_idx) {
            Some(Ok(batch)) => Ok(Some(PyDataFrame(batch.take()))),
            Some(Err(e)) => Err(PyPolarsErr::Polars(e).into()),
            None => Ok(None),
        }
    }

    /// Get the total number of batches in the file.
    ///
    /// Returns
    /// -------
    /// int
    ///     The total number of batches.
    pub fn blocks_total(&self) -> usize {
        self.reader.blocks_total()
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyResult<PyDataFrame>> {
        slf.reader
            .next()
            .map(|res| res.map(|batch| PyDataFrame(batch.take())).map_err(|e| PyPolarsErr::Polars(e).into()))
    }
}

#[pyclass(name="IpcCompression", eq, eq_int)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum PyIpcCompression {
    /// LZ4 (framed)
    LZ4,
    /// ZSTD
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

/// Writer for BSX data in Arrow IPC format.
///
/// Parameters
/// ----------
/// sink : str or file-like object
///     Path or file-like object to write to.
/// chr_names : list[str]
///     List of chromosome names.
/// compression : str, optional
///     Compression algorithm to use ('uncompressed', 'lz4', 'zstd'). Defaults to None (uncompressed).
/// custom_metadata : bytes, optional
///     Custom metadata to embed in the IPC file.
#[pyclass(name="BsxIpcWriter")]
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
        let file_like = PyFileLikeObject::with_requirements(sink, false, false, true, false)?;

        let writer = RsBsxIpcWriter::try_new(
            BufWriter::new(file_like),
            chr_names,
            compression.map(|x| x.into()),
            None,
        ).map_err(|e| PyIOError::new_err(e.to_string()))?;

        Ok(Self { writer: Some(writer) })
    }

    /// Create a writer using a FASTA index file (.fai) to get chromosome names.
    ///
    /// Parameters
    /// ----------
    /// sink : str or file-like object
    ///     Path or file-like object to write to.
    /// fai_path : str
    ///     Path to the FASTA index file.
    /// compression : str, optional
    ///     Compression algorithm ('uncompressed', 'lz4', 'zstd'). Defaults to None.
    /// custom_metadata : bytes, optional
    ///     Custom metadata.
    ///
    /// Returns
    /// -------
    /// BsxIpcWriter
    ///     A new writer instance.
    #[staticmethod]
    #[pyo3(signature = (sink, fai_path, compression=None))]
    pub fn from_sink_and_fai(
        sink: PyObject,
        fai_path: PathBuf,
        compression: Option<PyIpcCompression>,
    ) -> PyResult<Self> {
        let file_like = PyFileLikeObject::with_requirements(sink, false, false, true, false)?;

        let writer = RsBsxIpcWriter::try_from_sink_and_fai(
            BufWriter::new(file_like),
            fai_path,
            compression.map(|x| x.into()),
            None,
        ).map_err(|e| PyIOError::new_err(e.to_string()))?;

        Ok(Self { writer: Some(writer) })
    }

    /// Create a writer using a FASTA file to get chromosome names (will index if needed).
    ///
    /// Parameters
    /// ----------
    /// sink : str or file-like object
    ///     Path or file-like object to write to.
    /// fasta_path : str
    ///     Path to the FASTA file.
    /// compression : str, optional
    ///     Compression algorithm ('uncompressed', 'lz4', 'zstd'). Defaults to None.
    /// custom_metadata : bytes, optional
    ///     Custom metadata.
    ///
    /// Returns
    /// -------
    /// BsxIpcWriter
    ///     A new writer instance.
    #[staticmethod]
    #[pyo3(signature = (sink, fasta_path, compression=None))]
    pub fn from_sink_and_fasta(
        sink: PyObject,
        fasta_path: PathBuf,
        compression: Option<PyIpcCompression>,
    ) -> PyResult<Self> {
        let file_like = PyFileLikeObject::with_requirements(sink, false, false, true, false)?;

        let writer = RsBsxIpcWriter::try_from_sink_and_fasta(
            BufWriter::new(file_like),
            fasta_path,
            compression.map(|x| x.into()),
            None,
        ).map_err(|e| PyIOError::new_err(e.to_string()))?;

        Ok(Self { writer: Some(writer) })
    }

    /// Write an already encoded BSX batch (DataFrame).
    ///
    /// Parameters
    /// ----------
    /// batch : PyDataFrame
    ///     The encoded BSX batch (Polars DataFrame) to write.
    pub fn write_encoded_batch(&mut self, batch: PyDataFrame) -> PyResult<()> {
        let encoded_batch = unsafe { EncodedBsxBatch::new_unchecked(batch.0) };
        self.writer
            .as_mut()
            .ok_or_else(|| PyIOError::new_err("Writer is closed"))?
            .write_encoded_batch(encoded_batch)
            .map_err(|e| PyPolarsErr::Polars(e).into())
    }

    /// Encode and write a BSX batch (DataFrame).
    ///
    /// Parameters
    /// ----------
    /// batch : PyDataFrame
    ///     The BSX batch (Polars DataFrame) to encode and write.
    pub fn write_batch(&mut self, batch: PyDataFrame) -> PyResult<()> {
        let bsx_batch = unsafe { BsxBatch::new_unchecked(batch.0) };
        self.writer
            .as_mut()
            .ok_or_else(|| PyIOError::new_err("Writer is closed"))?
            .write_batch(bsx_batch)
            .map_err(|e| PyPolarsErr::Polars(e).into())
    }

    /// Finalize the IPC file and close the writer.
    pub fn close(&mut self) -> PyResult<()> {
        if let Some(mut writer) = self.writer.take() {
            writer.close().map_err(|e| PyPolarsErr::Polars(e))?;
        }
        Ok(())
    }

    pub fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    pub fn __exit__(&mut self, _exc_type: PyObject, _exc_value: PyObject, _traceback: PyObject) -> PyResult<()> {
        self.close()
    }
}
