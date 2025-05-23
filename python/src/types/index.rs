use std::fs::File;
use std::io::{BufReader, BufWriter};

use bsxplorer2::data_structs::coords::Contig;
use bsxplorer2::data_structs::Strand;
use bsxplorer2::io::bsx::BatchIndex;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::coords::PyContig;

/// Python wrapper for BatchIndex
#[pyclass(name = "BatchIndex")]
#[derive(Debug, Clone)]
pub struct PyBatchIndex {
    inner: BatchIndex,
}

impl From<BatchIndex> for PyBatchIndex {
    fn from(index: BatchIndex) -> Self {
        PyBatchIndex { inner: index }
    }
}

impl From<PyBatchIndex> for BatchIndex {
    fn from(py_index: PyBatchIndex) -> Self {
        py_index.inner
    }
}

impl PyBatchIndex {
    pub fn inner(&self) -> &BatchIndex {
        &self.inner
    }
}

#[pymethods]
impl PyBatchIndex {
    #[new]
    fn new() -> Self {
        PyBatchIndex {
            inner: BatchIndex::new(),
        }
    }

    /// Insert a contig and its corresponding batch index
    #[pyo3(
        text_signature = "($self, seqname: str, start: int, end: int, batch_idx: int)"
    )]
    fn insert(
        &mut self,
        seqname: String,
        start: u32,
        end: u32,
        batch_idx: usize,
    ) -> PyResult<()> {
        let contig = Contig::new(seqname.into(), start, end, Strand::None);
        self.inner.insert(contig, batch_idx);
        Ok(())
    }

    /// Find batch indices that overlap with a given contig
    #[pyo3(text_signature = "($self, seqname: str, start: int, end: int)")]
    fn find(
        &self,
        seqname: String,
        start: u32,
        end: u32,
    ) -> PyResult<Option<Vec<usize>>> {
        let contig = Contig::new(seqname.into(), start, end, Strand::None);
        Ok(self.inner.find(&contig))
    }

    /// Get chromosome order
    fn get_chr_order(&self) -> PyResult<Vec<String>> {
        Ok(self
            .inner
            .get_chr_order()
            .iter()
            .cloned()
            .map(|v| v.to_string())
            .collect())
    }

    /// Get index of a chromosome in the order
    fn get_chr_index(
        &self,
        chr: String,
    ) -> PyResult<Option<usize>> {
        Ok(self.inner.get_chr_index(&chr.into()))
    }

    /// Save index to a file
    #[pyo3(text_signature = "($self, filename)")]
    fn save(
        &self,
        filename: String,
    ) -> PyResult<()> {
        let file =
            File::create(filename).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut writer = BufWriter::new(file);
        self.inner
            .clone()
            .to_file(&mut writer)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(())
    }

    /// Load index from a file
    #[staticmethod]
    #[pyo3(text_signature = "(filename)")]
    fn load(filename: String) -> PyResult<Self> {
        let file =
            File::open(filename).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut reader = BufReader::new(file);
        let inner = BatchIndex::from_file(&mut reader)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(PyBatchIndex { inner })
    }

    pub fn sort(
        &self,
        contigs: Vec<PyContig>,
    ) -> Vec<PyContig> {
        self.inner
            .sort(contigs.into_iter().map(|contig| {
                Contig::new(
                    contig.seqname.into(),
                    contig.start,
                    contig.end,
                    contig.strand.into(),
                )
            }))
            .map(|contig| contig.into())
            .collect()
    }
}
