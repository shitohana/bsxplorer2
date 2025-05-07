use std::path::PathBuf;

use bsxplorer2::data_structs::annotation::{AnnotStore, GffEntry};
use pyo3::exceptions::{PyFileNotFoundError, PyValueError};
use pyo3::prelude::*;

use super::index::PyBatchIndex;
use crate::types::coords::PyContig;

#[pyclass(name = "AnnotStore")]
pub struct PyAnnotStore {
    inner: AnnotStore,
}

impl From<AnnotStore> for PyAnnotStore {
    fn from(inner: AnnotStore) -> Self { Self { inner } }
}

impl From<PyAnnotStore> for AnnotStore {
    fn from(py_annot_store: PyAnnotStore) -> Self { py_annot_store.inner }
}

#[pymethods]
impl PyAnnotStore {
    #[new]
    pub fn new() -> Self {
        Self {
            inner: AnnotStore::new(),
        }
    }

    #[staticmethod]
    pub fn from_gff(path: PathBuf) -> PyResult<Self> { todo!() }

    #[staticmethod]
    pub fn from_bed(path: PathBuf) -> PyResult<Self> {
        use std::fs::File;

        let file = File::open(&path).map_err(|e| {
            PyFileNotFoundError::new_err(format!(
                "Failed to open BED file: {}",
                e
            ))
        })?;

        use bio::io::bed;
        let mut reader = bed::Reader::new(file);
        let mut annot_store = AnnotStore::new();

        for record in reader.records() {
            let record = record.map_err(|e| {
                PyValueError::new_err(format!(
                    "Error reading BED record: {}",
                    e
                ))
            })?;

            let entry = GffEntry::from(record);
            annot_store.insert(entry);
        }

        Ok(PyAnnotStore { inner: annot_store })
    }

    pub fn len(&self) -> usize { self.inner.len() }

    pub fn is_empty(&self) -> bool { self.inner.len() == 0 }

    pub fn add_upstream(
        &mut self,
        length: u32,
    ) -> PyResult<()> {
        self.inner
            .add_upstream(|_| true, length);
        Ok(())
    }

    pub fn add_downstream(
        &mut self,
        length: u32,
    ) -> PyResult<()> {
        self.inner
            .add_downstream(|_| true, length);
        Ok(())
    }

    pub fn add_flanks(
        &mut self,
        length: u32,
    ) -> PyResult<()> {
        self.inner.add_flanks(|_| true, length);
        Ok(())
    }

    pub fn iter_sorted(
        &self,
        index: &PyBatchIndex,
    ) -> PyResult<Vec<PyContig>> {
        let batch_index = &index.inner();
        let results: Vec<PyContig> = self
            .inner
            .iter_sorted(batch_index)
            .map(|(_, entry)| {
                let contig = entry.contig.clone();
                PyContig::new(
                    contig.seqname().to_string(),
                    contig.start(),
                    contig.end(),
                    contig.strand().to_string().as_str(),
                )
                .unwrap()
            })
            .collect();

        Ok(results)
    }

    pub fn get_entry_ids(&self) -> Vec<String> {
        self.inner
            .id_map()
            .keys()
            .map(|id| id.to_string())
            .collect()
    }

    pub fn __len__(&self) -> usize { self.inner.len() }

    pub fn __repr__(&self) -> String {
        format!("AnnotStore(entries={})", self.inner.len())
    }
}
