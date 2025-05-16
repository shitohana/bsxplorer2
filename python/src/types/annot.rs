use std::collections::HashMap;
use std::path::PathBuf;

use bsxplorer2::data_structs::annotation::{AnnotStore,
                                           GffEntry,
                                           GffEntryAttributes,
                                           RawGffEntry};
use bsxplorer2::data_structs::coords::Contig;
use pyo3::exceptions::{PyFileNotFoundError, PyIOError, PyValueError};
use pyo3::prelude::*;

use super::index::PyBatchIndex;
use crate::types::coords::PyContig;

#[pyclass(name = "AnnotStore")]
pub struct PyAnnotStore {
    inner: AnnotStore,
}

impl From<AnnotStore> for PyAnnotStore {
    fn from(inner: AnnotStore) -> Self {
        Self { inner }
    }
}

impl From<PyAnnotStore> for AnnotStore {
    fn from(py_annot_store: PyAnnotStore) -> Self {
        py_annot_store.inner
    }
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
    pub fn from_gff(path: PathBuf) -> PyResult<Self> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .flexible(true)
            .ascii()
            .from_path(path)
            .map_err(|e| {
                PyIOError::new_err(format!("Failed to read GFF file: {}", e))
            })?;

        let mut entries = Vec::new();
        for result in reader.deserialize() {
            let entry: RawGffEntry =
                result.map_err(|e| PyValueError::new_err(format!("{}", e)))?;
            entries.push(GffEntry::try_from(entry)?);
        }

        let annot = AnnotStore::from_iter(entries);
        Ok(Self { inner: annot })
    }

    #[staticmethod]
    pub fn from_bed(path: PathBuf) -> PyResult<Self> {
        use std::fs::File;

        let file = File::open(&path).map_err(|e| {
            PyFileNotFoundError::new_err(format!("Failed to open BED file: {}", e))
        })?;

        use bio::io::bed;
        let mut reader = bed::Reader::new(file);
        let mut annot_store = AnnotStore::new();

        for record in reader.records() {
            let record = record.map_err(|e| {
                PyValueError::new_err(format!("Error reading BED record: {}", e))
            })?;

            let entry = GffEntry::from(record);
            annot_store.insert(entry);
        }

        Ok(PyAnnotStore { inner: annot_store })
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.len() == 0
    }

    pub fn insert(
        &mut self,
        entry: PyGffEntry,
    ) -> PyResult<()> {
        self.inner.insert(entry.into());
        Ok(())
    }

    pub fn remove(
        &mut self,
        id: String,
    ) -> Option<PyGffEntry> {
        self.inner
            .remove(&id.into())
            .map(|entry| PyGffEntry::from(entry))
    }

    pub fn get_children(
        &self,
        id: String,
    ) -> Option<Vec<String>> {
        self.inner
            .get_children(&id.into())
            .map(|entry| entry.into_iter().map(|child| child.to_string()).collect())
    }

    pub fn get_parents(
        &self,
        id: String,
    ) -> Option<Vec<String>> {
        self.inner
            .get_parents(&id.into())
            .map(|entry| entry.into_iter().map(|parent| parent.to_string()).collect())
    }

    pub fn id_map(&self) -> HashMap<String, PyGffEntry> {
        self.inner
            .id_map()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.clone().into()))
            .collect()
    }

    pub fn parent_map(&self) -> HashMap<String, Vec<String>> {
        self.inner
            .parent_map()
            .iter_all()
            .map(|(k, v)| {
                (
                    k.to_string(),
                    v.into_iter().map(|parent| parent.to_string()).collect(),
                )
            })
            .collect()
    }

    pub fn children_map(&self) -> HashMap<String, Vec<String>> {
        self.inner
            .children_map()
            .iter_all()
            .map(|(k, v)| {
                (
                    k.to_string(),
                    v.into_iter().map(|child| child.to_string()).collect(),
                )
            })
            .collect()
    }

    pub fn add_upstream(
        &mut self,
        length: u32,
    ) -> PyResult<()> {
        self.inner.add_upstream(|_| true, length);
        Ok(())
    }

    pub fn add_downstream(
        &mut self,
        length: u32,
    ) -> PyResult<()> {
        self.inner.add_downstream(|_| true, length);
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

    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    pub fn __repr__(&self) -> String {
        format!("AnnotStore(entries={})", self.inner.len())
    }
}
#[pyclass(name = "GffEntry")]
#[derive(Debug, Clone)]
pub struct PyGffEntry {
    inner: GffEntry,
}

impl From<GffEntry> for PyGffEntry {
    fn from(inner: GffEntry) -> Self {
        Self { inner }
    }
}

impl From<PyGffEntry> for GffEntry {
    fn from(py_gff_entry: PyGffEntry) -> Self {
        py_gff_entry.inner
    }
}

#[pymethods]
impl PyGffEntry {
    #[new]
    #[pyo3(signature = (contig, source=None, feature_type=None, score=None, phase=None, id=None))]
    fn from_contig(
        contig: &PyContig,
        source: Option<String>,
        feature_type: Option<String>,
        score: Option<f64>,
        phase: Option<u8>,
        id: Option<String>,
    ) -> PyResult<Self> {
        // Convert PyContig to Contig<ArcStr, u32>
        let rust_contig: Contig<String, u32> = contig.into();
        let arc_contig = Contig::new(
            rust_contig.seqname().into(),
            rust_contig.start(),
            rust_contig.end(),
            rust_contig.strand(),
        );

        // Create attributes with optional ID
        let mut attributes = GffEntryAttributes::default();
        if let Some(id_str) = id {
            attributes = attributes.with_id(Some(id_str));
        }

        // Convert Optional<String> to Optional<ArcStr>
        let arc_source = source.map(|s| s.into());
        let arc_feature_type = feature_type.map(|s| s.into());

        // Create GffEntry
        let entry = GffEntry::new(
            arc_contig,
            arc_source,
            arc_feature_type,
            score,
            phase,
            Some(attributes),
        );

        Ok(Self { inner: entry })
    }

    #[getter]
    fn get_id(&self) -> String {
        self.inner.id.to_string()
    }

    #[getter]
    fn get_contig(&self) -> PyResult<PyContig> {
        let contig = &self.inner.contig;
        PyContig::new(
            contig.seqname().to_string(),
            contig.start(),
            contig.end(),
            contig.strand().to_string().as_str(),
        )
    }

    #[getter]
    fn get_source(&self) -> String {
        self.inner.source.to_string()
    }

    #[getter]
    fn get_feature_type(&self) -> String {
        self.inner.feature_type.to_string()
    }

    #[getter]
    fn get_score(&self) -> Option<f64> {
        self.inner.score
    }

    #[getter]
    fn get_phase(&self) -> Option<u8> {
        self.inner.phase
    }

    fn __repr__(&self) -> String {
        format!(
            "GffEntry(id={}, feature={}, {}:{}-{})",
            self.inner.id,
            self.inner.feature_type,
            self.inner.contig.seqname(),
            self.inner.contig.start(),
            self.inner.contig.end()
        )
    }
}
