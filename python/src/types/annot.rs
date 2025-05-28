use std::path::PathBuf;

use bsxplorer2::data_structs::annotation::{
    GffEntry,
    GffEntryAttributes,
    HcAnnotStore,
    RawGffEntry,
};
use bsxplorer2::data_structs::coords::Contig;
use pyo3::exceptions::{
    PyFileNotFoundError,
    PyIOError,
    PyValueError,
};
use pyo3::prelude::*;

use crate::types::coords::PyContig;

#[pyclass(name = "AnnotStore")]
pub struct PyAnnotStore {
    inner: HcAnnotStore,
}

impl From<HcAnnotStore> for PyAnnotStore {
    fn from(inner: HcAnnotStore) -> Self {
        Self { inner }
    }
}

impl From<PyAnnotStore> for HcAnnotStore {
    fn from(py_annot_store: PyAnnotStore) -> Self {
        py_annot_store.inner
    }
}

#[pymethods]
impl PyAnnotStore {
    #[new]
    pub fn new() -> Self {
        Self {
            inner: HcAnnotStore::new(),
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

        let annot = HcAnnotStore::from_iter(entries);
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
        let mut annot_store = HcAnnotStore::new();

        for record in reader.records() {
            let record = record.map_err(|e| {
                PyValueError::new_err(format!("Error reading BED record: {}", e))
            })?;

            let entry = GffEntry::from(record);
            annot_store.insert(entry).map_err(|e| {
                PyValueError::new_err(format!("Failed to insert entry: {}", e))
            })?;
        }

        Ok(PyAnnotStore { inner: annot_store })
    }

    #[staticmethod]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            inner: HcAnnotStore::with_capacity(capacity),
        }
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    pub fn insert(
        &mut self,
        entry: PyGffEntry,
    ) -> PyResult<()> {
        self.inner.insert(entry.into()).map_err(|e| {
            PyValueError::new_err(format!("Failed to insert entry: {}", e))
        })?;
        Ok(())
    }

    pub fn get_entry(
        &self,
        id: String,
    ) -> Option<PyGffEntry> {
        self.inner
            .get_entry(&id.into())
            .map(|entry| PyGffEntry::from(entry.clone()))
    }

    pub fn get_entries_regex(
        &self,
        pattern: String,
    ) -> Vec<PyGffEntry> {
        self.inner
            .get_entries_regex(&pattern).unwrap()
            .into_iter()
            .map(|entry| PyGffEntry::from(entry.clone()))
            .collect()
    }

    pub fn get_children(
        &self,
        id: String,
    ) -> Option<Vec<String>> {
        self.inner
            .get_children(&id.into())
            .map(|entry| entry.into_iter().map(|child| child.to_string()).collect())
    }

    pub fn get_parent(
        &self,
        id: String,
    ) -> Option<String> {
        self.inner
            .get_parent(&id.into())
            .map(|parent| parent.to_string())
    }

    pub fn genomic_query(
        &self,
        contig: &PyContig,
    ) -> Option<Vec<String>> {
        let rust_contig: Contig = contig.clone().into();
        self.inner
            .genomic_query(&rust_contig)
            .map(|ids| ids.into_iter().map(|id| id.to_string()).collect())
    }

    pub fn get_feature_types(&self) -> Vec<String> {
        self.inner
            .get_feature_types()
            .into_iter()
            .map(|s| s.to_string())
            .collect()
    }

    pub fn add_flank(
        &mut self,
        selector: PyObject,
        flank: i32,
        prefix: String,
    ) -> PyResult<()> {
        Python::with_gil(|py| {
            let selector_fn = |entry: &GffEntry| -> bool {
                // Create a PyGffEntry from the Rust GffEntry
                let py_entry = PyGffEntry::from(entry.clone());

                // Call the Python selector function with the PyGffEntry
                match selector.call1(py, (py_entry,)) {
                    Ok(result) => {
                        // Convert the Python result to a boolean
                        match result.extract::<bool>(py) {
                            Ok(bool_result) => bool_result,
                            Err(_) => {
                                // If extraction fails, default to false
                                false
                            },
                        }
                    },
                    Err(_) => {
                        // If the call fails, default to false
                        false
                    },
                }
            };

            self.inner.add_flank(selector_fn, flank, &prefix);
            Ok(())
        })
    }

    pub fn sort_self(&mut self) -> PyResult<()> {
        self.inner
            .sort_self()
            .map_err(|e| PyValueError::new_err(format!("Failed to sort: {}", e)))
    }

    pub fn iter(&self) -> PyAnnotStoreIterator {
        PyAnnotStoreIterator {
            entries: self.inner.iter().cloned().collect(),
            index:   0,
        }
    }

    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    pub fn __repr__(&self) -> String {
        format!("AnnotStore(entries={})", self.inner.len())
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<PyAnnotStoreIterator>> {
        let iter = PyAnnotStoreIterator {
            entries: slf.inner.iter().cloned().collect(),
            index:   0,
        };
        Py::new(slf.py(), iter)
    }
}

#[pyclass]
pub struct PyAnnotStoreIterator {
    entries: Vec<GffEntry>,
    index:   usize,
}

#[pymethods]
impl PyAnnotStoreIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyGffEntry> {
        if slf.index < slf.entries.len() {
            let entry = slf.entries[slf.index].clone();
            slf.index += 1;
            Some(PyGffEntry::from(entry))
        }
        else {
            None
        }
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
        let rust_contig: Contig = contig.clone().into();

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
            rust_contig,
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

    #[getter]
    fn get_attributes(&self) -> PyGffEntryAttributes {
        PyGffEntryAttributes::from(self.inner.attributes.clone())
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

#[pyclass(name = "GffEntryAttributes")]
#[derive(Debug, Clone)]
pub struct PyGffEntryAttributes {
    inner: GffEntryAttributes,
}

impl From<GffEntryAttributes> for PyGffEntryAttributes {
    fn from(inner: GffEntryAttributes) -> Self {
        Self { inner }
    }
}

impl From<PyGffEntryAttributes> for GffEntryAttributes {
    fn from(py_attrs: PyGffEntryAttributes) -> Self {
        py_attrs.inner
    }
}

#[pymethods]
impl PyGffEntryAttributes {
    #[new]
    fn new() -> Self {
        Self {
            inner: GffEntryAttributes::default(),
        }
    }

    #[getter]
    fn get_id(&self) -> Option<String> {
        self.inner.id().map(|s| s.to_string())
    }

    #[getter]
    fn get_name(&self) -> Option<Vec<String>> {
        self.inner
            .name()
            .map(|names| names.iter().map(|s| s.to_string()).collect())
    }

    #[getter]
    fn get_alias(&self) -> Option<Vec<String>> {
        self.inner
            .alias()
            .map(|aliases| aliases.iter().map(|s| s.to_string()).collect())
    }

    #[getter]
    fn get_parent(&self) -> Option<Vec<String>> {
        self.inner
            .parent()
            .map(|parents| parents.iter().map(|s| s.to_string()).collect())
    }

    fn __repr__(&self) -> String {
        format!("GffEntryAttributes({})", self.inner)
    }
}
