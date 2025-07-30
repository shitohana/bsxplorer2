use std::collections::{
    HashMap,
    HashSet,
};
use std::path::PathBuf;

use bsxplorer2::data_structs::annotation::{
    EntryId,
    GffEntry,
    GffEntryAttributes,
    HcAnnotStore,
};
use bsxplorer2::data_structs::coords::Contig;
use pyo3::exceptions::{
    PyFileNotFoundError,
    PyRuntimeError,
    PyValueError,
};
use pyo3::prelude::*;
use slotmap::{
    Key,
    KeyData,
};

use crate::types::coords::PyContig;

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

    // Add a getter for the 'other' hashmap
    #[getter]
    fn get_other(&self) -> HashMap<String, String> {
        HashMap::from_iter(self.inner.other.clone().into_iter())
    }

    fn __repr__(&self) -> String {
        format!("GffEntryAttributes({})", self.inner)
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
        id: Option<String>, /* This is the GFF attribute ID (BsxSmallStr), not
                             * EntryId (u64) */
    ) -> PyResult<Self> {
        let rust_contig: Contig = contig.clone().into();

        let mut attributes = GffEntryAttributes::default();
        if let Some(id_str) = id {
            attributes = attributes.with_id(id_str);
        }

        let arc_source = source.map(|s| s.into());
        let arc_feature_type = feature_type.map(|s| s.into());

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
    // This getter returns the GFF attribute ID string, not the internal EntryId
    // (u64)
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
            "GffEntry(id='{}', feature='{}', {}:{}-{})", /* Use single quotes for id
                                                          * string */
            self.inner.id,
            self.inner.feature_type,
            self.inner.contig.seqname(),
            self.inner.contig.start(),
            self.inner.contig.end()
        )
    }
}

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
        let file = std::fs::File::open(&path).map_err(|e| {
            PyFileNotFoundError::new_err(format!("Failed to open GFF file: {}", e))
        })?;
        let annot_store = HcAnnotStore::from_gff(file).map_err(|e| {
            PyValueError::new_err(format!("Failed to parse GFF: {}", e))
        })?;
        Ok(Self { inner: annot_store })
    }

    #[staticmethod]
    pub fn from_bed(path: PathBuf) -> PyResult<Self> {
        let file = std::fs::File::open(&path).map_err(|e| {
            PyFileNotFoundError::new_err(format!("Failed to open BED file: {}", e))
        })?;
        let annot_store = HcAnnotStore::from_bed(file).map_err(|e| {
            PyValueError::new_err(format!("Failed to parse BED: {}", e))
        })?;

        Ok(PyAnnotStore { inner: annot_store })
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    // Returns the u64 EntryId of the inserted entry
    pub fn insert(
        &mut self,
        entry: PyGffEntry,
    ) -> PyResult<u64> {
        let id = self.inner.insert(entry.into());
        Ok(id.data().as_ffi())
    }

    // Retrieves an entry by its internal EntryId (u64)
    pub fn get_entry(
        &self,
        id_u64: u64,
    ) -> Option<PyGffEntry> {
        let entry_id = EntryId::from(KeyData::from_ffi(id_u64));
        self.inner
            .get(entry_id)
            .map(|entry| PyGffEntry::from(entry.clone()))
    }

    // Retrieves entries whose GFF attribute ID matches a regex pattern (uses
    // GffEntry.id: BsxSmallStr)
    pub fn get_entries_regex(
        &self,
        pattern: String,
    ) -> PyResult<Vec<PyGffEntry>> {
        self.inner
            .get_entries_regex(&pattern)
            .map(|entries| {
                entries
                    .into_iter()
                    .map(|entry| PyGffEntry::from(entry.clone()))
                    .collect()
            })
            .map_err(|e| PyValueError::new_err(format!("Invalid regex pattern: {}", e)))
    }

    // Initializes the internal tree structure
    pub fn init_tree(&mut self) -> PyResult<()> {
        self.inner.init_tree().map_err(|e| {
            PyRuntimeError::new_err(format!("Failed to initialize tree: {}", e))
        })
    }

    // Initializes the internal interval map
    pub fn init_imap(&mut self) -> PyResult<()> {
        self.inner.init_imap(); // init_imap does not return Result
        Ok(())
    }

    // Retrieves the internal EntryIds (u64) of direct children of a given EntryId
    // (u64)
    pub fn get_children(
        &self,
        id_u64: u64,
    ) -> PyResult<Option<Vec<u64>>> {
        let entry_id = EntryId::from(KeyData::from_ffi(id_u64));
        self.inner
            .get_children(entry_id)
            .map(|opt_ids| {
                opt_ids.map(|ids| {
                    ids.into_iter()
                        .map(|child_id| child_id.data().as_ffi())
                        .collect()
                })
            })
            .map_err(|e| {
                PyRuntimeError::new_err(format!("Failed to get children: {}", e))
            })
    }

    // Retrieves the internal EntryId (u64) of the direct parent of a given EntryId
    // (u64)
    pub fn get_parent(
        &self,
        id_u64: u64,
    ) -> PyResult<Option<u64>> {
        let entry_id = EntryId::from(KeyData::from_ffi(id_u64));
        self.inner
            .get_parent(&entry_id) // get_parent takes &EntryId
            .map(|opt_id| opt_id.map(|id| id.data().as_ffi()))
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to get parent: {}", e)))
    }

    // Performs a genomic query and returns a list of matching internal EntryIds
    // (u64)
    pub fn genomic_query(
        &self,
        contig: &PyContig,
    ) -> PyResult<Vec<u64>> {
        let rust_contig: Contig = contig.clone().into();
        self.inner
            .genomic_query(&rust_contig)
            .map(|ids| ids.into_iter().map(|id| id.data().as_ffi()).collect())
            .map_err(|e| {
                PyRuntimeError::new_err(format!("Genomic query failed: {}", e))
            })
    }

    // Gets feature types present in the store mapped to internal EntryIds (u64)
    pub fn get_feature_types(&self) -> HashMap<String, Vec<u64>> {
        self.inner
            .get_feature_types()
            .into_iter()
            .map(|(feat_type, ids)| {
                let ids_u64 = ids.into_iter().map(|id| id.data().as_ffi()).collect();
                (feat_type.to_string(), ids_u64)
            })
            .collect()
    }

    // Adds flanking regions to entries selected by their internal EntryId (u64)
    pub fn add_flanks(
        &mut self,
        ids_u64: Vec<u64>,
        flank: i32,
        prefix: String,
    ) -> PyResult<()> {
        // Convert u64 list to HashSet of EntryId
        let id_set: HashSet<EntryId> = ids_u64
            .into_iter()
            .map(|id_u64| EntryId::from(KeyData::from_ffi(id_u64)))
            .collect();

        // Iterate over existing entries, check if their ID is in the set, and add
        // flanks This replicates the logic from the Rust add_flank but selects
        // using the provided id_set.

        let selected_entries: Vec<(EntryId, GffEntry)> = self
            .inner
            .iter() // yields (EntryId, &GffEntry)
            .filter(|(id, _entry)| id_set.contains(id)) // filter by internal EntryId
            .map(|(id, entry)| (id.clone(), entry.clone())) // clone EntryId and GffEntry
            .collect();

        for (id, parent_entry) in selected_entries {
            let (start, end) = if flank > 0 {
                // Flank downstream (after end)
                (
                    parent_entry.contig.end_gpos(),
                    parent_entry.contig.end_gpos().shift(flank as isize),
                )
            }
            else {
                // Flank upstream (before start)
                (
                    parent_entry.contig.start_gpos().shift(flank as isize),
                    parent_entry.contig.start_gpos(),
                )
            };

            // Ensure start <= end for the range
            let (start, end) = if start <= end {
                (start, end)
            }
            else {
                (end, start)
            };

            let mut feature_type = prefix.to_string();
            feature_type.push_str(parent_entry.feature_type.as_str());

            // Create a new unique ID for the flank entry (this is the GFF attribute ID
            // string) We don't need to predict the EntryId here, it's
            // assigned on insert.
            let flank_attribute_id_str =
                format!("{}_flank_{}", parent_entry.id(), flank);

            let flank_entry = GffEntry::new(
                (start..end).into(),
                None,
                Some(feature_type.into()),
                None,
                None,
                Some(
                    GffEntryAttributes::default()
                        .with_id(flank_attribute_id_str) // Set the unique attribute ID
                        .with_parent(vec![parent_entry.id().clone()]), /* Set GFF parent attribute to original entry's GFF attribute ID */
                ),
            );

            // Use insert which handles adding to all internal structures and assigns
            // the new EntryId
            self.inner.insert(flank_entry);
        }

        Ok(())
    }

    // Remove the sort_self method as it doesn't exist in HcAnnotStore

    // The iterator will yield (u64, PyGffEntry) tuples
    pub fn iter(&self) -> PyResult<PyAnnotStoreIterator> {
        let iter = PyAnnotStoreIterator {
            entries: self
                .inner
                .iter()
                .map(|(id, entry)| (id.data().as_ffi(), entry.clone()))
                .collect(),
            index:   0,
        };
        Ok(iter)
    }

    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    pub fn __repr__(&self) -> String {
        format!("AnnotStore(entries={})", self.inner.len())
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<PyAnnotStoreIterator>> {
        let iter = PyAnnotStoreIterator {
            entries: slf
                .inner
                .iter()
                .map(|(id, entry)| (id.data().as_ffi(), entry.clone()))
                .collect(),
            index:   0,
        };
        Py::new(slf.py(), iter)
    }
}

#[pyclass]
// Iterator now yields (u64 EntryId, GffEntry)
pub struct PyAnnotStoreIterator {
    entries: Vec<(u64, GffEntry)>,
    index:   usize,
}

#[pymethods]
impl PyAnnotStoreIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    // Yields a tuple: (u64 EntryId, PyGffEntry)
    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<(u64, PyGffEntry)> {
        if slf.index < slf.entries.len() {
            let (id_u64, entry) = slf.entries[slf.index].clone();
            slf.index += 1;
            Some((id_u64, PyGffEntry::from(entry)))
        }
        else {
            None
        }
    }
}
