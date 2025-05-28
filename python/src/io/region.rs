use std::fs::File;

use bsxplorer2::data_structs::batch::BsxBatch;
use bsxplorer2::io::bsx::{
    BatchIndex,
    BsxFileReader,
    RegionReader,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3_polars::error::PyPolarsErr;
use pyo3_polars::PyDataFrame;

use super::bsx::PyBsxFileReader;
use crate::types::batch::PyBsxBatch;
use crate::types::coords::PyContig;
use crate::types::index::PyBatchIndex;
use crate::types::utils::{
    PyContext,
    PyStrand,
};
use crate::utils::FileOrFileLike;

#[derive(Clone)]
#[pyclass(name = "FilterOperation")]
pub enum PyFilterOperation {
    PosLt { value: u32 },
    PosGt { value: u32 },
    CoverageGt { value: u32 },
    Strand { value: PyStrand },
    Context { value: PyContext },
}

fn apply_filters(
    batch: BsxBatch,
    filters: &[PyFilterOperation],
) -> BsxBatch {
    if filters.is_empty() {
        return batch;
    }

    let mut lazy_batch = batch.lazy();

    for filter in filters {
        lazy_batch = match filter {
            PyFilterOperation::PosLt { value } => lazy_batch.filter_pos_lt(*value),
            PyFilterOperation::PosGt { value } => lazy_batch.filter_pos_gt(*value),
            PyFilterOperation::CoverageGt { value } => {
                lazy_batch.filter_coverage_gt(*value)
            },
            PyFilterOperation::Strand { value } => {
                lazy_batch.filter_strand(value.clone().into())
            },
            PyFilterOperation::Context { value } => {
                lazy_batch.filter_context(value.clone().into())
            },
        };
    }

    // This should not fail since we're only applying filters
    lazy_batch.collect().unwrap()
}

#[pyclass(unsendable, name = "RegionReader")]
pub struct PyRegionReader {
    inner:   Option<RegionReader>,
    filters: Vec<PyFilterOperation>,
}

#[pyclass(unsendable, name = "RegionReaderIterator")]
pub struct PyRegionReaderIterator {
    inner:         RegionReader,
    contigs:       Vec<PyContig>,
    current_index: usize,
    filters:       Vec<PyFilterOperation>,
}

#[pymethods]
impl PyRegionReader {
    #[new]
    fn new(file: FileOrFileLike) -> PyResult<Self> {
        let mut reader = match file {
            FileOrFileLike::File(path) => BsxFileReader::try_new(File::open(path)?)?,
            FileOrFileLike::ROnlyFileLike(handle) => BsxFileReader::try_new(handle)?,
            FileOrFileLike::RWFileLike(handle) => BsxFileReader::try_new(handle)?,
        };
        let index: BatchIndex =
            BatchIndex::from_reader(&mut reader).map_err(|e| PyPolarsErr::Polars(e))?;

        Ok(Self {
            inner:   Some(RegionReader::new(reader, index, None)),
            filters: Vec::new(),
        })
    }

    #[classmethod]
    fn from_reader(
        _cls: &Bound<'_, pyo3::types::PyType>,
        reader: PyBsxFileReader,
    ) -> PyResult<Self> {
        let inner = RegionReader::from_reader(reader.into_inner())?;

        Ok(Self {
            inner:   Some(inner),
            filters: Vec::new(),
        })
    }

    fn query(
        &mut self,
        contig: PyContig,
    ) -> PyResult<Option<PyDataFrame>> {
        let reader = self.inner.as_mut().ok_or_else(|| {
            PyValueError::new_err("RegionReader has been moved to an iterator")
        })?;

        let result = reader.query(contig.into(), None);

        match result {
            Ok(Some(batch)) => {
                let final_batch = apply_filters(batch, &self.filters);
                Ok(Some(PyDataFrame(final_batch.into_inner())))
            },
            Ok(None) => Ok(None),
            Err(e) => Err(PyErr::from(e)),
        }
    }

    fn reset(&mut self) {
        if let Some(ref mut reader) = self.inner {
            reader.reset();
        }
    }

    fn chr_order(&self) -> PyResult<Vec<String>> {
        let reader = self.inner.as_ref().ok_or_else(|| {
            PyValueError::new_err("RegionReader has been moved to an iterator")
        })?;

        Ok(reader
            .index()
            .get_chr_order()
            .iter()
            .cloned()
            .map(|chr| chr.to_string())
            .collect())
    }

    fn add_filter(
        &mut self,
        filter: PyFilterOperation,
    ) {
        self.filters.push(filter);
    }

    fn clear_filters(&mut self) {
        self.filters.clear();
    }

    fn filter_pos_lt(
        &mut self,
        value: u32,
    ) {
        self.add_filter(PyFilterOperation::PosLt { value });
    }

    fn filter_pos_gt(
        &mut self,
        value: u32,
    ) {
        self.add_filter(PyFilterOperation::PosGt { value });
    }

    fn filter_coverage_gt(
        &mut self,
        value: u32,
    ) {
        self.add_filter(PyFilterOperation::CoverageGt { value });
    }

    fn filter_strand(
        &mut self,
        value: PyStrand,
    ) {
        self.add_filter(PyFilterOperation::Strand { value });
    }

    fn filter_context(
        &mut self,
        value: PyContext,
    ) {
        self.add_filter(PyFilterOperation::Context { value });
    }

    fn index(&self) -> PyResult<PyBatchIndex> {
        let reader = self.inner.as_ref().ok_or_else(|| {
            PyValueError::new_err("RegionReader has been moved to an iterator")
        })?;

        Ok(reader.index().clone().into())
    }

    fn iter_contigs(
        &mut self,
        contigs: Vec<PyContig>,
    ) -> PyResult<PyRegionReaderIterator> {
        let reader = self.inner.take().ok_or_else(|| {
            PyValueError::new_err("RegionReader has already been moved to an iterator")
        })?;

        Ok(PyRegionReaderIterator {
            inner: reader,
            contigs,
            current_index: 0,
            filters: self.filters.clone(),
        })
    }
}

#[pymethods]
impl PyRegionReaderIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> PyResult<Option<PyBsxBatch>> {
        if self.current_index < self.contigs.len() {
            let contig = self.contigs[self.current_index].clone();
            self.current_index += 1;

            let result = self.inner.query(contig.into(), None);
            match result {
                Ok(Some(batch)) => {
                    let final_batch = apply_filters(batch, &self.filters);
                    Ok(Some(final_batch.into()))
                },
                Ok(None) => self.__next__(), // Skip to next contig if no data
                Err(e) => Err(PyErr::from(e)),
            }
        }
        else {
            Ok(None) // StopIteration
        }
    }
}
