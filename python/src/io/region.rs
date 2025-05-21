use bsxplorer2::data_structs::batch::BsxBatch;
use bsxplorer2::io::bsx::{BatchIndex, BsxFileReader, RegionReader};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use pyo3_polars::PyDataFrame;

use crate::types::batch::PyBsxBatch;
use crate::types::coords::PyContig;
use crate::types::index::PyBatchIndex;
use crate::utils::{FileOrFileLike, ReadHandle};

#[pyclass(unsendable, name = "RegionReader")]
pub struct PyRegionReader {
    inner: RegionReader<Box<dyn ReadHandle>>,
}

#[pymethods]
impl PyRegionReader {
    #[new]
    fn new(file: FileOrFileLike) -> PyResult<Self> {
        let file = file.get_reader()?;
        let mut bsx_reader: BsxFileReader<Box<_>> = BsxFileReader::new(file);
        let index: BatchIndex = bsx_reader.index().map_err(PyErr::from)?.clone();

        Ok(Self {
            inner: RegionReader::new(bsx_reader, index, None),
        })
    }

    fn query(
        &mut self,
        contig: PyContig,
        postprocess_fn: Py<PyAny>,
    ) -> PyResult<Option<PyDataFrame>> {
        Python::with_gil(|py| {
            if !postprocess_fn.clone_ref(py).into_bound(py).is_callable() {
                return Err(PyValueError::new_err(
                    "Preprocess function must be callable",
                ));
            }
            Ok(())
        })?;

        // Create a Box<dyn Fn> to match the expected function type
        let rust_postprocess =
            Box::new(move |entry: BsxBatch| -> anyhow::Result<BsxBatch> {
                let py_batch_res = Python::with_gil(|py| -> PyResult<PyBsxBatch> {
                    let py_batch = PyBsxBatch::from(entry);
                    let args = PyTuple::new_bound(py, &[py_batch.into_py(py)]);
                    let result = postprocess_fn.call1(py, args)?;
                    let py_batch = result.extract::<PyBsxBatch>(py)?;
                    Ok(py_batch)
                });
                py_batch_res
                    .map_err(|e| anyhow::anyhow!("{}", e))
                    .map(|batch| batch.into())
            });

        let result = self.inner.query(contig.into(), None);

        match result {
            Ok(Some(batch)) => Ok(Some(PyDataFrame(batch.into_inner()))),
            Ok(None) => Ok(None),
            Err(e) => Err(PyErr::from(e)),
        }
    }

    fn reset(&mut self) {
        self.inner.reset();
    }

    fn chr_order(&self) -> Vec<String> {
        self.inner
            .index()
            .get_chr_order()
            .iter()
            .cloned()
            .map(|chr| chr.to_string())
            .collect()
    }

    fn set_preprocess_fn(
        &mut self,
        preprocess_fn: Py<PyAny>,
    ) -> PyResult<()> {
        Python::with_gil(|py| {
            if !preprocess_fn.clone_ref(py).into_bound(py).is_callable() {
                return Err(PyValueError::new_err(
                    "Preprocess function must be callable",
                ));
            }
            Ok(())
        })?;

        // Create a Box<dyn Fn> to match the expected function type
        let rust_preprocess_fn =
            Box::new(move |entry: BsxBatch| -> anyhow::Result<BsxBatch> {
                let py_batch_res = Python::with_gil(|py| -> PyResult<PyBsxBatch> {
                    let py_batch = PyBsxBatch::from(entry);
                    let args = PyTuple::new_bound(py, &[py_batch.into_py(py)]);
                    let result = preprocess_fn.call1(py, args)?;
                    let py_batch = result.extract::<PyBsxBatch>(py)?;
                    Ok(py_batch)
                });
                py_batch_res
                    .map_err(|e| anyhow::anyhow!("{}", e))
                    .map(|batch| batch.into())
            });

        self.inner.set_preprocess_fn(Some(rust_preprocess_fn));
        Ok(())
    }

    fn index(&self) -> PyBatchIndex {
        self.inner.index().clone().into()
    }
}
