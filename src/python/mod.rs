#[cfg(feature = "python")]
pub mod python {
    use crate::bsxipc::{BSXFile, BSXIterator, BSXKey};
    use crate::read::ReportType;
    use crate::utils::types::{Context, Strand};
    use polars::frame::DataFrame;
    use polars::prelude::SchemaRef;
    use pyo3::exceptions::{PyNotImplementedError, PyRuntimeError};
    use pyo3::prelude::*;
    use pyo3_polars::{PyDataFrame, PySchema};
    use std::fs::File;
    
    
    #[pyclass]
    pub struct BSXFilePy {
        inner: BSXFile<File>,
    }

    #[pymethods]
    impl BSXFilePy {
        #[new]
        pub fn new(path: String, projection: Option<Vec<usize>>) -> Self {
            let file = File::open(path).expect("Failed to open file");
            let inner = BSXFile::new(file, projection).expect("Failed to open BSX file");
            Self { inner }
        }

        pub fn find_index(
            &self,
            chr: String,
            strand: String,
            context: String,
            start: Option<u64>,
            end: Option<u64>,
        ) -> Option<Vec<u32>> {
            let key = BSXKey(
                chr,
                Strand::from_str(strand.as_str()),
                Context::from_str(context.as_str()),
            );
            let range = (start, end);

            self.inner.find_index(key, range)
        }

        pub fn add_index(&self, annotation: PyDataFrame, context: String) -> PyResult<PyDataFrame> {
            let context = Context::from_str(context.as_str());
            let mut annotation: DataFrame = annotation.into();

            match self
                .inner
                .add_index(&mut annotation, Some(context))
            {
                Ok(df) => Ok(df.into()),
                Err(_) => Err(PyRuntimeError::new_err("Failed to add BSX file annotation")),
            }
        }

        #[getter]
        pub fn get_batch_num(&self) -> PyResult<usize> {
            self.inner.get_batch_num().into()
        }

        #[getter]
        pub fn get_polars_schema(&self) -> PyResult<PySchema> {
            let schema = self.inner.get_polars_schema();
            let pyschema = PySchema(SchemaRef::new(schema));
            pyschema.into()
        }
    }

    #[pyclass]
    pub struct BSXIteratorPy {
        inner: BSXIterator<File>,
    }

    #[pymethods]
    impl BSXIteratorPy {
        #[new]
        pub fn new(bsxfile_py: BSXFilePy) -> Self {
            let inner = bsxfile_py.inner.into_iter();
            Self { inner }
        }

        pub fn get_batch(&mut self, index: u32) -> PyResult<PyDataFrame> {
            let res = self.inner.get_batch(index);
            match res {
                Ok(df) => Ok(df.into()),
                Err(_) => Err(PyRuntimeError::new_err("Failed to get batch")),
            }
        }

        pub fn get_region(
            &mut self,
            indices: &[u32],
            start: u64,
            end: u64,
        ) -> PyResult<PyDataFrame> {
            let res = self
                .inner
                .get_region(indices, start, end);
            match res {
                Ok(df) => df.into(),
                Err(_) => Err(PyRuntimeError::new_err("Failed to get region")),
            }
        }

        pub fn iter_regions(
            &mut self,
            indexed_annotation: PyDataFrame,
        ) -> impl Iterator<Item = PyResult<PyDataFrame>> {
            let df: DataFrame = indexed_annotation.into();
            let iterator: impl Iterator<Item = PyResult<PyDataFrame>> = self
                .inner
                .iter_regions(&df)
                .map(|df| -> PyResult<PyDataFrame> {
                    match df {
                        Ok(df) => Ok(df.into()),
                        Err(_) => Err(PyRuntimeError::new_err("Failed to get region df")),
                    }
                });
            iterator
        }
    }

    #[pyfunction]
    fn convert_report(
        input_path: String,
        report_type: String,
        output_path: String,
        fasta_path: String,
        chunk_size: Option<usize>,
        low_memory: Option<bool>,
        n_threads: Option<usize>,
    ) -> PyResult<String> {
        let report_type = match report_type.to_lowercase().as_str() {
            "bismark" => ReportType::BISMARK,
            "bedgraph" => ReportType::BEDGRAPH,
            "coverage" => ReportType::COVERAGE,
            "cgmap" => ReportType::CGMAP,
            _ => {return Err(PyNotImplementedError::new_err("{report_type} is not yet implemented"))}
        };
        report_type.convert_report(
            input_path.as_str(),
            output_path.as_str(),
            fasta_path.as_str(),
            chunk_size,
            low_memory,
            n_threads,
        );
        Ok(output_path)
    }
}
