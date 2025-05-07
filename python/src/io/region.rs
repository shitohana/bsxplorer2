use std::fs::File;

use bsxplorer2::data_structs::coords::Contig;
use bsxplorer2::data_structs::enums::Strand;
use bsxplorer2::io::bsx::{BatchIndex, BsxFileReader, RegionReader};
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

#[pyclass(unsendable, name = "RegionReader")]
pub struct PyRegionReader {
    reader: RegionReader<File, String, u32>,
}

#[pymethods]
impl PyRegionReader {
    #[new]
    fn new(path: String) -> PyResult<Self> {
        let file = File::open(path).map_err(PyErr::from)?;
        let mut bsx_reader: BsxFileReader<File> = BsxFileReader::new(file);
        let index: BatchIndex<String, u32> = bsx_reader
            .index()
            .map_err(PyErr::from)?
            .clone();

        Ok(Self {
            reader: RegionReader::new(bsx_reader, index, None),
        })
    }

    fn query(
        &mut self,
        seqname: String,
        start: u32,
        end: u32,
    ) -> PyResult<Option<PyDataFrame>> {
        let contig = Contig::new(seqname, start, end, Strand::None);

        let result = self.reader.query(contig, None);

        match result {
            Ok(Some(batch)) => Ok(Some(PyDataFrame(batch.into()))),
            Ok(None) => Ok(None),
            Err(e) => Err(PyErr::from(e)),
        }
    }

    fn reset(&mut self) { self.reader.reset(); }

    fn chr_order(&self) -> Vec<String> {
        self.reader
            .index()
            .get_chr_order()
            .iter()
            .cloned()
            .collect()
    }
}
