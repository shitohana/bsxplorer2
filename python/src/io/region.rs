use std::fs::File;

use bsxplorer2::data_structs::coords::Contig;
use bsxplorer2::data_structs::enums::Strand;
use bsxplorer2::io::bsx::{BatchIndex, BsxFileReader, RegionReader};
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

/// RegionReader is a reader for BSX files that operates on a specific region of
/// the genome. (Python wrapper)
#[pyclass(unsendable, name = "RegionReader")] // Mark as unsendable as it holds File/Reader which is not Sync/Send
pub struct PyRegionReader {
    reader: RegionReader<File, String, u32>,
}

#[pymethods]
impl PyRegionReader {
    #[new]
    fn new(path: String) -> PyResult<Self> {
        let file = File::open(path).map_err(PyErr::from)?;
        // The BsxFileReader needs Read + Seek. std::fs::File provides this.
        let mut bsx_reader: BsxFileReader<File> = BsxFileReader::new(file);
        // Indexing requires Read + Seek. The index() method consumes &mut self,
        // so we need the reader to be mutable.
        let index: BatchIndex<String, u32> = bsx_reader
            .index()
            .map_err(PyErr::from)?
            .clone();

        Ok(Self {
            // The RegionReader also needs Read + Seek. We pass the
            // BsxFileReader itself. Since BsxFileReader wraps the
            // File, it fulfills the requirement. The generic
            // parameters S and P are fixed to String and u32.
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

        // The Rust query method takes an optional postprocess_fn (a fn
        // pointer). We omit this from the Python interface for
        // simplicity.
        let result = self.reader.query(contig, None);

        match result {
            Ok(Some(batch)) => {
                // EncodedBsxBatch wraps a DataFrame. Convert it to PyDataFrame.
                Ok(Some(PyDataFrame(batch.into())))
            },
            Ok(None) => Ok(None),
            Err(e) => Err(PyErr::from(e)),
        }
    }

    /// Resets the cache.
    fn reset(&mut self) { self.reader.reset(); }

    /// Returns the chromosome order as a list of strings.
    fn chr_order(&self) -> Vec<String> {
        self.reader
            .index()
            .chr_order()
            .iter()
            .cloned()
            .collect()
    }
}
