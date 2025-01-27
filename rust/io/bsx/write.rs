use crate::bsx_batch::{BsxBatch, BsxBatchMethods, EncodedBsxBatch};
use crate::utils::get_categorical_dtype;
#[cfg(feature = "python")]
use crate::utils::wrap_polars_result;
use itertools::Itertools;
use polars::datatypes::DataType;
use polars::error::PolarsResult;
use polars::export::arrow::datatypes::Metadata;
use polars::prelude::{IpcCompression, IpcWriterOptions, Schema};
#[cfg(feature = "python")]
use pyo3_polars::error::PyPolarsErr;
use std::error::Error;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;

pub struct BsxIpcWriter<W>
where
    W: Write,
{
    writer: polars::io::ipc::BatchedWriter<W>,
    schema: Schema,
}

impl<W> BsxIpcWriter<W>
where
    W: Write,
{
    pub fn try_new(
        sink: W,
        chr_names: Vec<String>,
        compression: Option<IpcCompression>,
        custom_metadata: Option<Arc<Metadata>>,
    ) -> Result<BsxIpcWriter<W>, Box<dyn Error>> {
        let opts = IpcWriterOptions {
            compression,
            maintain_order: true,
        };
        let chr_dtype = get_categorical_dtype(chr_names);
        let schema = EncodedBsxBatch::get_schema(&chr_dtype);

        let mut writer = opts.to_writer(sink);
        if let Some(metadata) = custom_metadata {
            writer.set_custom_schema_metadata(metadata);
        }
        let batched_writer = writer.batched(&schema)?;
        Ok(Self {
            writer: batched_writer,
            schema,
        })
    }

    pub fn try_from_sink_and_fai(
        sink: W,
        fai_path: PathBuf,
        compression: Option<IpcCompression>,
        custom_metadata: Option<Arc<Metadata>>,
    ) -> Result<Self, Box<dyn Error>> {
        let index = bio::io::fasta::Index::from_file(&fai_path)?;
        let chr_names = index
            .sequences()
            .into_iter()
            .map(|seq| seq.name)
            .collect_vec();
        Self::try_new(sink, chr_names, compression, custom_metadata)
    }

    pub fn write_encoded_batch(&mut self, batch: EncodedBsxBatch) -> PolarsResult<()> {
        self.writer.write_batch(batch.data())
    }

    pub fn write_batch(&mut self, batch: BsxBatch) -> PolarsResult<()> {
        let encoded = EncodedBsxBatch::encode(batch, self.get_chr_dtype())?;
        self.write_encoded_batch(encoded)
    }

    pub fn close(&mut self) -> PolarsResult<()> {
        self.writer.finish()
    }

    pub fn get_chr_dtype(&self) -> &DataType {
        self.schema.get("chr").unwrap()
    }
}

impl<W: Write> Drop for BsxIpcWriter<W> {
    fn drop(&mut self) {
        let _ = self.close();
    }
}

#[cfg(feature = "python")]
mod python {
    use super::*;
    use crate::utils::wrap_box_result;
    use crate::wrap_polars_result;
    use pyo3::exceptions::PyIOError;
    use pyo3::prelude::*;
    use std::fs::File;
    use std::io::BufWriter;

    #[pyclass(name = "IpcCompression")]
    #[derive(Clone)]
    pub enum PyIpcCompression {
        ZSTD,
        LZ4,
    }

    impl From<PyIpcCompression> for IpcCompression {
        fn from(value: PyIpcCompression) -> Self {
            match value {
                PyIpcCompression::ZSTD => IpcCompression::ZSTD,
                PyIpcCompression::LZ4 => IpcCompression::LZ4,
            }
        }
    }

    #[pyclass(name = "BsxWriter")]
    pub struct PyBsxIpcWriter {
        inner: BsxIpcWriter<BufWriter<File>>,
    }

    #[pymethods]
    impl PyBsxIpcWriter {
        #[new]
        fn new(
            path: String,
            fai_path: String,
            compression: Option<PyIpcCompression>,
        ) -> PyResult<Self> {
            let res = BsxIpcWriter::try_from_sink_and_fai(
                BufWriter::new(File::create(PathBuf::from(path))?),
                fai_path.into(),
                compression.map(|c| c.into()),
                None,
            );
            let inner = wrap_box_result!(PyIOError, res)?;
            Ok(Self { inner })
        }

        fn write_batch(&mut self, batch: BsxBatch) -> PyResult<()> {
            wrap_polars_result!(self.inner.write_batch(batch))
        }

        fn write_encoded_batch(&mut self, batch: EncodedBsxBatch) -> PyResult<()> {
            wrap_polars_result!(self.inner.write_encoded_batch(batch))
        }

        fn close(&mut self) -> PyResult<()> {
            wrap_polars_result!(self.inner.close())
        }
    }
}

#[cfg(feature = "python")]
pub use python::{PyBsxIpcWriter, PyIpcCompression};
