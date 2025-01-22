use crate::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use crate::utils::get_categorical_dtype;
use itertools::Itertools;
use polars::datatypes::DataType;
use polars::error::PolarsResult;
use polars::prelude::{IpcCompression, IpcWriterOptions, Schema};
use std::error::Error;
use std::io::Write;
use std::path::PathBuf;

// TODO Add marker for no context information to invalidate Nulls in context
//      column
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
    pub fn try_new(sink: W, chr_names: Vec<String>) -> Result<BsxIpcWriter<W>, Box<dyn Error>> {
        let opts = IpcWriterOptions {
            compression: Some(IpcCompression::default()),
            maintain_order: true,
        };
        let chr_dtype = get_categorical_dtype(chr_names);
        let schema = EncodedBsxBatch::get_schema(&chr_dtype);

        let writer = opts.to_writer(sink).batched(&schema)?;
        Ok(Self { writer, schema })
    }

    pub fn try_from_sink_and_fai(sink: W, fai_path: PathBuf) -> Result<Self, Box<dyn Error>> {
        let index = bio::io::fasta::Index::from_file(&fai_path)?;
        let chr_names = index
            .sequences()
            .into_iter()
            .map(|seq| seq.name)
            .collect_vec();
        Self::try_new(sink, chr_names)
    }

    pub fn write_batch(&mut self, batch: EncodedBsxBatch) -> PolarsResult<()> {
        self.writer.write_batch(batch.data())
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
