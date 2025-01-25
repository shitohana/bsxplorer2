use crate::bsx_batch::{BsxBatch, BsxBatchMethods, EncodedBsxBatch};
use crate::utils::get_categorical_dtype;
use itertools::Itertools;
use polars::datatypes::DataType;
use polars::error::PolarsResult;
use polars::prelude::{IpcCompression, IpcWriterOptions, Schema};
use std::error::Error;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;
use polars::export::arrow::datatypes::Metadata;

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
        Ok(Self { writer: batched_writer, schema })
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
