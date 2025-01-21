pub mod ipc;
pub mod read;

pub mod writer {
    use crate::io::report::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
    use itertools::Itertools;
    use polars::datatypes::CategoricalOrdering;
    use polars::io::ipc::{BatchedWriter, IpcWriterOptions};
    use polars::prelude::*;
    use std::error::Error;
    use std::io::Write;
    use std::path::PathBuf;
    use std::sync::Arc;
    use crate::utils::get_categorical_dtype;
    
    // TODO Add marker for no context information to invalidate Nulls in context
    //      column
    pub struct BsxIpcWriter<W>
    where
        W: Write 
    {
        writer: BatchedWriter<W>,
        schema: Schema
    }

    impl<W> BsxIpcWriter<W> where W: Write {
        pub fn try_new(sink: W, chr_names: Vec<String>) -> Result<BsxIpcWriter<W>, Box<dyn Error>> {
            let opts = IpcWriterOptions {
                compression: Some(IpcCompression::default()),
                maintain_order: true,
            };
            let chr_dtype = get_categorical_dtype(chr_names);
            let schema = EncodedBsxBatch::get_schema(&chr_dtype);

            let writer = opts.to_writer(sink).batched(&schema)?;
            Ok(Self {writer, schema})
        }
        
        pub fn try_from_sink_and_fai(
            sink: W,
            fai_path: PathBuf,
        ) -> Result<Self, Box<dyn Error>> {
            let index = bio::io::fasta::Index::from_file(&fai_path)?;
            let chr_names= index.sequences().into_iter().map(|seq| seq.name).collect_vec();
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
}