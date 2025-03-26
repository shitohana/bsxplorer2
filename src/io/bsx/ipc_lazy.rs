use std::io::{Read, Seek, SeekFrom};

use anyhow::anyhow;
use polars::export::arrow::io::ipc::read::{read_file_dictionaries,
                                           read_file_metadata,
                                           Dictionaries,
                                           FileMetadata,
                                           OutOfSpecKind};
use polars::prelude::*;
use polars_arrow_format::ipc::planus::ReadAsRoot;

pub struct LazyIpcFileReader<R: Read + Seek> {
    handle:          R,
    metadata:        FileMetadata,
    dictionaries:    Option<Dictionaries>,
    current_block:   usize,
    blocks_total:    usize,
    data_scratch:    Vec<u8>,
    message_scratch: Vec<u8>,
}

impl<R: Read + Seek> LazyIpcFileReader<R> {
    pub fn try_new(mut handle: R) -> anyhow::Result<Self> {
        let metadata = read_file_metadata(&mut handle).map_err(|e| {
            anyhow!("Could not read IPC file metadata").context(e)
        })?;
        let blocks_total = metadata.blocks.len();

        Ok(Self {
            handle,
            metadata,
            blocks_total,
            current_block: 0,
            dictionaries: Default::default(),
            data_scratch: Default::default(),
            message_scratch: Default::default(),
        })
    }

    pub fn arrow_schema(&self) -> &ArrowSchema { &self.metadata.schema }

    pub fn polars_schema(&self) -> Schema {
        Schema::from_arrow_schema(self.arrow_schema())
    }

    pub fn metadata(&self) -> &FileMetadata { &self.metadata }

    pub fn blocks_total(&self) -> usize { self.blocks_total }

    pub fn current_block(&self) -> usize { self.current_block }

    fn read_dictionaries(&mut self) -> PolarsResult<()> {
        if self.dictionaries.is_none() {
            self.dictionaries = Some(read_file_dictionaries(
                &mut self.handle,
                &self.metadata,
                &mut self.data_scratch,
            )?);
        };
        Ok(())
    }

    fn get_message_from_block_offset(
        &mut self,
        offset: u64,
    ) -> anyhow::Result<polars_arrow_format::ipc::MessageRef> {
        self.handle
            .seek(SeekFrom::Start(offset))?;
        let mut meta_buf = [0; 4];
        self.handle.read_exact(&mut meta_buf)?;
        if meta_buf == CONTINUATION_MARKER {
            // continuation marker encountered, read message next
            self.handle.read_exact(&mut meta_buf)?;
        }
        let meta_len: usize = i32::from_le_bytes(meta_buf)
            .try_into()
            .map_err(|_| anyhow!(OutOfSpecKind::UnexpectedNegativeInteger))?;

        self.message_scratch.clear();
        self.message_scratch
            .try_reserve(meta_len)?;
        self.handle
            .by_ref()
            .take(meta_len as u64)
            .read_to_end(&mut self.message_scratch)?;

        Ok(polars_arrow_format::ipc::MessageRef::read_as_root(
            &self.message_scratch,
        )
        .map_err(|err| anyhow!(OutOfSpecKind::InvalidFlatbufferMessage(err)))?)
    }

    fn get_record_batch(
        message: polars_arrow_format::ipc::MessageRef
    ) -> anyhow::Result<polars_arrow_format::ipc::RecordBatchRef> {
        let header = message
            .header()
            .map_err(|err| {
                polars_err!(oos = OutOfSpecKind::InvalidFlatbufferHeader(err))
            })?
            .ok_or_else(|| {
                polars_err!(oos = OutOfSpecKind::MissingMessageHeader)
            })?;
        match header {
            polars_arrow_format::ipc::MessageHeaderRef::RecordBatch(batch) => {
                Ok(batch)
            },
            _ => Err(anyhow!(OutOfSpecKind::UnexpectedMessageType)),
        }
    }
}

const CONTINUATION_MARKER: [u8; 4] = [0xFF; 4];
