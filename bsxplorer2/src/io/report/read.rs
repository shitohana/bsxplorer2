use crate::data_structs::batch::{
    BatchType, BsxBatchBuilder, BsxBatchMethods, BsxTypeTag, LazyBsxBatch,
};
use crate::io::read_chrom;
use crate::io::report::schema::ReportTypeSchema;
use anyhow::bail;
use bio::io::fasta::{Reader as FastaReader, Record as FastaRecord};

use crate::data_structs::context_data::ContextData;
#[cfg(feature = "compression")]
use crate::io::compression::Compression;
use crate::utils::get_categorical_dtype;
use crossbeam::channel::Receiver;
use polars::io::mmap::MmapBytesReader;
use polars::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::BufReader;
use std::marker::PhantomData;
use std::path::PathBuf;
use std::thread::JoinHandle;

/// A wrapper around BatchedCsvReader that manages ownership of the reader
pub struct OwnedBatchedCsvReader<F>
where
    F: MmapBytesReader,
{
    #[allow(dead_code)]
    // this exists because we need to keep ownership
    /// Schema of the CSV file
    pub schema: SchemaRef,
    /// The batched reader for the CSV file
    pub batched_reader: BatchedCsvReader<'static>,
    // keep ownership
    /// Original CSV reader
    pub _reader: CsvReader<F>,
}

impl<F> OwnedBatchedCsvReader<F>
where
    F: MmapBytesReader + 'static,
{
    /// Create a new OwnedBatchedCsvReader from a CsvReader and schema
    pub(crate) fn new(
        mut reader: CsvReader<F>,
        schema: Arc<Schema>,
    ) -> Self {
        let batched_reader = reader
            .batched_borrowed()
            .expect("Could not create batched CSV reader.");
        let batched_reader: BatchedCsvReader<'static> =
            unsafe { std::mem::transmute(batched_reader) };
        Self {
            batched_reader,
            _reader: reader,
            schema,
        }
    }
}

impl<F> OwnedBatchedCsvReader<F>
where
    F: MmapBytesReader + 'static,
{
    /// Read the next n batches from the CSV file
    pub fn next_batches(
        &mut self,
        n: usize,
    ) -> PolarsResult<Option<Vec<DataFrame>>> {
        self.batched_reader.next_batches(n)
    }
}

pub struct ReportReaderBuilder<B: BsxBatchMethods + BsxTypeTag> {
    report_type: ReportTypeSchema,
    chunk_size: usize,
    fasta_path: Option<PathBuf>,
    fai_path: Option<PathBuf>,
    batch_size: usize,
    n_threads: Option<usize>,
    low_memory: bool,
    queue_len: usize,
    #[cfg(feature = "compression")]
    compression: Option<Compression>,
    _batch_type: PhantomData<B>,
}

impl<B: BsxBatchMethods + BsxTypeTag> Default for ReportReaderBuilder<B> {
    fn default() -> Self {
        Self {
            report_type: ReportTypeSchema::Bismark,
            chunk_size: 10000,
            fasta_path: None,
            fai_path: None,
            batch_size: 100000,
            n_threads: None,
            low_memory: false,
            queue_len: 1000,
            #[cfg(feature = "compression")]
            compression: None,
            _batch_type: PhantomData,
        }
    }
}

impl<B: BsxBatchMethods + BsxTypeTag> ReportReaderBuilder<B> {
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn with_report_type(
        mut self,
        report_type: ReportTypeSchema,
    ) -> Self {
        self.report_type = report_type;
        self
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn with_chunk_size(
        mut self,
        chunk_size: usize,
    ) -> Self {
        self.chunk_size = chunk_size;
        self
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn with_fasta_path(
        mut self,
        fasta_path: PathBuf,
    ) -> Self {
        self.fasta_path = Some(fasta_path);
        self
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn with_fai_path(
        mut self,
        fai_path: PathBuf,
    ) -> Self {
        self.fai_path = Some(fai_path);
        self
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn with_batch_size(
        mut self,
        batch_size: usize,
    ) -> Self {
        self.batch_size = batch_size;
        self
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn with_n_threads(
        mut self,
        n_threads: usize,
    ) -> Self {
        self.n_threads = Some(n_threads);
        self
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn with_low_memory(
        mut self,
        low_memory: bool,
    ) -> Self {
        self.low_memory = low_memory;
        self
    }

    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn with_queue_len(
        mut self,
        queue_len: usize,
    ) -> Self {
        self.queue_len = queue_len;
        self
    }

    #[cfg(feature = "compression")]
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn with_compression(
        mut self,
        compression: Compression,
    ) -> Self {
        self.compression = Some(compression);
        self
    }
}

impl<B: BsxBatchMethods + BsxTypeTag> ReportReaderBuilder<B> {
    fn get_chr_dtype(&self) -> anyhow::Result<Option<DataType>> {
        let chroms = match (self.fasta_path.as_ref(), self.fai_path.as_ref()) {
            (_, Some(fai)) => Some(read_chrom(fai, true)?),
            (Some(fasta), None) => Some(read_chrom(fasta, false)?),
            (None, None) => return Ok(None),
        };

        let dtype = chroms.map(|v| get_categorical_dtype(v));
        Ok(dtype)
    }

    fn get_fasta_iterator(
        &self
    ) -> anyhow::Result<Option<Box<dyn Iterator<Item = FastaRecord>>>> {
        if let Some(fasta_path) = self.fasta_path.as_ref() {
            let reader =
                FastaReader::new(BufReader::new(File::open(fasta_path)?));
            let iterator = Some(Box::new(
                reader
                    .records()
                    .map(|r| r.expect("Failed to read record")),
            )
                as Box<dyn Iterator<Item = FastaRecord>>);
            Ok(iterator)
        } else {
            Ok(None)
        }
    }

    fn get_file_handle(
        &self,
        path: PathBuf,
    ) -> anyhow::Result<Box<dyn MmapBytesReader>> {
        #[cfg(feature = "compression")]
        {
            if let Some(compression) = &self.compression {
                let file = File::open(path)?;
                return Ok(compression.get_decoder(file)?);
            }
        }

        // No compression or compression feature not enabled
        Ok(Box::new(File::open(path)?))
    }

    fn get_csv_reader(
        &self,
        handle: Box<dyn MmapBytesReader>,
    ) -> anyhow::Result<OwnedBatchedCsvReader<Box<dyn MmapBytesReader>>> {
        let csv_reader = self
            .report_type
            .read_options()
            .with_n_threads(self.n_threads)
            .with_low_memory(self.low_memory)
            .with_chunk_size(self.batch_size)
            .into_reader_with_file_handle(handle);

        let owned_batched = OwnedBatchedCsvReader::new(
            csv_reader,
            self.report_type.schema().into(),
        );
        Ok(owned_batched)
    }

    pub fn build_from_handle(
        self,
        handle: Box<dyn MmapBytesReader>,
    ) -> anyhow::Result<ReportReader<B>> {
        use ReportTypeSchema as RS;

        let chr_dtype = self.get_chr_dtype()?;
        if matches!(B::type_enum(), BatchType::Encoded) && chr_dtype.is_none() {
            bail!("Either Fasta path or Fai path must be specified for Encoded reading")
        }
        let fasta_reader = self.get_fasta_iterator()?;
        if matches!(self.report_type, RS::Coverage | RS::BedGraph)
            && fasta_reader.is_none()
        {
            bail!(
                "Fasta path must be specified for Bedgraph or Coverage reading"
            )
        }

        let mut csv_reader = self.get_csv_reader(handle)?;
        let (sender, receiver) = crossbeam::channel::bounded(self.queue_len);

        let join_handle = std::thread::spawn(move || {
            while let Ok(Some(mut df)) = csv_reader.next_batches(1) {
                if df.len() != 1 {
                    panic!("Unexpected batch count {}", df.len());
                }
                sender
                    .send(df.pop().unwrap())
                    .expect("Failed to send data");
            }
        });

        let reader = ReportReader {
            _join_handle: join_handle,
            data_receiver: receiver,
            report_type: self.report_type,
            cached_batch: BTreeMap::new(),
            seen_chr: HashMap::new(),
            chunk_size: self.chunk_size,
            fasta_reader,
            cached_chr: None,
            chr_dtype,
        };
        Ok(reader)
    }

    pub fn build(
        self,
        path: PathBuf,
    ) -> anyhow::Result<ReportReader<B>> {
        let handle = self.get_file_handle(path)?;
        self.build_from_handle(handle)
    }
}

pub struct ReportReader<B: BsxBatchMethods + BsxTypeTag> {
    _join_handle: JoinHandle<()>,
    data_receiver: Receiver<DataFrame>,
    report_type: ReportTypeSchema,
    cached_batch: BTreeMap<usize, B>,
    seen_chr: HashMap<String, usize>,
    chunk_size: usize,

    fasta_reader: Option<Box<dyn Iterator<Item = FastaRecord>>>,
    cached_chr: Option<(String, ContextData)>,
    chr_dtype: Option<DataType>,
}

impl<B: BsxBatchMethods + BsxTypeTag> ReportReader<B> {
    fn take_cached(
        &mut self,
        force: bool,
    ) -> Option<B> {
        // No cache at all
        if self.cached_batch.is_empty() {
            return None;
        }
        // If there is more than one chromosome cached, force yield to avoid blocking on others
        let should_force = force || self.cached_batch.len() > 1;
        // As we've checked, that cache is not empty, we can take first
        let (idx, batch) = self.cached_batch.pop_first()?;
        // Split by chunk size
        let (first, second) = batch.split_at(self.chunk_size);
        // If cached length was less than chunk_size
        if second.is_empty() && !should_force {
            // Return None, if we do not force yield and expect data continuation
            self.cached_batch.insert(idx, first);
            return None;
        // If cached length was greater than chunk size
        } else if !second.is_empty() {
            // Add the remainder back to cache
            self.cached_batch.insert(idx, second);
        }
        // If we have not returned None yet we can yield first
        Some(first)
    }

    fn align_batch(
        &mut self,
        batch: B,
        is_final: bool,
    ) -> anyhow::Result<B> {
        // If some cached chr data
        if let Some((chr, mut cached_data)) = Option::take(&mut self.cached_chr)
        {
            // If cached different chromosome raise
            if batch.chr_val()? != chr {
                bail!("Chromosome mismatch")
            // Else align
            } else {
                // If it is a final batch - drain cached chr
                let (context_data, new_cache) = if is_final {
                    (cached_data, None)
                // If not final, take till end of batch
                } else {
                    let drained =
                        cached_data.drain_until(batch.end_pos().unwrap());
                    (drained, Some((chr, cached_data)))
                };
                // Update cache (leftover if not final else None)
                self.cached_chr = new_cache;
                let chr_val = batch.chr_val()?.to_string();
                // Align batch
                let aligned = LazyBsxBatch::from(batch)
                    .align_with_contexts(context_data, &chr_val)
                    .collect()?;
                Ok(aligned)
            }
        // If cache is empty read next chromosome
        } else {
            // Try read
            if let Some(new_sequence) = self
                .fasta_reader
                .as_mut()
                .unwrap()
                .next()
            {
                // Process sequence
                let mut new_context_data = ContextData::empty();
                new_context_data.read_sequence(new_sequence.seq(), 1);
                self.cached_chr =
                    Some((new_sequence.id().to_string(), new_context_data));
                // Try aligning again
                Ok(self
                    .align_batch(batch, is_final)
                    .map_err(|e| anyhow::anyhow!(e))?)
            // Everything is read, but more sequence requested. Raise
            } else {
                bail!("Sequence has already been fully written")
            }
        }
    }

    fn fill_cache(&mut self) -> anyhow::Result<()> {
        // Try receive another batch
        let new_batch = self.data_receiver.recv().map_err(|e| {
            anyhow::anyhow!("Channel closed unexpectedly while reading a batch")
                .context(e)
        })?;
        let last_one =
            self.data_receiver.is_empty() && self._join_handle.is_finished();
        // Partition by chromosome
        let partitioned = new_batch
            .partition_by_stable([self.report_type.chr_col()], true)
            .map_err(|e| {
                anyhow::anyhow!("Failed to partition batch").context(e)
            })?;
        // 1 if single chromosome or > 1 if multiple
        let n_partitioned = partitioned.len();

        for (index, data) in partitioned.into_iter().enumerate() {
            // Convert to BSX batch format
            let mut bsx_batch = BsxBatchBuilder::all_checks()
                .with_report_type(self.report_type)
                .with_check_single_chr(false)
                .with_chr_dtype(self.chr_dtype.clone()) // Will raise if not specified
                .build::<B>(data)?;

            // Align batch if required
            if self.report_type.need_align() {
                let is_final = (index < n_partitioned - 1) || last_one;
                bsx_batch = self.align_batch(bsx_batch, is_final)?;
            }

            // Add chromosome to seen and retrieve chr index
            let batch_chr = bsx_batch.chr_val()?.to_owned();
            let seen_chr_len = self.seen_chr.len();
            let chr_idx = self
                .seen_chr
                .entry(batch_chr)
                .or_insert(seen_chr_len);
            // Update cached batch
            self.cached_batch
                .entry(*chr_idx)
                .and_modify(|b| {
                    b.extend(&bsx_batch)
                        .expect("vstack failed");
                })
                .or_insert(bsx_batch);
        }

        Ok(())
    }

    pub fn set_fasta_reader(
        &mut self,
        fasta_reader: Option<Box<dyn Iterator<Item = FastaRecord>>>,
    ) {
        self.fasta_reader = fasta_reader;
    }
}

impl<B: BsxBatchMethods + BsxTypeTag> Iterator for ReportReader<B> {
    type Item = anyhow::Result<B>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Try to get a cached batch first
            if let Some(batch) = self.take_cached(false) {
                return Some(Ok(batch));
            }
            // If the reader thread is still alive or we can expect more data
            if !self._join_handle.is_finished()
                || !self.data_receiver.is_empty()
            {
                if let Err(e) = self.fill_cache() {
                    return Some(Err(e));
                }
                // After filling, try to take again
                continue;
            }
            // Reader is done and no more data is expected. Drain remaining cache.
            if let Some(batch) = self.take_cached(true) {
                return Some(Ok(batch));
            }
            // All done
            return None;
        }
    }
}
