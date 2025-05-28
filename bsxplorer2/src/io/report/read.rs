use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::{
    Path,
    PathBuf,
};
use std::thread::JoinHandle;

use anyhow::bail;
use bio::io::fasta::{
    Reader as FastaReader,
    Record as FastaRecord,
};
use crossbeam::channel::Receiver;
use hashbrown::HashMap;
use noodles::fasta::fai::{
    Reader as FaiReader,
    Record,
};
use noodles::fasta::index as index_fasta;
use polars::io::mmap::MmapBytesReader;
use polars::prelude::*;
use rayon::prelude::*;

use crate::prelude::*;
use crate::utils::{
    self,
    get_categorical_dtype,
    THREAD_POOL,
};
use crate::with_field_fn;

pub(crate) fn read_chrom_names<P: AsRef<Path>>(
    path: P,
    is_index: bool,
) -> std::io::Result<Vec<String>> {
    let index = if is_index {
        FaiReader::new(BufReader::new(File::open(path)?)).read_index()?
    }
    else {
        index_fasta(path)?
    };
    let records: Vec<Record> = index.into();
    Ok(records
        .into_iter()
        .map(|r| String::from_utf8_lossy(r.name()).to_string())
        .collect())
}

/// A wrapper around BatchedCsvReader that manages ownership of the reader
pub struct OwnedBatchedCsvReader<F>
where
    F: MmapBytesReader, {
    #[allow(dead_code)]
    // this exists because we need to keep ownership
    /// Schema of the CSV file
    pub schema:         SchemaRef,
    /// The batched reader for the CSV file
    pub batched_reader: BatchedCsvReader<'static>,
    // keep ownership
    /// Original CSV reader
    pub _reader:        CsvReader<F>,
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

/// Builder for `ReportReader` to configure its behavior.
pub struct ReportReaderBuilder {
    report_type: ReportType,
    chunk_size:  usize,
    fasta_path:  Option<PathBuf>,
    fai_path:    Option<PathBuf>,
    batch_size:  usize,
    low_memory:  bool,
    #[cfg(feature = "compression")]
    compression: Option<Compression>,
}

impl Default for ReportReaderBuilder {
    fn default() -> Self {
        Self {
            report_type: ReportType::Bismark,
            chunk_size: 10000,
            fasta_path: None,
            fai_path: None,
            batch_size: 100000,
            low_memory: false,
            #[cfg(feature = "compression")]
            compression: None,
        }
    }
}

impl ReportReaderBuilder {
    with_field_fn!(report_type, ReportType);

    with_field_fn!(chunk_size, usize);

    with_field_fn!(fasta_path, Option<PathBuf>);

    with_field_fn!(fai_path, Option<PathBuf>);

    with_field_fn!(batch_size, usize);

    with_field_fn!(low_memory, bool);

    #[cfg(feature = "compression")]
    with_field_fn!(compression, Option<Compression>);
}

impl ReportReaderBuilder {
    /// Determines the data type of the chromosome column based on FASTA/FAI.
    fn get_chr_dtype(&self) -> anyhow::Result<Option<DataType>> {
        let chroms = match (self.fasta_path.as_ref(), self.fai_path.as_ref()) {
            (_, Some(fai)) => Some(read_chrom_names(fai, true)?),
            (Some(fasta), None) => Some(read_chrom_names(fasta, false)?),
            (None, None) => return Ok(None),
        };

        let dtype = chroms.map(get_categorical_dtype);
        Ok(dtype)
    }

    /// Creates an iterator over the FASTA records.
    fn get_fasta_iterator(
        &self
    ) -> anyhow::Result<Option<Box<dyn Iterator<Item = FastaRecord>>>> {
        if let Some(fasta_path) = self.fasta_path.as_ref() {
            let reader = FastaReader::new(BufReader::new(File::open(fasta_path)?));
            let iterator = Some(Box::new(
                reader.records().map(|r| r.expect("Failed to read record")),
            ) as Box<dyn Iterator<Item = FastaRecord>>);
            Ok(iterator)
        }
        else {
            Ok(None)
        }
    }

    /// Opens the file and returns a `MmapBytesReader`.
    fn get_file_handle(
        &self,
        path: PathBuf,
    ) -> anyhow::Result<Box<dyn MmapBytesReader>> {
        #[cfg(feature = "compression")]
        {
            if let Some(compression) = &self.compression {
                let file = File::open(path)?;
                return compression.get_decoder(file);
            }
        }

        // No compression or compression feature not enabled
        Ok(Box::new(File::open(path)?))
    }

    /// Creates a `OwnedBatchedCsvReader` from a `MmapBytesReader`.
    fn get_csv_reader(
        &self,
        handle: Box<dyn MmapBytesReader>,
    ) -> anyhow::Result<OwnedBatchedCsvReader<Box<dyn MmapBytesReader>>> {
        let csv_reader = self
            .report_type
            .read_options()
            .with_low_memory(self.low_memory)
            .with_chunk_size(self.batch_size)
            .into_reader_with_file_handle(handle);

        let owned_batched =
            OwnedBatchedCsvReader::new(csv_reader, self.report_type.schema().into());
        Ok(owned_batched)
    }

    /// Builds a `ReportReader` from a `MmapBytesReader`.
    pub fn build_from_handle(
        self,
        handle: Box<dyn MmapBytesReader>,
    ) -> anyhow::Result<ReportReader> {
        use ReportType as RS;

        let chr_dtype = self.get_chr_dtype()?;
        if chr_dtype.is_none() {
            bail!("Either Fasta path or Fai path must be specified")
        }
        let fasta_reader = self.get_fasta_iterator()?;
        if matches!(self.report_type, RS::Coverage | RS::BedGraph)
            && fasta_reader.is_none()
        {
            bail!("Fasta path must be specified for Bedgraph or Coverage reading")
        }

        let mut csv_reader = self.get_csv_reader(handle)?;

        let threads = utils::n_threads();
        let (sender, receiver) = crossbeam::channel::bounded(threads * 2);

        let join_handle = std::thread::spawn(move || {
            while let Ok(Some(mut df)) = csv_reader.next_batches(threads) {
                df.drain(..)
                    .for_each(|df| sender.send(df).expect("Failed to send data"));
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

    /// Builds a `ReportReader` from a file path.
    pub fn build(
        self,
        path: PathBuf,
    ) -> anyhow::Result<ReportReader> {
        let handle = self.get_file_handle(path)?;
        self.build_from_handle(handle)
    }
}

/// Reads report data and yields batches.
pub struct ReportReader {
    _join_handle:  JoinHandle<()>,
    data_receiver: Receiver<DataFrame>,
    report_type:   ReportType,
    cached_batch:  BTreeMap<usize, BsxBatch>,
    seen_chr:      HashMap<BsxSmallStr, usize>,
    chunk_size:    usize,

    fasta_reader: Option<Box<dyn Iterator<Item = FastaRecord>>>,
    cached_chr:   Option<(BsxSmallStr, ContextData)>,
    chr_dtype:    Option<DataType>,
}

impl ReportReader {
    /// Tries to take a cached batch.
    fn take_cached(
        &mut self,
        force: bool,
    ) -> Option<BsxBatch> {
        // No cache at all
        if self.cached_batch.is_empty() {
            return None;
        }
        // If there is more than one chromosome cached, force yield to avoid
        // blocking on others
        let should_force = force || self.cached_batch.len() > 1;
        // As we've checked, that cache is not empty, we can take first
        let (idx, batch) = self.cached_batch.pop_first()?;
        // Split by chunk size
        let (first, second) = batch.split_at(self.chunk_size);
        // If cached length was less than chunk_size
        if second.is_empty() && !should_force {
            // Return None, if we do not force yield and expect data
            // continuation
            self.cached_batch.insert(idx, first);
            return None;
        // If cached length was greater than chunk size
        }
        else if !second.is_empty() {
            // Add the remainder back to cache
            self.cached_batch.insert(idx, second);
        }
        // If we have not returned None yet we can yield first
        Some(first)
    }

    /// Aligns a batch with context data.
    fn align_batch(
        &mut self,
        batch: BsxBatch,
        is_final: bool,
    ) -> anyhow::Result<BsxBatch> {
        // If some cached chr data
        if let Some((chr, mut cached_data)) = Option::take(&mut self.cached_chr) {
            // If cached different chromosome raise
            if batch.seqname().unwrap_or_default() != chr.as_str() {
                bail!("Chromosome mismatch")
            // Else align
            }
            else {
                // If it is a final batch - drain cached chr
                let (context_data, new_cache) = if is_final {
                    (cached_data, None)
                // If not final, take till end of batch
                }
                else {
                    let drained =
                        cached_data.drain_until(batch.last_pos().unwrap_or(0));
                    (drained, Some((chr, cached_data)))
                };
                // Update cache (leftover if not final else None)
                self.cached_chr = new_cache;
                let _chr_val = batch.seqname().unwrap_or_default().to_string();
                // Align batch
                let aligned = batch.add_context_data(context_data)?;
                Ok(aligned)
            }
        // If cache is empty read next chromosome
        }
        else {
            // Try read
            if let Some(new_sequence) = self.fasta_reader.as_mut().unwrap().next() {
                // Process sequence
                let mut new_context_data = ContextData::empty();
                new_context_data.read_sequence(new_sequence.seq(), 1);
                self.cached_chr = Some((new_sequence.id().into(), new_context_data));
                // Try aligning again
                Ok(self
                    .align_batch(batch, is_final)
                    .map_err(|e| anyhow::anyhow!(e))?)
            // Everything is read, but more sequence requested. Raise
            }
            else {
                bail!("Sequence has already been fully written")
            }
        }
    }

    /// Fills the cache with more data.
    fn fill_cache(&mut self) -> anyhow::Result<()> {
        let n_threads = THREAD_POOL.current_num_threads();
        let mut new_batches = Vec::with_capacity(n_threads);
        let mut stream_ended = false;

        while !self.data_receiver.is_empty() && new_batches.len() < n_threads {
            let new_batch = self.data_receiver.recv()?;
            stream_ended =
                self.data_receiver.is_empty() && self._join_handle.is_finished();

            new_batches.push(new_batch);
        }

        let converted = THREAD_POOL.install(|| {
            new_batches
                .into_par_iter()
                .flat_map(|df| {
                    df.partition_by_stable([self.report_type.chr_col()], true)
                })
                .flatten()
                .filter(|df| !df.is_empty())
                .map(|df| {
                    BsxBatchBuilder::all_checks()
                        .with_check_single_chr(false)
                        .with_chr_dtype(self.chr_dtype.clone())
                        .build_from_report(df, self.report_type)
                })
                .collect::<Result<Vec<_>, _>>()
        })?;

        let mut is_last = vec![false; converted.len()];
        let mut prev_chr = self.cached_chr.as_ref().map(|v| v.0.clone());
        let mut prev_end = 0;
        for i in 0..converted.len() {
            let batch_contig = converted[i].as_contig().unwrap();

            if i == 0
                && prev_chr.is_some()
                && prev_chr.as_ref() != Some(batch_contig.seqname())
            {
                is_last[i] = true;
            }
            else if i > 0 && prev_chr.as_ref() != Some(batch_contig.seqname()) {
                is_last[i - 1] = true;
            }
            else if i == converted.len() - 1 && stream_ended {
                is_last[i] = true;
            }

            if prev_chr.as_ref() == Some(batch_contig.seqname())
                && prev_end >= batch_contig.start()
            {
                return Err(anyhow::anyhow!(
                    "Overlapping records found at position {}",
                    batch_contig
                ));
            }

            prev_chr = Some(batch_contig.seqname().clone());
            prev_end = batch_contig.end();
        }

        for (mut bsx_batch, last_one) in converted.into_iter().zip(is_last.into_iter())
        {
            // Align batch if required
            if self.report_type.need_align() {
                bsx_batch = self.align_batch(bsx_batch, last_one)?;
            }

            // Add chromosome to seen and retrieve chr index
            let batch_chr = bsx_batch.seqname().unwrap_or_default().into();
            let seen_chr_len = self.seen_chr.len();
            let chr_idx = self.seen_chr.entry(batch_chr).or_insert(seen_chr_len);
            // Update cached batch
            self.cached_batch
                .entry(*chr_idx)
                .and_modify(|b| {
                    unsafe { b.extend_unchecked(&bsx_batch).expect("vstack failed") };
                })
                .or_insert(bsx_batch);
        }

        self.cached_batch
            .values_mut()
            .for_each(|batch| batch.rechunk());

        Ok(())
    }

    /// Sets the FASTA reader.
    pub fn set_fasta_reader(
        &mut self,
        fasta_reader: Option<Box<dyn Iterator<Item = FastaRecord>>>,
    ) {
        self.fasta_reader = fasta_reader;
    }
}

impl Iterator for ReportReader {
    type Item = anyhow::Result<BsxBatch>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Try to get a cached batch first
            if let Some(batch) = self.take_cached(false) {
                return Some(Ok(batch));
            }
            // If the reader thread is still alive or we can expect more data
            if !self._join_handle.is_finished() || !self.data_receiver.is_empty() {
                if let Err(e) = self.fill_cache() {
                    return Some(Err(e));
                }
                // After filling, try to take again
                continue;
            }
            // Reader is done and no more data is expected. Drain remaining
            // cache.
            if let Some(batch) = self.take_cached(true) {
                return Some(Ok(batch));
            }
            // All done
            return None;
        }
    }
}
