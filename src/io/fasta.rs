#[cfg(feature = "compression")]
use crate::io::compression::Compression;
use noodles::fasta::fai::io::Reader as FaiReader;
use noodles::fasta::fai::Index as FaIndex;
use noodles::fasta::io::indexed_reader::{
    Builder as FastaBuilder, IndexedReader as FastaIndexedReader,
};
use polars::io::mmap::MmapBytesReader;
use std::collections::VecDeque;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

pub(crate) fn get_reader(
    path: PathBuf,
    #[cfg(feature = "compression")]
    compression: Compression,
    index_path: Option<PathBuf>,
) -> anyhow::Result<FastaIndexedReader<BufReader<Box<dyn MmapBytesReader>>>> {
    let index = if let Some(index_path) = index_path {
        Some(FaiReader::new(BufReader::new(File::open(index_path)?)).read_index()?)
    } else {
        None
    };

    let mut builder = FastaBuilder::default();
    if let Some(index) = index {
        builder = builder.set_index(index);
    }

    #[cfg(feature = "compression")]
    let fasta_handle = BufReader::new(compression.get_decoder(File::open(path)?));
    #[cfg(not(feature = "compression"))]
    let fasta_handle = BufReader::new(Box::new(File::open(path)?) as Box<dyn MmapBytesReader>);
    builder
        .build_from_reader(fasta_handle)
        .map_err(|e| e.into())
}

pub(crate) struct FastaReadPlanEntry {
    name: String,
    length: u32,
    read: u32,
}

impl FastaReadPlanEntry {
    fn new(
        name: String,
        length: u32,
    ) -> Self {
        Self {
            name,
            length,
            read: 0,
        }
    }

    fn go_to(
        &mut self,
        position: u32,
    ) -> anyhow::Result<()> {
        if position >= self.read {
            anyhow::bail!("Trying to set position, which has already been read")
        } else if position > self.length {
            anyhow::bail!("Trying to set position out of bounds");
        }
        self.read = position;
        Ok(())
    }
}

pub(crate) struct FastaReadPlan {
    plan: VecDeque<FastaReadPlanEntry>,
}

impl FastaReadPlan {
    pub(crate) fn new(index: FaIndex) -> Self {
        let mut plan = VecDeque::new();
        for entry in index.as_ref().to_vec().into_iter() {
            plan.push_back(FastaReadPlanEntry::new(
                String::from_utf8(entry.name().to_vec())
                    .expect("Failed to read string"),
                entry.length() as u32,
            ));
        }

        Self { plan }
    }

    pub(crate) fn go_to<R: AsRef<str>>(
        &mut self,
        chr: R,
        position: u32,
    ) -> anyhow::Result<()> {
        let current = self
            .plan
            .front_mut()
            .ok_or(anyhow::anyhow!("Plan is empty"))?;
        if chr.as_ref() != current.name.as_str() {
            anyhow::bail!(
                "Chromosome name differs from expected. Check if input data is sorted as Fasta"
            )
        };
        current.go_to(position)
    }

    pub(crate) fn pop_front(&mut self) -> Option<FastaReadPlanEntry> {
        self.plan.pop_front()
    }
}
