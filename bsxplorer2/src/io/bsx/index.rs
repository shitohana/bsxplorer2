use std::io::{
    Read,
    Write,
};
use indexmap::IndexSet;
use itertools::Itertools;
use polars::error::PolarsResult;
use serde::{
    Deserialize,
    Serialize,
};

use super::BsxFileReader;
use crate::data_structs::coords::{
    Contig, ContigIntervalMap
};
use crate::data_structs::typedef::BsxSmallStr;

/// Index for batches in a BSX file.
#[derive(Debug, Clone, Serialize, Deserialize)]
/// Index for batches in a BSX file.
///
/// This index maps genomic regions (represented as `Contig` objects) to the
/// batch indices within the file where data for those regions can be found.
/// It also maintains an ordered list of chromosome names.
pub struct BatchIndex {
    map:       ContigIntervalMap<usize>,
    chr_order: IndexSet<BsxSmallStr>,
}

impl FromIterator<(Contig, usize)> for BatchIndex {
    fn from_iter<I: IntoIterator<Item = (Contig, usize)>>(iter: I) -> Self {
        let mut index = Self::new();
        for (contig, batch_idx) in iter {
            index.insert(contig, batch_idx);
        }
        index
    }
}

impl BatchIndex {
    /// Serializes the `BatchIndex` to a writer using bincode.
    ///
    /// # Arguments
    ///
    /// * `writer`: The writer to serialize the index to.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on success, or a `bincode::error::EncodeError` on
    /// failure.
    pub fn to_file<W: Write>(
        self,
        writer: &mut W,
    ) -> Result<(), bincode::error::EncodeError> {
        let config = bincode::config::standard();
        bincode::serde::encode_into_std_write(self, writer, config)?;
        Ok(())
    }

    /// Deserializes a `BatchIndex` from a reader using bincode.
    ///
    /// # Arguments
    ///
    /// * `reader`: The reader to deserialize the index from.
    ///
    /// # Returns
    ///
    /// Returns the deserialized `BatchIndex` on success, or a
    /// `bincode::error::DecodeError` on failure.
    pub fn from_file<R: Read>(
        reader: &mut R
    ) -> Result<Self, bincode::error::DecodeError> {
        let config = bincode::config::standard();
        bincode::serde::decode_from_std_read(reader, config)
    }

    /// Builds a `BatchIndex` by iterating through batches in a `BsxFileReader`.
    ///
    /// This method reads the contig information from each batch header and
    /// builds the interval tree and chromosome order.
    ///
    /// # Arguments
    ///
    /// * `reader`: The `BsxFileReader` to build the index from.
    ///
    /// # Returns
    ///
    /// Returns the built `BatchIndex` on success, or a `PolarsResult` error on
    /// failure.
    pub fn from_reader(reader: &mut BsxFileReader) -> PolarsResult<Self> {
        let contigs = reader
            .iter()
            .enumerate()
            .map(|(batch_idx, batch)| {
                batch.map(|b| (b.as_contig().unwrap(), batch_idx))
            })
            .collect::<PolarsResult<Vec<_>>>()?;
        Ok(Self::from_iter(contigs))
    }
}

impl Default for BatchIndex {
    fn default() -> Self {
        Self::new()
    }
}

impl BatchIndex {
    /// Creates a new empty `BatchIndex`.
    pub fn new() -> Self {
        Self {
            map:       ContigIntervalMap::new(),
            chr_order: IndexSet::new(),
        }
    }

    /// Inserts a contig and its corresponding batch index into the index.
    ///
    /// The contig's chromosome name is added to the chromosome order if not
    /// already present. The contig's interval is inserted into the interval
    /// tree for its chromosome.
    ///
    /// # Arguments
    ///
    /// * `contig`: The genomic contig to insert.
    /// * `batch_idx`: The index of the batch containing this contig.
    pub fn insert(
        &mut self,
        contig: Contig,
        batch_idx: usize,
    ) {
        self.chr_order.insert(contig.seqname().clone());

        self.map.insert(contig, batch_idx);
    }

    /// Sorts a set of contigs according to the internal chromosome order and
    /// start position.
    ///
    /// Contigs on chromosomes not present in the `chr_order` will be placed
    /// first in an arbitrary order relative to each other, then sorted by
    /// start.
    ///
    /// # Arguments
    ///
    /// * `contigs`: An iterator over the contigs to sort.
    ///
    /// # Returns
    ///
    /// An iterator yielding the sorted contigs.
    pub fn sort<I>(
        &self,
        contigs: I,
    ) -> impl Iterator<Item = Contig>
    where
        I: IntoIterator<Item = Contig>, {
        contigs
            .into_iter()
            .map(|contig| {
                (
                    // Use 0 for chromosomes not found, effectively putting them first
                    self.chr_order.get_index_of(contig.seqname()).unwrap_or(0),
                    contig,
                )
            })
            .sorted_by(|(left_chr, left_contig), (right_chr, right_contig)| {
                left_chr
                    .cmp(right_chr)
                    .then(left_contig.start().cmp(&right_contig.start()))
            })
            .map(|(_, contig)| contig)
    }

    /// Finds the batch indices that overlap with a given contig (query region).
    ///
    /// Searches the interval tree for the chromosome matching the query
    /// contig's sequence name to find overlapping intervals and returns
    /// their associated batch indices.
    ///
    /// # Arguments
    ///
    /// * `contig`: The query contig representing the genomic region to search.
    ///
    /// # Returns
    ///
    /// An `Option` containing a `Vec` of batch indices if overlapping intervals
    /// are found for the specified chromosome, otherwise `None`.
    pub fn find(
        &self,
        contig: &Contig,
    ) -> Option<Vec<usize>> {
        self.map.find(contig).map(|v| v.into_iter().cloned().collect())
    }

    /// Returns the chromosome order stored in the index.
    ///
    /// This order reflects the sequence in which chromosomes were encountered
    /// or explicitly inserted.
    pub fn get_chr_order(&self) -> &IndexSet<BsxSmallStr> {
        &self.chr_order
    }

    /// Gets the index of a chromosome name within the stored chromosome order.
    ///
    /// # Arguments
    ///
    /// * `chr`: The chromosome name to look up.
    ///
    /// # Returns
    ///
    /// An `Option` containing the zero-based index of the chromosome name
    /// in the order, or `None` if the chromosome name is not found.
    pub fn get_chr_index(
        &self,
        chr: &BsxSmallStr,
    ) -> Option<usize> {
        self.chr_order.get_index_of(chr)
    }

    /// Returns a reference to the underlying map storing the interval trees
    /// per chromosome.
    pub fn map(&self) -> &ContigIntervalMap<usize> {
        &self.map
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_structs::Strand;

    #[test]
    fn test_insert_and_find() {
        let mut index: BatchIndex = BatchIndex::new();

        let contig1 = Contig::new("chr1".into(), 1, 100, Strand::None);
        let contig2 = Contig::new("chr1".into(), 50, 150, Strand::None);
        let contig3 = Contig::new("chr2".into(), 1, 100, Strand::None);

        index.insert(contig1.clone(), 1);
        index.insert(contig2.clone(), 2);
        index.insert(contig3.clone(), 3);

        // Test finding overlapping contigs
        let query_contig1 = Contig::new("chr1".into(), 20, 60, Strand::None);
        let result1 = index.find(&query_contig1).unwrap();
        assert_eq!(result1, vec![1, 2]);

        // Test finding non-overlapping contigs
        let query_contig2 = Contig::new("chr1".into(), 200, 300, Strand::None);
        let result2 = index.find(&query_contig2);
        assert_eq!(result2, None);

        // Test finding contigs on different chromosomes
        let query_contig3 = Contig::new("chr2".into(), 50, 60, Strand::None);
        let result3 = index.find(&query_contig3).unwrap();
        assert_eq!(result3, vec![3]);
    }

    #[test]
    fn test_sort() {
        let mut index: BatchIndex = BatchIndex::new();
        index.chr_order.insert("chr2".into());
        index.chr_order.insert("chr1".into());

        let contig1 = Contig::new("chr1".into(), 50, 150, Strand::None);
        let contig2 = Contig::new("chr1".into(), 1, 100, Strand::None);
        let contig3 = Contig::new("chr2".into(), 1, 100, Strand::None);
        let contig4 = Contig::new("chr2".into(), 50, 150, Strand::None);

        let contigs = vec![
            contig1.clone(),
            contig2.clone(),
            contig3.clone(),
            contig4.clone(),
        ];

        let sorted_contigs: Vec<_> = index.sort(contigs).collect();

        assert_eq!(sorted_contigs[0], contig3);
        assert_eq!(sorted_contigs[1], contig4);
        assert_eq!(sorted_contigs[2], contig2);
        assert_eq!(sorted_contigs[3], contig1);
    }
}
