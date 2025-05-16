use std::io::{Read, Write};
use std::ops::Range;

use bio::data_structures::interval_tree::IntervalTree;
use hashbrown::HashMap;
use indexmap::IndexSet;
use itertools::Itertools;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};

use crate::data_structs::coords::{Contig, GenomicPosition};
use crate::data_structs::typedef::{SeqNameStr, SeqPosNum};

/// Index for batches in a BSX file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchIndex<S, P>
where
    S: SeqNameStr,
    P: SeqPosNum, {
    map:       HashMap<S, IntervalTree<GenomicPosition<S, P>, usize>>,
    chr_order: IndexSet<S>,
}

impl<S, P> BatchIndex<S, P>
where
    S: SeqNameStr + Serialize,
    P: SeqPosNum + Serialize,
{
    pub fn to_file<W: Write>(
        self,
        writer: &mut W,
    ) -> Result<(), bincode::error::EncodeError> {
        let config = bincode::config::standard();
        bincode::serde::encode_into_std_write(self, writer, config)?;
        Ok(())
    }
}

impl<S, P> BatchIndex<S, P>
where
    S: SeqNameStr + for<'de> Deserialize<'de> + DeserializeOwned,
    P: SeqPosNum + for<'de> Deserialize<'de> + DeserializeOwned,
{
    pub fn from_file<R: Read>(
        reader: &mut R
    ) -> Result<Self, bincode::error::DecodeError> {
        let config = bincode::config::standard();
        bincode::serde::decode_from_std_read(reader, config)
    }
}

impl<S, P> Default for BatchIndex<S, P>
where
    S: SeqNameStr,
    P: SeqPosNum,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<S, P> BatchIndex<S, P>
where
    S: SeqNameStr,
    P: SeqPosNum,
{
    pub fn new() -> Self {
        Self {
            map:       HashMap::new(),
            chr_order: IndexSet::new(),
        }
    }

    /// Insert a contig and its corresponding batch index.
    pub fn insert(
        &mut self,
        contig: Contig<S, P>,
        batch_idx: usize,
    ) {
        self.chr_order.insert(contig.seqname().clone());

        self.map
            .entry(contig.seqname().to_owned())
            .and_modify(|tree| {
                tree.insert(Range::<_>::from(contig.clone()), batch_idx);
            })
            .or_insert_with(|| {
                let mut tree = IntervalTree::new();
                tree.insert(Range::<_>::from(contig), batch_idx);
                tree
            });
    }

    /// Sort a set of contigs according to the chromosome order and start
    /// position.
    pub fn sort<I>(
        &self,
        contigs: I,
    ) -> impl Iterator<Item = Contig<S, P>>
    where
        I: IntoIterator<Item = Contig<S, P>>, {
        contigs
            .into_iter()
            .map(|contig| {
                (
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

    /// Find the batch indices that overlap with a given contig.
    pub fn find(
        &self,
        contig: &Contig<S, P>,
    ) -> Option<Vec<usize>> {
        if let Some(tree) = self.map.get(contig.seqname()) {
            let batches = tree
                .find(Range::<_>::from(contig.clone()))
                .map(|entry| *entry.data())
                .collect_vec();
            if batches.is_empty() {
                None
            }
            else {
                Some(batches)
            }
        }
        else {
            None
        }
    }

    /// Returns the chromosome order.
    pub fn get_chr_order(&self) -> &IndexSet<S> {
        &self.chr_order
    }

    pub fn get_chr_index(
        &self,
        chr: &S,
    ) -> Option<usize> {
        self.chr_order.get_index_of(chr)
    }

    /// Returns the underlying map.
    pub fn map(&self) -> &HashMap<S, IntervalTree<GenomicPosition<S, P>, usize>> {
        &self.map
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_structs::enums::Strand;

    #[test]
    fn test_insert_and_find() {
        let mut index: BatchIndex<String, u64> = BatchIndex::new();

        let contig1 = Contig::new("chr1".to_string(), 1, 100, Strand::None);
        let contig2 = Contig::new("chr1".to_string(), 50, 150, Strand::None);
        let contig3 = Contig::new("chr2".to_string(), 1, 100, Strand::None);

        index.insert(contig1.clone(), 1);
        index.insert(contig2.clone(), 2);
        index.insert(contig3.clone(), 3);

        // Test finding overlapping contigs
        let query_contig1 = Contig::new("chr1".to_string(), 20, 60, Strand::None);
        let result1 = index.find(&query_contig1).unwrap();
        assert_eq!(result1, vec![1, 2]);

        // Test finding non-overlapping contigs
        let query_contig2 = Contig::new("chr1".to_string(), 200, 300, Strand::None);
        let result2 = index.find(&query_contig2);
        assert_eq!(result2, None);

        // Test finding contigs on different chromosomes
        let query_contig3 = Contig::new("chr2".to_string(), 50, 60, Strand::None);
        let result3 = index.find(&query_contig3).unwrap();
        assert_eq!(result3, vec![3]);
    }

    #[test]
    fn test_sort() {
        let mut index: BatchIndex<String, u64> = BatchIndex::new();
        index.chr_order.insert("chr2".to_string());
        index.chr_order.insert("chr1".to_string());

        let contig1 = Contig::new("chr1".to_string(), 50, 150, Strand::None);
        let contig2 = Contig::new("chr1".to_string(), 1, 100, Strand::None);
        let contig3 = Contig::new("chr2".to_string(), 1, 100, Strand::None);
        let contig4 = Contig::new("chr2".to_string(), 50, 150, Strand::None);

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
