use std::hash::Hash;
use std::ops::Range;

use bio::data_structures::interval_tree::IntervalTree;
use hashbrown::HashMap;
use indexmap::IndexSet;
use itertools::Itertools;
use num::{PrimInt, Unsigned};
use serde::{Deserialize, Serialize};

use crate::data_structs::coords::{Contig, GenomicPosition};

/// Index for batches in a BSX file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchIndex<S, P>
where
    S: AsRef<str> + Clone + Hash + Eq,
    P: Unsigned + PrimInt, {
    map:       HashMap<S, IntervalTree<GenomicPosition<S, P>, usize>>,
    chr_order: IndexSet<S>,
}

impl<S, P> Default for BatchIndex<S, P>
where
    S: AsRef<str> + Clone + Hash + Eq,
    P: Unsigned + PrimInt,
{
    fn default() -> Self { Self::new() }
}

impl<S, P> BatchIndex<S, P>
where
    S: AsRef<str> + Clone + Eq + Hash,
    P: Unsigned + PrimInt,
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
        self.chr_order
            .insert(contig.seqname().clone());

        self.map
            .entry(contig.seqname())
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
                    self.chr_order
                        .get_index_of(&contig.seqname())
                        .unwrap_or(0),
                    contig,
                )
            })
            .sorted_by(|(left_chr, left_contig), (right_chr, right_contig)| {
                left_chr.cmp(right_chr).then(
                    left_contig
                        .start()
                        .cmp(&right_contig.start()),
                )
            })
            .map(|(_, contig)| contig)
    }

    /// Find the batch indices that overlap with a given contig.
    pub fn find(
        &self,
        contig: &Contig<S, P>,
    ) -> Option<Vec<usize>> {
        if let Some(tree) = self.map.get(&contig.seqname()) {
            let batches = tree
                .find(Range::<_>::from(contig.clone()))
                .map(|entry| *entry.data())
                .collect_vec();
            Some(batches)
        }
        else {
            None
        }
    }

    /// Returns the chromosome order.
    pub fn chr_order(&self) -> &IndexSet<S> { &self.chr_order }

    /// Returns the underlying map.
    pub fn map(
        &self
    ) -> &HashMap<S, IntervalTree<GenomicPosition<S, P>, usize>> {
        &self.map
    }
}
