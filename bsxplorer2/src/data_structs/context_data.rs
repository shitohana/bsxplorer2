use crate::data_structs::region::GenomicPosition;
use crate::utils::types::PosNum;
use itertools::Itertools;

use polars::df;
use polars::error::PolarsResult;
use polars::frame::DataFrame;

/// Structure to hold context data_structs for positions
pub struct ContextData<N>
where
    N: PosNum, {
    /// Genomic positions
    positions: Vec<N>,
    /// Methylation contexts (CpG=Some(true), CHG=Some(false), CHH=None)
    contexts:  Vec<Option<bool>>,
    /// Strand information (true=forward, false=reverse)
    strands:   Vec<bool>,
}

impl<N: PosNum> ContextData<N> {
    /// Get the number of positions in this context data_structs
    pub fn len(&self) -> usize { self.contexts.len() }

    /// Check if the context data_structs is empty
    pub fn is_empty(&self) -> bool { self.contexts.is_empty() }

    /// Create a new empty ContextData
    pub(crate) fn new() -> Self {
        ContextData {
            positions: Vec::new(),
            contexts:  Vec::new(),
            strands:   Vec::new(),
        }
    }

    /// Filter the context data_structs based on a predicate function on
    /// positions
    pub(crate) fn filter<F: Fn(N) -> bool>(
        &self,
        predicate: F,
    ) -> Self {
        let mut new_self = Self::new();
        for (idx, pos) in self.positions.iter().enumerate() {
            if predicate(*pos) {
                new_self.add_row(*pos, self.contexts[idx], self.strands[idx])
            }
        }
        new_self
    }

    /// Get the column names for this context data_structs
    #[allow(dead_code)]
    pub(crate) fn col_names() -> &'static [&'static str] {
        &["position", "context", "strand"]
    }

    /// Get the name of the position column
    pub(crate) fn position_col() -> &'static str { "position" }

    /// Create a new ContextData with the given capacity
    pub(crate) fn with_capacity(capacity: usize) -> Self {
        ContextData {
            positions: Vec::with_capacity(capacity),
            contexts:  Vec::with_capacity(capacity),
            strands:   Vec::with_capacity(capacity),
        }
    }

    /// Add a row to the context data_structs
    pub(crate) fn add_row(
        &mut self,
        position: N,
        context: Option<bool>,
        strand: bool,
    ) {
        self.positions.push(position);
        self.strands.push(strand);
        self.contexts.push(context);
    }

    /// Create context data_structs from a DNA sequence starting at the given
    /// position
    pub fn from_sequence(
        seq: &[u8],
        start: GenomicPosition<N>,
    ) -> Self {
        let start_pos = start.position();
        let fw_bound: usize = seq.len() - 2;
        let rv_bound: usize = 2;

        let mut new = Self::with_capacity(seq.len());

        let ascii_seq = seq.to_ascii_uppercase();

        'seq_iter: for (index, nuc) in ascii_seq.iter().enumerate() {
            let forward = match nuc {
                b'C' => true,
                b'G' => false,
                _ => continue 'seq_iter,
            };
            let context = if forward {
                if index >= fw_bound {
                    continue 'seq_iter;
                };

                if ascii_seq[index + 1] == b'G' {
                    Some(true)
                }
                else if ascii_seq[index + 2] == b'G' {
                    Some(false)
                }
                else {
                    None
                }
            }
            else {
                if index <= rv_bound {
                    continue 'seq_iter;
                };

                if ascii_seq[index - 1] == b'C' {
                    Some(true)
                }
                else if ascii_seq[index - 2] == b'C' {
                    Some(false)
                }
                else {
                    None
                }
            };

            new.add_row(
                start_pos
                    + N::from(index).unwrap_or_else(|| {
                        panic!(
                            "Failed to convert index {} to position type",
                            index
                        )
                    }),
                context,
                forward,
            );
        }

        new.shrink_to_fit();
        new
    }

    /// Shrink the capacity of the internal vectors to fit their contents
    pub(crate) fn shrink_to_fit(&mut self) {
        self.positions.shrink_to_fit();
        self.contexts.shrink_to_fit();
        self.strands.shrink_to_fit();
    }

    /// Convert the context data_structs to a DataFrame
    pub fn into_dataframe(self) -> PolarsResult<DataFrame> {
        df![
            "position" => self.positions.iter().map(|x| x.to_u64().unwrap()).collect_vec(),
            "context" => self.contexts,
            "strand" => self.strands
        ]
    }

    /// Get a reference to the positions
    pub fn positions(&self) -> &Vec<N> { &self.positions }

    /// Get a reference to the contexts
    pub fn contexts(&self) -> &Vec<Option<bool>> { &self.contexts }

    /// Get a reference to the strands
    pub fn strands(&self) -> &Vec<bool> { &self.strands }
}
