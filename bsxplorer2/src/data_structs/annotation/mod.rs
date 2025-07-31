/// Contains data structures and logic for handling genomic annotations.
///
/// This module provides the core types for representing and storing genomic
/// features, primarily based on the GFF3 specification but also compatible
/// with BED.
///
/// The main components are:
/// - [`GffEntry`]: Represents a single feature/annotation entry from a file
///   like GFF or BED. It holds coordinate information, feature type, and
///   attributes.
/// - [`GffEntryAttributes`]: A structured representation of the key-value pairs
///   found in the GFF attributes column.
/// - [`RawGffEntry`]: An intermediate structure used for deserializing raw
///   lines from a GFF file.
/// - [`HcAnnotStore`]: A high-level data structure that stores multiple
///   [`GffEntry`] instances. It provides efficient ways to query entries by
///   genomic location using an interval tree and by parent-child relationships
///   using a tree structure.
///
/// The module also includes functionality for parsing GFF and BED files
/// into the [`HcAnnotStore`].
mod annot_store;
mod gff_entry;

pub use annot_store::{
    EntryId,
    EntryTree,
    HcAnnotStore
};
pub use gff_entry::{
    GffEntry,
    GffEntryAttributes,
    RawGffEntry,
};

#[cfg(test)]
mod tests;
