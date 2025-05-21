mod annot_store;
mod gff_entry;

pub use annot_store::AnnotStore;
pub use gff_entry::{GffEntry, GffEntryAttributes, RawGffEntry};

#[cfg(test)]
mod tests;
