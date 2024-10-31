use polars::prelude::SerWriter;
use polars_core::prelude::*;
use std::fs::File;
use std::io::{BufReader, Read, Seek};

/// Struct for sequence data
pub struct SeqData {
    /// Names of the records from FASTA
    pub chr_names: Vec<String>,
    /// Binary offsets of chroms
    pub chr_offsets: Vec<u64>,
    /// Number of G/C-nucleotides
    pub gc_content: u64,
    /// Number of non-G/C nucleotides
    pub non_gc_content: u64,
}

impl SeqData {
    fn new() -> Self {
        SeqData {
            chr_names: Vec::new(),
            chr_offsets: Vec::new(),
            gc_content: 0,
            non_gc_content: 0,
        }
    }
}

/// Scans fasta to determine
///
/// - Number of nucleotides
/// - G/C content count
/// - Chromosome names
/// - Sequence binary offsets
pub fn scan_fasta(path: &str) -> SeqData {
    let mut reader: BufReader<File> = BufReader::new(File::open(path).unwrap());
    let mut buf = [0_u8; 1];

    let mut reading_name = false;
    let mut reading_tag = false;
    let mut name = String::new();
    let mut seq_offset: u64;

    let mut seq_data = SeqData::new();

    while reader.read_exact(&mut buf).is_ok() {
        let c = u8::from_le_bytes(buf);

        // Start of tag
        if c == 0x3e {
            reading_name = true
        };
        // End of name
        if reading_name && c == 0x20 {
            reading_tag = true;
            reading_name = false;
            seq_data.chr_names.push(name.clone());
            name.clear()
        }
        // End of tag
        if reading_tag && c == 0x0a {
            reading_tag = false;
            // Save seq offset
            seq_offset = reader
                .stream_position()
                .expect("Could not get current position!");
            seq_data.chr_offsets.push(seq_offset);
        }

        // Push name char (if name continues)
        if reading_name && !reading_tag {
            if c != 0x3e {
                name.push(c as char)
            }
        }
        // Add GC content
        else if !reading_name && !reading_tag {
            // If not newline
            match c {
                // T|t|G|g
                0x54 | 0x74 | 0x47 | 0x67 => seq_data.gc_content += 1,
                // \n
                0x0a => {}
                // Other
                _ => seq_data.non_gc_content += 1,
            }
        }
    }
    seq_data
}
