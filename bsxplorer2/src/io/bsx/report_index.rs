use std::{fs::File, io, path::Path, sync::Arc};
use memmap2::Mmap;
use memchr::memchr;

/// An index entry mapping a genomic chunk to file offsets
#[derive(Debug, Clone)]
pub struct IndexEntry {
    /// Chromosome name (interned)
    pub chrom: Arc<String>,
    /// Start position of chunk (1-based)
    pub start: u32,
    /// End position of chunk (inclusive)
    pub end: u32,
    /// Byte offset in file where chunk begins
    pub offset_start: u64,
    /// Byte offset in file where chunk ends (exclusive)
    pub offset_end: u64,
}

/// Builds an index for a methylation report ASCII CSV file with optional header and skipped rows.
/// Uses memory-mapping and zero-copy parsing for speed.
///
/// # Arguments
/// * `path` - input file path
/// * `target_chunk_size` - byte-threshold at which to split chunks
/// * `has_header` - whether the first line is a header (skip it)
/// * `skip_rows` - number of additional lines to skip after the header
pub fn build_index<P: AsRef<Path>>(
    path: P,
    target_chunk_size: usize,
    has_header: bool,
    skip_rows: usize,
) -> io::Result<Vec<IndexEntry>> {
    // Memory-map file
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let data = &mmap[..];
    let mut offset: u64 = 0;
    let total_len = data.len() as u64;

    // Skip header line if requested
    if has_header {
        if let Some(pos) = memchr(b'\n', &data[offset as usize..]) {
            offset += (pos + 1) as u64;
        } else {
            return Ok(Vec::new());
        }
    }
    // Skip additional rows
    for _ in 0..skip_rows {
        if let Some(pos) = memchr(b'\n', &data[offset as usize..]) {
            offset += (pos + 1) as u64;
        } else {
            break;
        }
    }

    let mut entries = Vec::new();
    let mut chunk_start_offset = offset;
    let mut chunk_size = 0;
    let mut current_chrom: Arc<String> = Arc::new(String::new());
    let mut chunk_start_pos = 0;
    let mut previous_pos = 0;
    let mut first_record = true;

    // Iterate each data line
    while offset < total_len {
        let line_start = offset as usize;
        let rem = &data[line_start..];
        // Find newline
        let nl = match memchr(b'\n', rem) {
            Some(n) => n,
            None => break,
        };
        let line_end = line_start + nl;
        let line = &data[line_start..line_end];
        // Advance offset past this line
        offset = (line_end + 1) as u64;
        let line_len = nl + 1;

        // Parse chromosome and position fields
        let c1 = memchr(b',', line)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing comma"))?;
        let c2 = memchr(b',', &line[c1 + 1..])
            .map(|p| p + c1 + 1)
            .unwrap_or(c1 + 1);
        let chr_bytes = &line[..c1];
        let pos = fast_parse_u32(&line[c1 + 1..c2]);

        // On first record, initialize chunk & chromosome
        if first_record {
            let chrom_str = unsafe { std::str::from_utf8_unchecked(chr_bytes) };
            current_chrom = Arc::new(chrom_str.to_string());
            chunk_start_pos = pos;
            first_record = false;
        } else if chr_bytes != current_chrom.as_bytes() {
            // Chromosome changed: flush previous chunk
            entries.push(IndexEntry {
                chrom: Arc::clone(&current_chrom),
                start: chunk_start_pos,
                end: previous_pos,
                offset_start: chunk_start_offset,
                offset_end: line_start as u64,
            });
            // Start new chunk
            let chrom_str = unsafe { std::str::from_utf8_unchecked(chr_bytes) };
            current_chrom = Arc::new(chrom_str.to_string());
            chunk_start_pos = pos;
            chunk_start_offset = line_start as u64;
            chunk_size = 0;
        }

        // Accumulate size and split on threshold if more lines follow
        chunk_size += line_len;
        if chunk_size >= target_chunk_size && offset < total_len {
            entries.push(IndexEntry {
                chrom: Arc::clone(&current_chrom),
                start: chunk_start_pos,
                end: pos,
                offset_start: chunk_start_offset,
                offset_end: offset,
            });
            // Begin next chunk
            chunk_start_offset = offset;
            chunk_start_pos = pos;
            chunk_size = 0;
        }

        previous_pos = pos;
    }

    // Flush last chunk if any records processed
    if !first_record {
        entries.push(IndexEntry {
            chrom: Arc::clone(&current_chrom),
            start: chunk_start_pos,
            end: previous_pos,
            offset_start: chunk_start_offset,
            offset_end: offset,
        });
    }

    Ok(entries)
}

/// Parses ASCII-digit slice into u32 without checks
#[inline]
fn fast_parse_u32(data: &[u8]) -> u32 {
    let mut acc = 0;
    for &b in data {
        acc = acc * 10 + (b - b'0') as u32;
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Helper to write a test file with given lines
    fn write_file(lines: &[&str]) -> NamedTempFile {
        let mut tmp = NamedTempFile::new().unwrap();
        for &line in lines {
            write!(tmp, "{}\n", line).unwrap();
        }
        tmp
    }

    #[test]
    fn test_empty_file() {
        let tmp = NamedTempFile::new().unwrap();
        assert!(build_index(tmp.path(), 100, true, 0).unwrap().is_empty());
        assert!(build_index(tmp.path(), 100, false, 0).unwrap().is_empty());
    }

    #[test]
    fn test_single_record_with_and_without_header() {
        // With header
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "chrom,pos,other").unwrap();
        writeln!(tmp, "chrX,42,foo").unwrap();
        let idx = build_index(tmp.path(), 10, true, 0).unwrap();
        assert_eq!(idx.len(), 1);
        assert_eq!(&*idx[0].chrom, "chrX");
        assert_eq!(idx[0].start, 42);
        assert_eq!(idx[0].end, 42);
        // Without header
        let mut tmp2 = NamedTempFile::new().unwrap();
        writeln!(tmp2, "chrY,100,bar").unwrap();
        let idx2 = build_index(tmp2.path(), 10, false, 0).unwrap();
        assert_eq!(idx2.len(), 1);
        assert_eq!(&*idx2[0].chrom, "chrY");
    }

    #[test]
    fn test_skip_rows() {
        let lines = &["chrom,pos,other", "skip,1,a", "skip,2,b", "chrZ,5,x"];
        let tmp = write_file(lines);
        let idx = build_index(tmp.path(), 10, false, 2).unwrap();
        assert_eq!(idx.len(), 1);
        assert_eq!(&*idx[0].chrom, "chrZ");
        assert_eq!(idx[0].start, 5);
    }

    #[test]
    fn test_multiple_chunks_and_chromosomes() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "chrom,pos,other").unwrap();
        for i in 1..=10 { writeln!(tmp, "chr1,{},x", i).unwrap(); }
        for i in 1..=5  { writeln!(tmp, "chr2,{},y", i).unwrap(); }
        let idx = build_index(tmp.path(), 10, true, 0).unwrap();
        assert!(idx.iter().any(|e| &*e.chrom == "chr1"));
        assert!(idx.iter().any(|e| &*e.chrom == "chr2"));
    }

    #[test]
    fn test_exact_threshold_boundary() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "chrom,pos,other").unwrap();
        let r1 = "chrA,1,a";
        let r2 = "chrA,2,bb";
        writeln!(tmp, "{}", r1).unwrap();
        writeln!(tmp, "{}", r2).unwrap();
        let rec_bytes = r1.len() + 1 + r2.len() + 1;
        let idx = build_index(tmp.path(), rec_bytes, true, 0).unwrap();
        assert_eq!(idx.len(), 1);
    }

    #[test]
    fn test_large_dataset_accuracy() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "chrom,pos,other").unwrap();
        let chrs = ["chr1", "chr2", "chr3"];
        let count = 1_000;
        for &chr in &chrs {
            for i in 1..=count {
                writeln!(tmp, "{},{}", chr, i).unwrap();
            }
        }
        let idx = build_index(tmp.path(), 2048, true, 0).unwrap();

        // Map chromosome to total covered positions
        let mut coverage: HashMap<String, usize> = HashMap::new();
        for entry in idx {
            let len = (entry.end - entry.start + 1) as usize;
            *coverage.entry((*entry.chrom).clone()).or_default() += len;
        }
        for &chr in &chrs {
            assert_eq!(coverage.get(chr), Some(&count));
        }
    }
}
