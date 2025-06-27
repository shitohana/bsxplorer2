mod merge;
mod shrink;

pub(crate) use merge::MergeArgs;
pub(crate) use shrink::DimRedArgs;

pub(self) fn write_imap<W: std::io::Write>(
    writer: W,
    values: &bsxplorer2::prelude::ContigIntervalMap<
        bsxplorer2::tools::dimred::merge::EqFloat,
    >,
) -> anyhow::Result<()> {
    let file = std::io::BufWriter::new(writer);
    let mut writer = bio::io::bed::Writer::new(file);

    use itertools::Itertools;
    for chr in values.chr_names().iter().sorted() {
        let chr_imap = values.inner().get(chr).unwrap();
        for entry in chr_imap.intervals.iter() {
            let val: bsxplorer2::tools::dimred::merge::EqFloat = entry.val;
            let mut record = bio::io::bed::Record::new();
            record.set_start(entry.start as u64);
            record.set_end(entry.stop as u64);
            record.set_chrom(chr.as_str());
            record.set_score(val.to_string().as_str());

            writer.write(&record)?;
        }
    }
    Ok(())
}
