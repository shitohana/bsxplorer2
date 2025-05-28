use crate::data_structs::batch::BsxBatch;
use crate::data_structs::MethylationStats;

#[allow(unused)]
struct GenomeWideStats {
    stat: MethylationStats,
}

impl FromIterator<BsxBatch> for GenomeWideStats {
    fn from_iter<I: IntoIterator<Item = BsxBatch>>(iter: I) -> Self {
        let mut stat = MethylationStats::new();
        for batch in iter {
            stat.merge(&batch.get_methylation_stats());
        }
        stat.finalize_methylation();
        GenomeWideStats { stat }
    }
}
