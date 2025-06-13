mod frombsx;
mod r2r;
mod tobsx;

use std::path::PathBuf;

use bsxplorer2::prelude::*;
use clap::Args;
pub use frombsx::FromBsxConvert;
pub use r2r::R2RConvert;
pub use tobsx::ToBsxConvert;

#[derive(Debug, Clone, Args)]
pub(self) struct FromReportArgs {
    #[clap(short='f', long = "from", required = true, value_enum, default_value_t = ReportType::Bismark)]
    from_type:        ReportType,
    #[clap(short='F', long = "from-compression", required = true, value_enum, default_value_t = Compression::None)]
    from_compression: Compression,
    #[arg(
        long,
        default_value_t = false,
        help = "Use less RAM, but elongate computation."
    )]
    low_memory:       bool,
    #[arg(long = "fa", help = "Path to the reference sequence file.")]
    fasta_path:       Option<PathBuf>,
    #[arg(long, default_value_t = 2 << 20, help = "Size of raw batches.")]
    batch_size:       usize,
}
