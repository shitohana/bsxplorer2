use std::fs::File;
use std::io::{
    BufReader,
    Read,
};
use std::path::PathBuf;

use bsxplorer2::data_structs::coords::ContigIntervalMap;
use bsxplorer2::prelude::*;
use bsxplorer2::tools::dimred::merge::{
    merge_breakpoints,
    EqFloat,
    MergeType,
};
use bsxplorer2::utils::THREAD_POOL;
use clap::Args;
use rayon::prelude::*;

use crate::dimred::write_imap;
use crate::utils::{expand_wildcard, CliError};
use crate::PipelineCommand;

pub fn read_segments<R: Read>(handle: R) -> anyhow::Result<ContigIntervalMap<EqFloat>> {
    bio::io::bed::Reader::new(handle)
        .records()
        .into_iter()
        .map(|record_res| -> anyhow::Result<(Contig, EqFloat)> {
            let record = record_res?;
            let score = record
                .score()
                .ok_or(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "File should contain a valid score column",
                ))?
                .to_string()
                .parse::<EqFloat>()?;
            let contig = Contig::from(record);
            Ok((contig, score))
        })
        .collect::<anyhow::Result<ContigIntervalMap<EqFloat>>>()
}

#[derive(Debug, Clone, Args)]
pub struct MergeArgs {
    #[arg(required = true)]
    paths:      Vec<String>,
    #[arg(required = true, default_value_t = MergeType::Full)]
    merge_type: MergeType,
    #[arg(required = true, help = "Prefix for output files")]
    output:     String,
}

impl PipelineCommand for MergeArgs {
    fn run(&self) -> anyhow::Result<()> {
        let paths = self
            .paths
            .iter()
            .map(|s| expand_wildcard(s))
            .collect::<Result<Vec<_>, CliError>>()?
            .concat();

        let intervals = THREAD_POOL.install(
            || -> anyhow::Result<Vec<ContigIntervalMap<EqFloat>>> {
                paths
                    .iter()
                    .map(|p| {
                        File::open(p)
                            .map(BufReader::new)
                            .map_err(|e| CliError::CantOpenFile(p.to_owned(), e))
                    })
                    .collect::<Result<Vec<_>, CliError>>()?
                    .into_par_iter()
                    .map(|f| read_segments(f))
                    .collect::<anyhow::Result<Vec<_>>>()
            },
        )?;

        let merged = merge_breakpoints(&intervals, self.merge_type.clone());

        paths.par_iter().zip(merged.par_iter())
            .map(|(path, imap): (&PathBuf, &ContigIntervalMap<_>)| -> anyhow::Result<()> {
                write_imap(File::create(path)?, imap)
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        todo!()
    }
}
