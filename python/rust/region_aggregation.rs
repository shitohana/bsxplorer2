use std::fs::File;
use std::path::PathBuf;

use bsxplorer2::data_structs::Strand;
use bsxplorer2::io::bsx::{BsxFileReader, RegionReader};
use bsxplorer2::tools::region_aggregation::{
    aggregate_counts_for_regions, parse_optional_context, AggregationOptions,
    RegionCountResult, RegionSpec, StrandPolicy,
};
use pyo3::exceptions::{PyFileNotFoundError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyfunction]
#[pyo3(signature = (
    methylation_path,
    regions,
    sample_id=None,
    context=None,
    strand_policy="both",
    min_total=0,
    chunk_size=10000,
    include_empty_regions=true
))]
pub fn aggregate_region_counts_rust(
    py: Python<'_>,
    methylation_path: PathBuf,
    regions: &Bound<'_, PyAny>,
    sample_id: Option<String>,
    context: Option<String>,
    strand_policy: &str,
    min_total: u32,
    chunk_size: usize,
    include_empty_regions: bool,
) -> PyResult<Vec<PyObject>> {
    if !methylation_path.exists() {
        return Err(PyFileNotFoundError::new_err(format!(
            "methylation_path does not exist: {}",
            methylation_path.display()
        )));
    }

    let region_specs = parse_region_specs(regions)?;
    let options = AggregationOptions {
        context: parse_optional_context(context.as_deref())
            .map_err(|err| PyValueError::new_err(err.to_string()))?,
        strand_policy: StrandPolicy::parse(strand_policy)
            .map_err(|err| PyValueError::new_err(err.to_string()))?,
        min_total,
        chunk_size,
        include_empty_regions,
    };

    let file = File::open(&methylation_path).map_err(|err| {
        PyFileNotFoundError::new_err(format!(
            "Could not open methylation_path {}: {err}",
            methylation_path.display()
        ))
    })?;
    let reader = BsxFileReader::try_new(file)
        .map_err(|err| PyRuntimeError::new_err(err.to_string()))?;
    let mut region_reader = RegionReader::from_reader(reader)
        .map_err(|err| PyRuntimeError::new_err(err.to_string()))?;

    let results = aggregate_counts_for_regions(
        &mut region_reader,
        &region_specs,
        sample_id.as_deref(),
        &options,
    )
    .map_err(|err| PyRuntimeError::new_err(err.to_string()))?;

    results
        .iter()
        .map(|result| result_to_dict(py, result))
        .collect()
}

fn parse_region_specs(regions: &Bound<'_, PyAny>) -> PyResult<Vec<RegionSpec>> {
    let iterator = regions.iter().map_err(|_| {
        PyValueError::new_err("regions must be an iterable of dict-like rows")
    })?;
    let mut region_specs = Vec::new();

    for (idx, item) in iterator.enumerate() {
        let item = item?;
        let dict = item
            .downcast::<PyDict>()
            .map_err(|_| PyValueError::new_err("each region row must be a dict"))?;
        region_specs.push(parse_region_spec(dict, idx)?);
    }

    Ok(region_specs)
}

fn parse_region_spec(
    dict: &Bound<'_, PyDict>,
    idx: usize,
) -> PyResult<RegionSpec> {
    let region_id = get_string(dict, &["region_id", "name", "id"])?
        .unwrap_or_else(|| format!("region_{}", idx + 1));
    let seqname = get_string(dict, &["seqname", "chrom", "chr", "chromosome"])?
        .ok_or_else(|| PyValueError::new_err("region row is missing seqname/chrom"))?;
    let start = get_u32(dict, &["start"])?
        .ok_or_else(|| PyValueError::new_err("region row is missing start"))?;
    let end = get_u32(dict, &["end"])?
        .ok_or_else(|| PyValueError::new_err("region row is missing end"))?;

    if start > end {
        return Err(PyValueError::new_err(format!(
            "invalid coordinates for {region_id}: start {start} > end {end}"
        )));
    }

    let strand = get_string(dict, &["strand"])?
        .as_deref()
        .map(parse_optional_strand)
        .transpose()?;
    let context = parse_optional_context(get_string(dict, &["context"])?.as_deref())
        .map_err(|err| PyValueError::new_err(err.to_string()))?;

    Ok(RegionSpec {
        region_id,
        seqname,
        start,
        end,
        strand,
        context,
    })
}

fn get_string(
    dict: &Bound<'_, PyDict>,
    keys: &[&str],
) -> PyResult<Option<String>> {
    for key in keys {
        if let Some(value) = dict.get_item(*key)? {
            if value.is_none() {
                continue;
            }
            return value.extract::<String>().map(Some);
        }
    }
    Ok(None)
}

fn get_u32(
    dict: &Bound<'_, PyDict>,
    keys: &[&str],
) -> PyResult<Option<u32>> {
    for key in keys {
        if let Some(value) = dict.get_item(*key)? {
            if value.is_none() {
                continue;
            }
            let raw_value = value.extract::<u64>()?;
            if raw_value > u32::MAX as u64 {
                return Err(PyValueError::new_err(format!(
                    "{key} is too large for BSX2 coordinates: {raw_value}"
                )));
            }
            return Ok(Some(raw_value as u32));
        }
    }
    Ok(None)
}

fn parse_optional_strand(value: &str) -> PyResult<Strand> {
    match value {
        "+" | "plus" | "Plus" => Ok(Strand::Forward),
        "-" | "minus" | "Minus" => Ok(Strand::Reverse),
        "." | "" | "none" | "None" => Ok(Strand::None),
        other => Err(PyValueError::new_err(format!(
            "Unsupported strand value: {other}"
        ))),
    }
}

fn result_to_dict(
    py: Python<'_>,
    result: &RegionCountResult,
) -> PyResult<PyObject> {
    let row = PyDict::new_bound(py);
    row.set_item("region_id", &result.region_id)?;
    row.set_item("seqname", &result.seqname)?;
    row.set_item("chrom", &result.seqname)?;
    row.set_item("start", result.start)?;
    row.set_item("end", result.end)?;
    row.set_item("strand", &result.strand)?;
    row.set_item("context", &result.context)?;
    row.set_item("sample_id", &result.sample_id)?;
    row.set_item("mC", result.m_c)?;
    row.set_item("uC", result.u_c)?;
    row.set_item("total", result.total)?;
    row.set_item("n_cytosines", result.n_cytosines)?;
    match result.mean_methylation {
        Some(value) => row.set_item("mean_methylation", value)?,
        None => row.set_item("mean_methylation", py.None())?,
    };
    row.set_item("coverage_qc", &result.coverage_qc)?;
    Ok(row.into_py(py))
}
