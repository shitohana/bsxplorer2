from __future__ import annotations

import subprocess
import sys
import os
from pathlib import Path

import pandas as pd
import pytest

from bsx2.analysis.region_signal import (
    RegionSignalConfig,
    aggregate_region_signal,
    available_region_signal_backends,
)
from bsx2.analysis.region_signal_rust import rust_region_aggregator_available


REPO_ROOT = Path(__file__).resolve().parents[3]
SUBPROCESS_ENV = {**os.environ, "PYTHONPATH": str(REPO_ROOT / "python" / "src")}


def _regions() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "region_id": ["r1", "empty"],
            "chrom": ["chr1", "chr1"],
            "start": [1, 100],
            "end": [30, 120],
            "context": ["CG", "CG"],
            "strand": ["+", "+"],
        }
    )


def _counts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s1"],
            "chrom": ["chr1", "chr1", "chr1"],
            "position": [10, 20, 30],
            "strand": ["+", "-", "+"],
            "context": ["CG", "CHG", "CG"],
            "mC": [4, 2, 7],
            "uC": [6, 3, 3],
        }
    )


def test_backend_discovery_reports_pandas() -> None:
    backends = available_region_signal_backends()
    assert backends["pandas"] is True
    assert backends["default"] in {"rust", "pandas"}


def test_pandas_fallback_still_works() -> None:
    out = aggregate_region_signal(
        _regions(),
        _counts(),
        RegionSignalConfig(context="CG", backend="pandas"),
    )
    first = out[out["region_id"] == "r1"].iloc[0]
    assert first["mC"] == 11
    assert first["uC"] == 9
    assert first["total"] == 20
    assert first["coverage_qc"] == "ok"


def test_rust_backend_unavailable_is_clear() -> None:
    if rust_region_aggregator_available():
        pytest.skip("Rust binding is available in this environment")
    with pytest.raises(RuntimeError, match="Rust binding"):
        aggregate_region_signal(
            _regions(),
            None,
            RegionSignalConfig(backend="rust", methylation_path="missing.bsx"),
        )


def test_include_empty_regions_pandas() -> None:
    out = aggregate_region_signal(
        _regions(),
        _counts(),
        RegionSignalConfig(context="CG", backend="pandas", include_empty_regions=True),
    )
    empty = out[out["region_id"] == "empty"].iloc[0]
    assert empty["coverage_qc"] == "no_records"


def test_seqname_alias_map_pandas() -> None:
    regions = pd.DataFrame({"region_id": ["r1"], "chrom": ["chrA01"], "start": [1], "end": [3]})
    counts = pd.DataFrame({"sample_id": ["s1"], "chrom": ["A01"], "position": [2], "mC": [3], "uC": [1]})
    out = aggregate_region_signal(
        regions,
        counts,
        RegionSignalConfig(seqname_aliases={"chrA01": "A01"}, backend="pandas"),
    )
    assert out["mC"].iloc[0] == 3


def test_invalid_coordinates_raise_clear_error() -> None:
    regions = pd.DataFrame({"region_id": ["bad"], "chrom": ["chr1"], "start": [30], "end": [1]})
    with pytest.raises(ValueError):
        aggregate_region_signal(regions, _counts(), RegionSignalConfig(backend="pandas"))


def test_auto_backend_with_opposite_uses_pandas() -> None:
    out = aggregate_region_signal(
        _regions(),
        _counts(),
        RegionSignalConfig(backend="auto", strand_policy="opposite", methylation_path="would_prefer_rust.bsx"),
    )
    first = out[out["region_id"] == "r1"].iloc[0]
    assert first["mC"] == 2
    assert first["uC"] == 3


def test_rust_backend_with_opposite_raises_clear_error() -> None:
    with pytest.raises(ValueError, match="strand_policy='opposite' is not supported by Rust backend"):
        aggregate_region_signal(
            _regions(),
            _counts(),
            RegionSignalConfig(backend="rust", strand_policy="opposite", methylation_path="sample.bsx"),
        )


def test_cli_backend_pandas(tmp_path: Path) -> None:
    regions_path = tmp_path / "regions.tsv"
    counts_path = tmp_path / "counts.tsv"
    out_path = tmp_path / "out.tsv"
    qc_path = tmp_path / "qc.tsv"
    summary_path = tmp_path / "summary.md"
    _regions().to_csv(regions_path, sep="\t", index=False)
    _counts().to_csv(counts_path, sep="\t", index=False)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "dmr_validation_framework.workflows.region_signal_aggregation",
            "--regions",
            str(regions_path),
            "--counts",
            str(counts_path),
            "--out",
            str(out_path),
            "--qc-out",
            str(qc_path),
            "--summary-out",
            str(summary_path),
            "--backend",
            "pandas",
        ],
        cwd=REPO_ROOT,
        env=SUBPROCESS_ENV,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out_path.exists()
    assert "backend_used: pandas" in summary_path.read_text(encoding="utf-8")


def test_cli_backend_auto(tmp_path: Path) -> None:
    regions_path = tmp_path / "regions.tsv"
    counts_path = tmp_path / "counts.tsv"
    out_path = tmp_path / "out.tsv"
    qc_path = tmp_path / "qc.tsv"
    _regions().to_csv(regions_path, sep="\t", index=False)
    _counts().to_csv(counts_path, sep="\t", index=False)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "dmr_validation_framework.workflows.region_signal_aggregation",
            "--regions",
            str(regions_path),
            "--counts",
            str(counts_path),
            "--out",
            str(out_path),
            "--qc-out",
            str(qc_path),
            "--backend",
            "auto",
        ],
        cwd=REPO_ROOT,
        env=SUBPROCESS_ENV,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out_path.exists()


def test_cli_backend_rust_clear_message_when_unavailable(tmp_path: Path) -> None:
    if rust_region_aggregator_available():
        pytest.skip("Compiled Rust binding is available")
    regions_path = tmp_path / "regions.tsv"
    out_path = tmp_path / "out.tsv"
    qc_path = tmp_path / "qc.tsv"
    _regions().to_csv(regions_path, sep="\t", index=False)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "dmr_validation_framework.workflows.region_signal_aggregation",
            "--regions",
            str(regions_path),
            "--counts",
            str(tmp_path / "missing.bsx"),
            "--out",
            str(out_path),
            "--qc-out",
            str(qc_path),
            "--backend",
            "rust",
        ],
        cwd=REPO_ROOT,
        env=SUBPROCESS_ENV,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "Rust binding" in result.stderr or "methylation_path" in result.stderr
