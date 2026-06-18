from argparse import Namespace

import pandas as pd

from dmr_validation_framework.checks.overlap_thresholds import run
from dmr_validation_framework.core.io import load_caller_tables


def write_tsv(path, rows):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_load_caller_tables_prefers_native_over_support_matrix(tmp_path):
    write_tsv(
        tmp_path / "dmr_method_support_matrix.tsv",
        [
            {
                "region_id": f"r{i}",
                "chrom": "Chr1",
                "start": i * 100,
                "end": i * 100 + 100,
                "context": "CG",
                "dss_significant": "true",
                "methylkit_significant": "true",
                "dss_q_value": 0.001,
                "methylkit_q_value": 0.001,
            }
            for i in range(5)
        ],
    )
    write_tsv(
        tmp_path / "dss_results_canonical.tsv",
        [
            {"dmr_id": "dss_1", "chrom": "Chr1", "start": 0, "end": 100, "context": "CG", "q_value": 0.001},
            {"dmr_id": "dss_2", "chrom": "Chr1", "start": 200, "end": 300, "context": "CG", "q_value": 0.001},
        ],
    )
    write_tsv(
        tmp_path / "methylkit_results_canonical.tsv",
        [
            {"dmr_id": "mk_1", "chrom": "Chr1", "start": 20, "end": 80, "context": "CG", "q_value": 0.001},
            {"dmr_id": "mk_2", "chrom": "Chr1", "start": 210, "end": 260, "context": "CG", "q_value": 0.001},
        ],
    )

    tables = {table.name: table for table in load_caller_tables([tmp_path])}

    assert tables["DSS"].source_kind == "native"
    assert tables["methylKit"].source_kind == "native"
    assert len(tables["DSS"].data) == 2
    assert len(tables["methylKit"].data) == 2


def test_overlap_threshold_run_excludes_support_matrix(monkeypatch, tmp_path):
    write_tsv(
        tmp_path / "dmr_method_support_matrix.tsv",
        [
            {
                "region_id": f"support_{i}",
                "chrom": "Chr1",
                "start": i * 100,
                "end": i * 100 + 100,
                "context": "CG",
                "dss_significant": "true",
                "methylkit_significant": "true",
                "dss_q_value": 0.001,
                "methylkit_q_value": 0.001,
            }
            for i in range(10)
        ],
    )
    write_tsv(
        tmp_path / "dss_results_canonical.tsv",
        [
            {"dmr_id": "dss_1", "chrom": "Chr1", "start": 0, "end": 100, "context": "CG", "q_value": 0.001},
            {"dmr_id": "dss_2", "chrom": "Chr1", "start": 200, "end": 300, "context": "CG", "q_value": 0.001},
        ],
    )
    write_tsv(
        tmp_path / "methylkit_results_canonical.tsv",
        [
            {"dmr_id": "mk_1", "chrom": "Chr1", "start": 20, "end": 80, "context": "CG", "q_value": 0.001},
            {"dmr_id": "mk_2", "chrom": "Chr1", "start": 210, "end": 260, "context": "CG", "q_value": 0.001},
        ],
    )
    monkeypatch.setenv("DMR_VALIDATION_SEARCH_ROOTS", str(tmp_path))
    monkeypatch.chdir(tmp_path)
    out_dir = tmp_path / "validation"

    code = run(
        Namespace(
            out_dir=out_dir,
            thresholds="0.5,0.8",
            matching_policy="best_reciprocal",
        )
    )

    assert code == 0
    result = pd.read_csv(out_dir / "overlap_threshold_sensitivity.tsv", sep="\t")
    counts = dict(zip(result["threshold"], result["n_pairs"]))
    assert counts[0.5] == 2
    assert counts[0.8] == 0
    assert set(result["caller_a_source_kind"]) == {"native"}
    assert set(result["caller_b_source_kind"]) == {"native"}
