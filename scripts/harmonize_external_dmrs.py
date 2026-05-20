#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "src"))

from bsx2.analysis import adapter_for_caller, build_caller_support_matrix, write_schema_validation_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Import external DMR candidates into the canonical BSX2 schema.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--caller", required=True, choices=["generic_bed", "DSS", "methylKit", "dmrseq", "metilene"])
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--contrast-id", default="")
    parser.add_argument("--condition-a", default="")
    parser.add_argument("--condition-b", default="")
    parser.add_argument("--context")
    parser.add_argument("--counts-manifest")
    parser.add_argument("--dmr-region-count-tests")
    parser.add_argument("--dmr-evidence-scores")
    parser.add_argument("--annotation")
    parser.add_argument("--enable-evidence-join", default="false")
    parser.add_argument("--enable-beta-binomial-join", default="false")
    args = parser.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adapter_cls = adapter_for_caller(args.caller)
    adapter = adapter_cls(args.input, contrast_id=args.contrast_id, condition_a=args.condition_a, condition_b=args.condition_b, context=args.context)
    canonical = adapter.read()
    canonical.to_csv(out / "external_dmr_candidates_canonical.tsv", sep="\t", index=False)
    write_schema_validation_report(canonical, out / "external_dmr_schema_validation.tsv")
    pd.DataFrame({"warning": adapter.warnings or [""]}).to_csv(out / "external_dmr_import_warnings.tsv", sep="\t", index=False)
    canonical.to_csv(out / "dmr_evidence_harmonized.tsv", sep="\t", index=False)
    summary = pd.DataFrame([{"source_caller": args.caller, "n_candidates": len(canonical)}])
    summary.to_csv(out / "dmr_caller_comparison_summary.tsv", sep="\t", index=False)
    build_caller_support_matrix([canonical]).to_csv(out / "dmr_caller_support_matrix.tsv", sep="\t", index=False)
    (out / "dmr_external_method_audit.md").write_text("# External DMR Method Audit\n\nSchema harmonization only; external callers were not run.\n", encoding="utf-8")
    (out / "external_dmr_harmonization_manifest.json").write_text(json.dumps({"created_at": datetime.now(timezone.utc).isoformat(), "input": args.input, "caller": args.caller, "outputs": sorted(p.name for p in out.iterdir())}, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
