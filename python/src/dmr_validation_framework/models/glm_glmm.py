#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


from dmr_validation_framework.core.glm_glmm import (
    compare_beta_binomial_glm_vs_glmm,
    read_cpg_level_glmm_results,
    read_region_level_glm_results,
)


def sha256_file(path: Path) -> str | None:
    if not path.exists() or path.is_dir():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare aggregated beta-binomial GLM results with confirmatory CpG-level GLMM results."
    )
    parser.add_argument("--glm-results", required=True)
    parser.add_argument("--glmm-results", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--q-threshold", type=float, default=0.05)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    glm = read_region_level_glm_results(args.glm_results)
    glmm = read_cpg_level_glmm_results(args.glmm_results)
    comparison = compare_beta_binomial_glm_vs_glmm(glm, glmm, q_threshold=args.q_threshold)

    out_path = out_dir / "glm_vs_glmm_comparison.tsv"
    summary_path = out_dir / "glm_vs_glmm_comparison_summary.md"
    manifest_path = out_dir / "glm_vs_glmm_comparison_manifest.json"
    comparison.to_csv(out_path, sep="\t", index=False)
    summary_path.write_text(
        "\n".join([
            "# GLM vs CpG-level GLMM Comparison Summary",
            "",
            "The CpG-level GLMM is used as a confirmatory layer for selected DMR candidates, not as a genome-wide DMR caller.",
            "",
            f"- n_regions: {len(comparison)}",
            f"- confirmed_by_glmm: {int(comparison['confirmed_by_glmm'].sum())}",
            f"- glm_only_candidate: {int(comparison['glm_only_candidate'].sum())}",
            f"- glmm_more_conservative: {int(comparison['glmm_more_conservative'].sum())}",
            "",
        ]),
        encoding="utf-8",
    )
    manifest = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "glm_results": str(args.glm_results),
            "glmm_results": str(args.glmm_results),
        },
        "parameters": {"q_threshold": args.q_threshold},
        "outputs": {
            out_path.name: {"path": str(out_path), "sha256": sha256_file(out_path)},
            summary_path.name: {"path": str(summary_path), "sha256": sha256_file(summary_path)},
        },
        "note": "Comparison only; no DMR statistical model was changed.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0



def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
