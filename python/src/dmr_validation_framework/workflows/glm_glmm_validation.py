#!/usr/bin/env python
"""Run the controlled GLM-like vs CpG-level GLMM confirmation layer."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from dmr_validation_framework.core.columns import find_col, normalize_region_table
from dmr_validation_framework.core.io import read_table
from dmr_validation_framework.models import cpg_level_glmm, glm_glmm
from dmr_validation_framework.reports import glmm_plots
from dmr_validation_framework.reports.caller_summary import discover_path


DEFAULT_CONTEXTS = "CG,CHG,CHH"
DEFAULT_GLM_CALLER = "methylKit"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a small confirmatory GLM-like vs CpG-level GLMM layer from existing caller outputs. "
            "By default this uses methylKit DMRs as the aggregated GLM-like candidate set and validates "
            "top regions with glmmTMB on CpG-level count data."
        )
    )
    parser.add_argument("--root", required=True, help="Dataset root with metadata/sample_manifest.tsv.")
    parser.add_argument("--external-root", help="Directory with caller outputs. Defaults to <root>/external_callers.")
    parser.add_argument("--out-dir", help="Defaults to <root>/glm_glmm_validation.")
    parser.add_argument("--rscript", help="Optional explicit Rscript/Rscript.exe path for glmmTMB.")
    parser.add_argument("--contexts", default=DEFAULT_CONTEXTS)
    parser.add_argument("--glm-caller", default=DEFAULT_GLM_CALLER)
    parser.add_argument("--top-n-per-context", type=int, default=10)
    parser.add_argument("--q-threshold", type=float, default=0.10)
    parser.add_argument(
        "--candidate-q-threshold",
        type=float,
        help="Optional q-value filter before top-N selection. Defaults to --q-threshold when q-values exist.",
    )
    parser.add_argument("--condition-column", default="condition")
    parser.add_argument("--case-label", default="treatment")
    parser.add_argument("--control-label", default="control")
    parser.add_argument("--covariates", default="")
    parser.add_argument("--min-cpg", type=int, default=3)
    parser.add_argument("--min-replicates-per-group", type=int, default=2)
    parser.add_argument("--min-total", type=int, default=5)
    parser.add_argument("--max-zero-coverage-fraction", type=float, default=0.8)
    parser.add_argument("--coverage-set-mode", choices=("per_sample", "common"), default="per_sample")
    parser.add_argument("--min-common-cpgs", type=int, default=3)
    parser.add_argument("--no-coverage-qc", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def split_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def sha256_file(path: Path) -> str | None:
    if not path.exists() or path.is_dir():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_sample_manifest(root: Path) -> pd.DataFrame:
    path = root / "metadata" / "sample_manifest.tsv"
    if not path.exists():
        raise FileNotFoundError(f"sample manifest not found: {path}")
    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = [str(col).strip().lstrip("\ufeff").strip('"') for col in df.columns]
    missing = {"sample_id", "condition"} - set(df.columns)
    if missing:
        raise ValueError(f"sample_manifest.tsv missing required columns: {', '.join(sorted(missing))}")
    df["sample_id"] = df["sample_id"].astype(str)
    df["condition"] = df["condition"].astype(str)
    return df


def caller_slug(caller: str) -> str:
    lower = caller.lower().replace("_", "").replace("-", "")
    if lower == "methylkit":
        return "methylkit"
    if lower == "bsmooth":
        return "bsmooth"
    if lower == "combp":
        return "combp"
    return lower


def select_caller_candidates(
    external_root: Path,
    *,
    caller: str,
    contexts: list[str],
    top_n_per_context: int,
    q_threshold: float,
    candidate_q_threshold: float | None,
) -> pd.DataFrame:
    selected: list[pd.DataFrame] = []
    threshold = q_threshold if candidate_q_threshold is None else candidate_q_threshold
    for context in contexts:
        path = discover_path(external_root, caller, context)
        if path is None or not path.exists():
            raise FileNotFoundError(f"{caller} {context} caller table not found under {external_root}")
        raw = read_table(path)
        regions, notes = normalize_region_table(raw, source=caller)
        if regions.empty:
            raise ValueError(f"{caller} {context} table has no usable chrom/start/end intervals: {path}; notes={notes}")
        regions["context"] = context
        regions["source_caller"] = caller
        if "region_id" not in regions.columns:
            regions["region_id"] = [f"{caller}_{context}_{idx + 1}" for idx in range(len(regions))]
        regions["region_id"] = regions["region_id"].astype(str)
        if "q_value" in regions.columns and regions["q_value"].notna().any():
            filtered = regions[pd.to_numeric(regions["q_value"], errors="coerce") <= threshold].copy()
            if not filtered.empty:
                regions = filtered
        if "delta" in regions.columns:
            regions["_abs_delta_sort"] = pd.to_numeric(regions["delta"], errors="coerce").abs().fillna(-1.0)
        else:
            regions["_abs_delta_sort"] = -1.0
        if "q_value" in regions.columns:
            regions["_q_sort"] = pd.to_numeric(regions["q_value"], errors="coerce").fillna(float("inf"))
        else:
            regions["_q_sort"] = float("inf")
        regions = regions.sort_values(["_q_sort", "_abs_delta_sort"], ascending=[True, False])
        selected.append(regions.head(top_n_per_context).copy())
    if not selected:
        return pd.DataFrame()
    out = pd.concat(selected, ignore_index=True)
    keep = [
        col
        for col in (
            "region_id",
            "chrom",
            "start",
            "end",
            "context",
            "delta",
            "q_value",
            "source_caller",
            "caller",
        )
        if col in out.columns
    ]
    out = out[keep].copy()
    if "caller" in out.columns and "source_caller" not in out.columns:
        out = out.rename(columns={"caller": "source_caller"})
    if "source_caller" not in out.columns:
        out["source_caller"] = caller
    return out.drop_duplicates("region_id", keep="first").reset_index(drop=True)


def write_glm_like_results(regions: pd.DataFrame, out_path: Path) -> None:
    q_col = find_col(regions, ["q_value", "qvalue", "fdr", "padj"])
    delta_col = find_col(regions, ["delta", "delta_methylation"])
    out = pd.DataFrame(
        {
            "region_id": regions["region_id"].astype(str),
            "context": regions["context"].astype(str) if "context" in regions else "NA",
            "glm_status": "ok",
        }
    )
    out["glm_p_value"] = pd.NA
    out["glm_q_value"] = pd.to_numeric(regions[q_col], errors="coerce") if q_col else pd.NA
    out["glm_delta"] = pd.to_numeric(regions[delta_col], errors="coerce") if delta_col else pd.NA
    out.to_csv(out_path, sep="\t", index=False)


def methylkit_input_path(external_root: Path, context: str, sample_id: str) -> Path:
    path = external_root / "methylkit_input" / context / f"{sample_id}.bismark_coverage.tsv"
    if path.exists():
        return path
    matches = sorted((external_root / "methylkit_input" / context).glob(f"{sample_id}*.bismark_coverage.tsv"))
    if matches:
        return matches[0]
    raise FileNotFoundError(f"methylKit count input not found for sample={sample_id}, context={context}: {path}")


def extract_counts_from_methylkit_inputs(
    external_root: Path,
    regions: pd.DataFrame,
    manifest: pd.DataFrame,
    *,
    min_total: int,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    required = {"region_id", "chrom", "start", "end", "context"}
    missing = required - set(regions.columns)
    if missing:
        raise ValueError(f"selected regions missing required columns for CpG extraction: {', '.join(sorted(missing))}")
    region_table = regions.copy()
    region_table["start"] = pd.to_numeric(region_table["start"], errors="coerce").astype("Int64")
    region_table["end"] = pd.to_numeric(region_table["end"], errors="coerce").astype("Int64")
    region_table = region_table.dropna(subset=["start", "end"]).copy()
    for context, context_regions in region_table.groupby("context", sort=False):
        context = str(context)
        for sample in manifest["sample_id"].astype(str).tolist():
            path = methylkit_input_path(external_root, context, sample)
            counts = pd.read_csv(
                path,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "percent_methylation", "mC", "uC"],
            )
            counts["chrom"] = counts["chrom"].astype(str)
            counts["start"] = pd.to_numeric(counts["start"], errors="coerce")
            counts["mC"] = pd.to_numeric(counts["mC"], errors="coerce").fillna(0).astype(int)
            counts["uC"] = pd.to_numeric(counts["uC"], errors="coerce").fillna(0).astype(int)
            counts["total"] = counts["mC"] + counts["uC"]
            counts = counts[counts["total"] >= min_total].copy()
            if counts.empty:
                continue
            counts_by_chrom = {chrom: df for chrom, df in counts.groupby("chrom", sort=False)}
            for region in context_regions.itertuples(index=False):
                region_counts = counts_by_chrom.get(str(region.chrom))
                if region_counts is None or region_counts.empty:
                    continue
                mask = (region_counts["start"] >= int(region.start)) & (region_counts["start"] <= int(region.end))
                hit = region_counts.loc[mask, ["chrom", "start", "mC", "uC", "total"]].copy()
                if hit.empty:
                    continue
                hit["region_id"] = str(region.region_id)
                hit["context"] = context
                hit["sample_id"] = sample
                hit["pos"] = hit["start"].astype(int)
                hit["cpg_id"] = hit["chrom"].astype(str) + ":" + hit["pos"].astype(str)
                rows.append(hit[["region_id", "context", "chrom", "pos", "cpg_id", "sample_id", "mC", "uC", "total"]])
    if not rows:
        return pd.DataFrame(columns=["region_id", "context", "chrom", "pos", "cpg_id", "sample_id", "mC", "uC", "total"])
    return pd.concat(rows, ignore_index=True)


def write_counts_summary(counts: pd.DataFrame, out_path: Path) -> None:
    if counts.empty:
        summary = pd.DataFrame(
            [{"n_rows": 0, "n_regions": 0, "n_cpg": 0, "n_samples": 0, "mean_total": pd.NA}]
        )
    else:
        summary = pd.DataFrame(
            [
                {
                    "n_rows": len(counts),
                    "n_regions": counts["region_id"].nunique(),
                    "n_cpg": counts["cpg_id"].nunique(),
                    "n_samples": counts["sample_id"].nunique(),
                    "mean_total": pd.to_numeric(counts["total"], errors="coerce").mean(),
                }
            ]
        )
    summary.to_csv(out_path, sep="\t", index=False)


def write_manifest(out_dir: Path, args: argparse.Namespace, outputs: dict[str, str], status: str) -> None:
    payload = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "root": args.root,
        "external_root": args.external_root,
        "glm_caller": args.glm_caller,
        "contexts": split_csv(args.contexts),
        "parameters": {
            "top_n_per_context": args.top_n_per_context,
            "q_threshold": args.q_threshold,
            "candidate_q_threshold": args.candidate_q_threshold,
            "min_cpg": args.min_cpg,
            "min_total": args.min_total,
            "coverage_set_mode": args.coverage_set_mode,
        },
        "outputs": outputs,
        "output_hashes": {
            name: sha256_file(Path(path)) for name, path in outputs.items() if Path(path).exists() and Path(path).is_file()
        },
        "note": (
            "This workflow is a controlled confirmatory layer. It does not run a genome-wide GLMM caller; "
            "it validates selected caller candidates with CpG-level count models."
        ),
    }
    (out_dir / "glm_glmm_workflow_manifest.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> int:
    root = Path(args.root)
    external_root = Path(args.external_root) if args.external_root else root / "external_callers"
    out_dir = Path(args.out_dir) if args.out_dir else root / "glm_glmm_validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    contexts = split_csv(args.contexts)
    slug = caller_slug(args.glm_caller)

    manifest = read_sample_manifest(root)
    selected_path = out_dir / f"selected_{slug}_top{args.top_n_per_context}_per_context_regions.tsv"
    glm_path = out_dir / f"{slug}_as_glm_results.tsv"
    counts_path = out_dir / "region_cpg_counts.tsv"
    counts_summary_path = out_dir / "region_cpg_counts_region_summary.tsv"
    design_path = out_dir / "sample_design.tsv"
    glmm_dir = out_dir / "glmm"
    comparison_dir = out_dir / "glm_vs_glmm"
    plots_dir = out_dir / "glmm_plots"

    selected = select_caller_candidates(
        external_root,
        caller=args.glm_caller,
        contexts=contexts,
        top_n_per_context=args.top_n_per_context,
        q_threshold=args.q_threshold,
        candidate_q_threshold=args.candidate_q_threshold,
    )
    if selected.empty:
        raise ValueError("no GLM-like candidate regions selected")
    selected.to_csv(selected_path, sep="\t", index=False)
    write_glm_like_results(selected, glm_path)

    design = manifest.copy()
    design.to_csv(design_path, sep="\t", index=False)

    counts = extract_counts_from_methylkit_inputs(external_root, selected, manifest, min_total=args.min_total)
    counts.to_csv(counts_path, sep="\t", index=False)
    write_counts_summary(counts, counts_summary_path)

    if counts.empty:
        raise ValueError("no CpG-level counts overlap the selected regions")

    code = cpg_level_glmm.run(
        argparse.Namespace(
            region_cpg_counts=str(counts_path),
            design=str(design_path),
            rscript=args.rscript,
            dmr_evidence=str(selected_path),
            top_n=len(selected),
            condition_column=args.condition_column,
            case_label=args.case_label,
            control_label=args.control_label,
            covariates=args.covariates,
            min_cpg=args.min_cpg,
            min_replicates_per_group=args.min_replicates_per_group,
            min_total=args.min_total,
            max_zero_coverage_fraction=args.max_zero_coverage_fraction,
            qvalue_method="BH",
            coverage_set_mode=args.coverage_set_mode,
            min_coverage=args.min_total,
            min_covered_per_group=None,
            min_covered_A=None,
            min_covered_B=None,
            min_common_cpgs=args.min_common_cpgs,
            write_coverage_qc=not args.no_coverage_qc,
            out_dir=str(glmm_dir),
        )
    )
    if code != 0:
        return code

    code = glm_glmm.run(
        argparse.Namespace(
            glm_results=str(glm_path),
            glmm_results=str(glmm_dir / "cpg_level_glmm_results.tsv"),
            out_dir=str(comparison_dir),
            q_threshold=args.q_threshold,
        )
    )
    if code != 0:
        return code

    if not args.no_plots:
        code = glmm_plots.run(
            argparse.Namespace(
                comparison=str(comparison_dir / "glm_vs_glmm_comparison.tsv"),
                glmm_results=str(glmm_dir / "cpg_level_glmm_results.tsv"),
                region_cpg_counts=str(counts_path),
                out_dir=str(plots_dir),
                design=str(design_path),
                condition_column=args.condition_column,
                case_label=args.case_label,
                control_label=args.control_label,
                example_region_id=None,
                region_delta=None,
                q_threshold=args.q_threshold,
            )
        )
        if code != 0:
            return code

    outputs = {
        "selected_regions": str(selected_path),
        "glm_results": str(glm_path),
        "region_cpg_counts": str(counts_path),
        "region_cpg_counts_summary": str(counts_summary_path),
        "design": str(design_path),
        "glmm_results": str(glmm_dir / "cpg_level_glmm_results.tsv"),
        "glm_vs_glmm_comparison": str(comparison_dir / "glm_vs_glmm_comparison.tsv"),
        "glmm_plots": str(plots_dir),
    }
    write_manifest(out_dir, args, outputs, "completed")
    print(f"glm/glmm manifest: {out_dir / 'glm_glmm_workflow_manifest.json'}")
    print(f"glm/glmm comparison: {comparison_dir / 'glm_vs_glmm_comparison.tsv'}")
    if not args.no_plots:
        print(f"glm/glmm plots: {plots_dir}")
    return 0


def main(argv: list[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except Exception as exc:  # noqa: BLE001
        print(f"glm-glmm-validation failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
