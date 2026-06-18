#!/usr/bin/env python
"""Import external DMR caller outputs and score support for target regions.

The workflow intentionally does not run DSS, methylKit, dmrseq, metilene, or
other third-party callers. It is the framework layer between caller-native
outputs and the existing DMR validation/bundle machinery:

1. discover or accept external caller output tables;
2. normalize them through the existing canonical adapters;
3. overlap canonical caller intervals with a target DMR-like table;
4. write a bundle-ready external caller support matrix.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from bsx2.analysis.dmr_harmonization import direction_from_delta, safe_float, support_column
from bsx2.analysis.intervals import reciprocal_overlap
from dmr_validation_framework.core.columns import first_existing
from dmr_validation_framework.core.harmonization import adapter_for_caller, canonical_dmr_columns
from dmr_validation_framework.core.io import read_table

DEFAULT_CONTEXTS = ("CG", "CHG", "CHH")

DEFAULT_CALLER_PATTERNS: dict[str, tuple[str, ...]] = {
    "DSS": ("dss/dss_dmrs_{context}.tsv", "DSS/DSS_dmrs_{context}.tsv"),
    "methylKit": ("methylkit/methylkit_dmrs_{context}.tsv", "methylKit/methylkit_dmrs_{context}.tsv"),
    "metilene": ("metilene/metilene_dmrs_{context}.tsv",),
    "dmrseq": ("dmrseq/dmrseq_dmrs_{context}.tsv",),
    "BSmooth": ("bsmooth/bsmooth_dmrs_{context}.tsv", "bsseq/bsmooth_dmrs_{context}.tsv"),
    "DMRcate": ("dmrcate/dmrcate_dmrs_{context}.tsv", "DMRcate/dmrcate_dmrs_{context}.tsv"),
    "comb-p": (
        "combp/combp_dmrs_{context}.tsv",
        "combp_like/combp_like_dmrs_{context}.tsv",
        "comb-p/combp_dmrs_{context}.tsv",
    ),
}

SUMMARY_COLUMNS = [
    "n_target_regions",
    "n_regions_with_external_support",
    "fraction_regions_with_external_support",
    "max_callers_supporting",
]

INVENTORY_COLUMNS = ["caller", "context", "path", "source", "exists", "status", "n_rows", "error"]


@dataclass(frozen=True)
class CallerOutputSpec:
    caller: str
    path: Path
    context: str | None = None
    source: str = "explicit"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Harmonize external DMR caller outputs and score target-region support.",
    )
    parser.add_argument("--root", help="Dataset root. Used to derive defaults.")
    parser.add_argument(
        "--target-dmr-table",
        help=(
            "Target/internal DMR-like table. Defaults to "
            "<root>/strict_all_windows/processed_geo_cx/run/dmr_evidence_scores.tsv."
        ),
    )
    parser.add_argument(
        "--external-root",
        help="Directory containing caller outputs. Defaults to <root>/external_callers.",
    )
    parser.add_argument(
        "--out-dir",
        help="Output directory. Defaults to <root>/strict_all_windows/external_callers_harmonized.",
    )
    parser.add_argument(
        "--caller-output",
        action="append",
        default=[],
        metavar="CALLER[:CONTEXT]=PATH",
        help=(
            "Explicit caller output table. Repeat as needed. Examples: "
            "DSS:CG=.../dss_dmrs_CG.tsv or metilene=.../metilene_all.tsv."
        ),
    )
    parser.add_argument(
        "--caller-output-manifest",
        help="JSON or TSV manifest with caller/path/context columns.",
    )
    parser.add_argument(
        "--context",
        action="append",
        help="Methylation context to discover. Can be repeated or comma-separated. Defaults to CG,CHG,CHH.",
    )
    parser.add_argument("--contrast-id", default="control_vs_treatment")
    parser.add_argument("--condition-a", default="control")
    parser.add_argument("--condition-b", default="treatment")
    parser.add_argument("--overlap-threshold", type=float, default=0.50)
    parser.add_argument(
        "--include-candidate-only",
        action="store_true",
        help="Include target rows with evidence_class=candidate_only.",
    )
    parser.add_argument(
        "--max-target-regions",
        type=int,
        default=200_000,
        help="Maximum target rows to score after sorting/filtering.",
    )
    parser.add_argument(
        "--no-auto-discover",
        action="store_true",
        help="Use only explicit --caller-output/manifest entries.",
    )
    parser.add_argument(
        "--write-run-plan",
        action="store_true",
        help="Write a caller-output run plan documenting expected files and prerequisites.",
    )
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def contexts_from_args(values: list[str] | None) -> list[str]:
    if not values:
        return list(DEFAULT_CONTEXTS)
    contexts: list[str] = []
    for value in values:
        for token in value.split(","):
            token = token.strip().upper()
            if token and token not in contexts:
                contexts.append(token)
    return contexts


def default_paths(args: argparse.Namespace) -> tuple[Path | None, Path | None, Path | None, Path | None]:
    root = Path(args.root) if args.root else None
    strict_root = root / "strict_all_windows" if root else None
    target = Path(args.target_dmr_table) if args.target_dmr_table else None
    external_root = Path(args.external_root) if args.external_root else None
    out_dir = Path(args.out_dir) if args.out_dir else None
    if root is not None:
        target = target or strict_root / "processed_geo_cx" / "run" / "dmr_evidence_scores.tsv"
        external_root = external_root or root / "external_callers"
        out_dir = out_dir or strict_root / "external_callers_harmonized"
    return root, target, external_root, out_dir


def parse_caller_output(value: str) -> CallerOutputSpec:
    if "=" not in value:
        raise ValueError(f"--caller-output must use CALLER[:CONTEXT]=PATH syntax: {value}")
    lhs, path = value.split("=", 1)
    if ":" in lhs:
        caller, context = lhs.split(":", 1)
        context = context.strip().upper() or None
    else:
        caller, context = lhs, None
    caller = caller.strip()
    if not caller:
        raise ValueError(f"Missing caller in --caller-output: {value}")
    return CallerOutputSpec(caller=caller, context=context, path=Path(path), source="explicit")


def read_manifest(path: str | Path) -> list[CallerOutputSpec]:
    path = Path(path)
    specs: list[CallerOutputSpec] = []
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        entries = payload.get("caller_outputs", payload if isinstance(payload, list) else [])
        for entry in entries:
            specs.append(
                CallerOutputSpec(
                    caller=str(entry["caller"]),
                    context=str(entry["context"]).upper() if entry.get("context") else None,
                    path=Path(entry["path"]),
                    source=f"manifest:{path}",
                )
            )
    else:
        df = read_table(path)
        missing = {"caller", "path"} - set(df.columns)
        if missing:
            raise ValueError(f"Caller manifest misses columns: {sorted(missing)}")
        for row in df.to_dict("records"):
            context = row.get("context")
            specs.append(
                CallerOutputSpec(
                    caller=str(row["caller"]),
                    context=str(context).upper() if context is not None and not pd.isna(context) else None,
                    path=Path(str(row["path"])),
                    source=f"manifest:{path}",
                )
            )
    return specs


def discover_caller_outputs(external_root: Path, contexts: list[str]) -> list[CallerOutputSpec]:
    specs: list[CallerOutputSpec] = []
    for caller, patterns in DEFAULT_CALLER_PATTERNS.items():
        for context in contexts:
            for pattern in patterns:
                path = external_root / pattern.format(context=context)
                if path.exists():
                    specs.append(
                        CallerOutputSpec(caller=caller, context=context, path=path, source="auto_discover")
                    )
                    break
    return specs


def collect_output_specs(args: argparse.Namespace, external_root: Path | None, contexts: list[str]) -> list[CallerOutputSpec]:
    specs = [parse_caller_output(value) for value in args.caller_output]
    if args.caller_output_manifest:
        specs.extend(read_manifest(args.caller_output_manifest))
    if external_root is not None and not args.no_auto_discover:
        specs.extend(discover_caller_outputs(external_root, contexts))

    deduped: list[CallerOutputSpec] = []
    seen: set[tuple[str, str | None, Path]] = set()
    for spec in specs:
        key = (spec.caller, spec.context, spec.path)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(spec)
    return deduped


def import_external_tables(
    specs: list[CallerOutputSpec],
    *,
    contrast_id: str,
    condition_a: str,
    condition_b: str,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    tables: list[pd.DataFrame] = []
    inventory: list[dict[str, Any]] = []
    for spec in specs:
        row = {
            "caller": spec.caller,
            "context": spec.context or "",
            "path": str(spec.path),
            "source": spec.source,
            "exists": spec.path.exists(),
            "status": "missing",
            "n_rows": 0,
            "error": "",
        }
        if not spec.path.exists():
            inventory.append(row)
            continue
        try:
            adapter = adapter_for_caller(spec.caller)(
                spec.path,
                contrast_id=contrast_id,
                condition_a=condition_a,
                condition_b=condition_b,
                context=spec.context,
            )
            table = adapter.read()
            table = table.dropna(subset=["chrom", "start", "end"]).copy()
            table["start"] = pd.to_numeric(table["start"], errors="coerce")
            table["end"] = pd.to_numeric(table["end"], errors="coerce")
            table = table.dropna(subset=["start", "end"])
            table = table[table["start"] < table["end"]].reset_index(drop=True)
            row["status"] = "ok"
            row["n_rows"] = len(table)
            if adapter.warnings:
                row["status"] = "warning"
                row["error"] = ";".join(adapter.warnings)
            tables.append(table)
        except Exception as exc:  # noqa: BLE001
            row["status"] = "error"
            row["error"] = str(exc)
        inventory.append(row)
    if not tables:
        return pd.DataFrame(), inventory
    return pd.concat(tables, ignore_index=True), inventory


def load_target_regions(
    path: str | Path,
    *,
    include_candidate_only: bool,
    max_target_regions: int,
) -> pd.DataFrame:
    target = read_table(path)
    required = {"chrom", "start", "end", "context"}
    missing = required - set(target.columns)
    if missing:
        raise ValueError(f"Target DMR table misses columns: {sorted(missing)}")
    if not include_candidate_only and "evidence_class" in target.columns:
        target = target[target["evidence_class"].astype(str).ne("candidate_only")].copy()
    sort_cols = [col for col in ("evidence_rank", "class_rank", "context_rank", "chrom", "start") if col in target.columns]
    if sort_cols:
        target = target.sort_values(sort_cols)
    if max_target_regions > 0 and len(target) > max_target_regions:
        target = target.head(max_target_regions).copy()

    target["start"] = pd.to_numeric(target["start"], errors="coerce")
    target["end"] = pd.to_numeric(target["end"], errors="coerce")
    target = target.dropna(subset=["chrom", "start", "end", "context"])
    target = target[target["start"] < target["end"]].copy()
    if "region_id" not in target.columns:
        target["region_id"] = [
            f"{row.chrom}:{int(row.start)}-{int(row.end)}:{row.context}"
            for row in target.itertuples(index=False)
        ]
    return target.reset_index(drop=True)


def value_column(caller: str, suffix: str) -> str:
    return support_column(caller).removesuffix("_support") + suffix


def build_external_index(external: pd.DataFrame) -> dict[tuple[str, str, str], pd.DataFrame]:
    groups: dict[tuple[str, str, str], pd.DataFrame] = {}
    if external.empty:
        return groups
    external = external.copy()
    external["source_caller"] = external["source_caller"].astype(str)
    external["chrom"] = external["chrom"].astype(str)
    external["context"] = external["context"].astype(str)
    for key, group in external.groupby(["source_caller", "chrom", "context"], dropna=False):
        groups[tuple(str(part) for part in key)] = group.sort_values("start").reset_index(drop=True)
    return groups


def matching_external_records(
    group: pd.DataFrame,
    *,
    start: float,
    end: float,
    overlap_threshold: float,
) -> list[tuple[float, pd.Series]]:
    if group.empty:
        return []
    candidates = group[(group["start"] <= end) & (group["end"] >= start)]
    matches: list[tuple[float, pd.Series]] = []
    for _, row in candidates.iterrows():
        overlap = reciprocal_overlap(start, end, float(row["start"]), float(row["end"]))
        if overlap >= overlap_threshold:
            matches.append((overlap, row))
    return matches


def build_target_support_matrix(
    target: pd.DataFrame,
    external: pd.DataFrame,
    *,
    overlap_threshold: float,
) -> pd.DataFrame:
    callers = sorted(external["source_caller"].dropna().astype(str).unique()) if not external.empty else []
    external_index = build_external_index(external)
    rows: list[dict[str, Any]] = []
    delta_col = first_existing(target.columns, ("delta", "delta_methylation", "meth_diff", "mean_methylation_difference"))
    q_col = first_existing(target.columns, ("q_value", "q_value_safe", "q", "fdr", "padj"))
    carry_cols = [
        col
        for col in (
            "evidence_class",
            "evidence_rank",
            "support_delta",
            "support_q",
            "support_coverage",
            "full_replicate_support",
        )
        if col in target.columns
    ]

    for target_row in target.to_dict("records"):
        start = float(target_row["start"])
        end = float(target_row["end"])
        chrom = str(target_row["chrom"])
        context = str(target_row["context"])
        support_callers: list[str] = []
        overlaps: list[float] = []
        q_by_caller: dict[str, float] = {}
        delta_by_caller: dict[str, float] = {}
        direction_by_caller: dict[str, str] = {}

        for caller in callers:
            group = external_index.get((caller, chrom, context))
            if group is None:
                continue
            matches = matching_external_records(
                group,
                start=start,
                end=end,
                overlap_threshold=overlap_threshold,
            )
            if not matches:
                continue
            support_callers.append(caller)
            overlaps.extend(match[0] for match in matches)
            q_values = [safe_float(row.get("q_value")) for _, row in matches]
            q_values = [value for value in q_values if value is not None]
            deltas = [safe_float(row.get("delta")) for _, row in matches]
            deltas = [value for value in deltas if value is not None]
            if q_values:
                q_by_caller[caller] = min(q_values)
            if deltas:
                delta_by_caller[caller] = max(deltas, key=abs)
                direction_by_caller[caller] = direction_from_delta(delta_by_caller[caller])
            else:
                direction_by_caller[caller] = "unknown"

        directions = [direction for direction in direction_by_caller.values() if direction != "unknown"]
        direction_consensus = "unknown"
        if directions:
            direction_consensus = "same" if len(set(directions)) == 1 else "mixed"
        row_out: dict[str, Any] = {
            "region_id": target_row["region_id"],
            "chrom": chrom,
            "start": start,
            "end": end,
            "context": context,
            "target_delta": target_row.get(delta_col) if delta_col else pd.NA,
            "target_q_value": target_row.get(q_col) if q_col else pd.NA,
            "supporting_callers": ",".join(sorted(support_callers)),
            "n_callers_supporting": len(support_callers),
            "best_external_q_value": min(q_by_caller.values()) if q_by_caller else pd.NA,
            "n_callers_q05": sum(1 for value in q_by_caller.values() if value <= 0.05),
            "n_callers_q10": sum(1 for value in q_by_caller.values() if value <= 0.10),
            "max_external_abs_delta": max((abs(value) for value in delta_by_caller.values()), default=math.nan),
            "direction_consensus": direction_consensus,
            "caller_conflict_flag": direction_consensus == "mixed",
            "mean_reciprocal_overlap": pd.Series(overlaps).mean() if overlaps else pd.NA,
            "min_reciprocal_overlap": min(overlaps) if overlaps else pd.NA,
        }
        for col in carry_cols:
            row_out[col] = target_row.get(col)
        for caller in callers:
            row_out[support_column(caller)] = caller in support_callers
            row_out[value_column(caller, "_q_value")] = q_by_caller.get(caller, pd.NA)
            row_out[value_column(caller, "_delta")] = delta_by_caller.get(caller, pd.NA)
        rows.append(row_out)
    return pd.DataFrame(rows)


def write_run_plan(out_dir: Path, *, external_root: Path | None, contexts: list[str]) -> None:
    lines = [
        "# External DMR Caller Run Plan",
        "",
        "This framework step imports caller outputs; it does not execute third-party callers.",
        "",
        "Required external tools for caller execution:",
        "",
        "- R/Bioconductor: DSS, methylKit, dmrseq, bsseq/BSmooth, DMRcate/limma",
        "- metilene binary on PATH for metilene",
        "- per-site p-value table for comb-p-like regional merging",
        "",
        "Expected caller output locations for auto-discovery:",
        "",
    ]
    root = external_root or Path("<external-root>")
    for caller, patterns in DEFAULT_CALLER_PATTERNS.items():
        for context in contexts:
            lines.append(f"- {caller} {context}: `{root / patterns[0].format(context=context)}`")
    (out_dir / "external_caller_run_plan.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(
    *,
    out_dir: Path,
    target_path: Path,
    external_root: Path | None,
    specs: list[CallerOutputSpec],
    inventory: list[dict[str, Any]],
    external: pd.DataFrame,
    support: pd.DataFrame,
    args: argparse.Namespace,
    contexts: list[str],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    inventory_df = pd.DataFrame(inventory, columns=INVENTORY_COLUMNS)
    inventory_df.to_csv(out_dir / "external_caller_inventory.tsv", sep="\t", index=False)
    if external.empty:
        pd.DataFrame(columns=canonical_dmr_columns()).to_csv(
            out_dir / "external_dmr_candidates_canonical.tsv",
            sep="\t",
            index=False,
        )
    else:
        external.to_csv(out_dir / "external_dmr_candidates_canonical.tsv", sep="\t", index=False)
    support.to_csv(out_dir / "external_caller_support_matrix.tsv", sep="\t", index=False)
    summary = {
        "n_target_regions": int(len(support)),
        "n_regions_with_external_support": int((support.get("n_callers_supporting", pd.Series(dtype=int)) > 0).sum()),
        "fraction_regions_with_external_support": (
            float((support["n_callers_supporting"] > 0).mean()) if len(support) else 0.0
        ),
        "max_callers_supporting": int(support["n_callers_supporting"].max()) if len(support) else 0,
    }
    pd.DataFrame([{column: summary[column] for column in SUMMARY_COLUMNS}]).to_csv(
        out_dir / "external_caller_support_summary.tsv",
        sep="\t",
        index=False,
    )
    if args.write_run_plan:
        write_run_plan(out_dir, external_root=external_root, contexts=contexts)
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "target_dmr_table": str(target_path),
        "external_root": str(external_root) if external_root else None,
        "contexts": contexts,
        "overlap_threshold": args.overlap_threshold,
        "include_candidate_only": args.include_candidate_only,
        "max_target_regions": args.max_target_regions,
        "caller_outputs": [
            {
                "caller": spec.caller,
                "context": spec.context,
                "path": str(spec.path),
                "source": spec.source,
            }
            for spec in specs
        ],
        "summary": summary,
        "outputs": sorted(
            {
                "external_caller_harmonization_manifest.json",
                "external_caller_inventory.tsv",
                "external_caller_support_matrix.tsv",
                "external_caller_support_summary.tsv",
                "external_dmr_candidates_canonical.tsv",
                *(
                    ["external_caller_run_plan.md"]
                    if args.write_run_plan
                    else []
                ),
            }
        ),
        "method_note": (
            "External DMR callers are not executed by this workflow. The workflow imports already-produced "
            "caller outputs through canonical adapters and overlaps them with the target DMR table."
        ),
    }
    (out_dir / "external_caller_harmonization_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )


def run(args: argparse.Namespace) -> int:
    _, target_path, external_root, out_dir = default_paths(args)
    if target_path is None or out_dir is None:
        raise SystemExit("--root or both --target-dmr-table and --out-dir are required")
    contexts = contexts_from_args(args.context)
    specs = collect_output_specs(args, external_root, contexts)
    external, inventory = import_external_tables(
        specs,
        contrast_id=args.contrast_id,
        condition_a=args.condition_a,
        condition_b=args.condition_b,
    )
    target = load_target_regions(
        target_path,
        include_candidate_only=args.include_candidate_only,
        max_target_regions=args.max_target_regions,
    )
    support = build_target_support_matrix(
        target,
        external,
        overlap_threshold=args.overlap_threshold,
    )
    write_outputs(
        out_dir=out_dir,
        target_path=target_path,
        external_root=external_root,
        specs=specs,
        inventory=inventory,
        external=external,
        support=support,
        args=args,
        contexts=contexts,
    )
    print(f"external caller outputs discovered/imported: {len(specs)}")
    print(f"canonical external rows: {len(external)}")
    print(f"target regions scored: {len(support)}")
    print(f"wrote external caller support: {out_dir / 'external_caller_support_matrix.tsv'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
