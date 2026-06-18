"""Reusable IO, discovery, and table-loading helpers."""

from __future__ import annotations

import csv
import gzip
import html
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from dmr_validation_framework.core.columns import (
    ALIASES,
    canonical_caller_name,
    first_existing,
    normalize_region_table,
    significant_or_all,
)

DEFAULT_OUT_DIR = Path("outputs/validation_audit")
SEARCH_ROOTS_ENV = "DMR_VALIDATION_SEARCH_ROOTS"

SKIP_DIRS = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    "__pycache__",
    "target",
    "node_modules",
    "cargo_home",
    "conda_pkgs",
    "envs",
    "tools",
    "validation_audit",
    "validation_audit_precheck",
}


def default_search_roots() -> list[Path]:
    roots: list[Path] = []
    env_value = os.environ.get(SEARCH_ROOTS_ENV, "")
    if env_value:
        roots.extend(Path(token) for token in env_value.split(os.pathsep) if token)
    roots.extend([Path.cwd(), Path.cwd() / "outputs"])
    seen: set[Path] = set()
    out: list[Path] = []
    for root in roots:
        root = root.expanduser()
        if not root.exists():
            continue
        resolved = root.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        out.append(root)
    return out


def ensure_out_dir(path: Path | str | None = None) -> Path:
    out = Path(path) if path else DEFAULT_OUT_DIR
    out.mkdir(parents=True, exist_ok=True)
    return out


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def count_rows(path: Path, max_count: int = 200_000) -> str:
    try:
        with open_text(path) as fh:
            next(fh, None)
            n = 0
            for n, _ in enumerate(fh, start=1):
                if n >= max_count:
                    return f">={max_count}"
            return str(n)
    except Exception:
        return "NA"


def read_table(path: str | Path, nrows: int | None = None) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    suffixes = [s.lower() for s in path.suffixes]
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".feather", ".arrow"}:
        return pd.read_feather(path)
    sep = "," if suffix == ".csv" or ".csv" in suffixes else "	"
    try:
        return pd.read_csv(path, sep=sep, nrows=nrows, low_memory=False)
    except Exception:
        return pd.read_csv(path, sep=None, engine="python", nrows=nrows)


def read_table_auto(path: str | Path, nrows: int | None = None) -> pd.DataFrame:
    return read_table(path, nrows=nrows)


def write_table(path: str | Path, data: pd.DataFrame) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        data.to_parquet(path, index=False)
    elif suffix in {".feather", ".arrow"}:
        data.to_feather(path)
    elif suffix == ".csv":
        data.to_csv(path, index=False)
    else:
        data.to_csv(path, sep="	", index=False)

def find_candidate_files(
    roots: Sequence[Path] | None = None,
    max_depth: int = 8,
    max_files: int = 50_000,
) -> list[Path]:
    roots = list(roots or default_search_roots())
    found: list[Path] = []
    seen: set[Path] = set()
    allowed_suffixes = (".tsv", ".csv", ".txt", ".bed", ".gff", ".gff3", ".gz", ".parquet", ".feather", ".arrow")
    role_tokens = (
        "dmr",
        "glm",
        "glmm",
        "dmrseq",
        "radmeth",
        "moabs",
        "methylsig",
        "bsmooth",
        "bsseq",
        "biseq",
        "dmrcate",
        "methylasso",
        "hmm-dm",
        "hmmdm",
        "comb-p",
        "combp",
        "metilene",
        "occupancy",
        "mapped_dmr",
        "unmapped_dmr",
        "density",
        "gene",
        "genes",
        "annotation",
        "dss",
        "methylkit",
        "regional_evidence",
        "support_matrix",
    )
    for root in roots:
        root = root.resolve()
        if root in seen:
            continue
        seen.add(root)
        if root.is_file():
            found.append(root)
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            current = Path(dirpath)
            try:
                rel_depth = len(current.relative_to(root).parts)
            except ValueError:
                rel_depth = 0
            dirnames[:] = [
                d for d in dirnames if d not in SKIP_DIRS and not d.startswith(".")
            ]
            if rel_depth >= max_depth:
                dirnames[:] = []
            for filename in filenames:
                path = current / filename
                lower = str(path).lower()
                if not lower.endswith(allowed_suffixes):
                    continue
                if not any(token in lower for token in role_tokens):
                    continue
                found.append(path)
                if len(found) >= max_files:
                    return _unique_paths(found)
    return _unique_paths(found)


def _unique_paths(paths: Iterable[Path]) -> list[Path]:
    out: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        key = path.resolve()
        if key in seen:
            continue
        seen.add(key)
        out.append(path)
    return out

def classify_file_role(path: Path) -> str:
    lower = str(path).lower()
    name = path.name.lower()
    known_external_callers = (
        "dmrseq",
        "radmeth",
        "moabs",
        "methylsig",
        "bsmooth",
        "bsseq",
        "biseq",
        "dmrcate",
        "methylasso",
        "hmm-dm",
        "hmmdm",
        "comb-p",
        "combp",
        "metilene",
    )
    if "occupancy_count_matrix" in name:
        return "occupancy_count_matrix"
    if "occupancy_matrix" in name:
        return "occupancy_matrix"
    if "mapped_dmr_centers" in name:
        return "mapped_dmr_centers"
    if "unmapped_dmr_centers" in name:
        return "unmapped_dmr_centers"
    if "density" in name and "dmr" in lower:
        return "density_profile"
    if "glm_vs_glmm" in name:
        return "glm_vs_glmm_comparison"
    if "glmm" in lower and name.endswith((".tsv", ".csv", ".txt")):
        return "glmm_results"
    if ("aggregated_glm" in lower or "beta_binomial_glm" in lower) and name.endswith(
        (".tsv", ".csv", ".txt")
    ):
        return "glm_results"
    if "support_matrix" in name or "method_support" in name:
        return "method_support_matrix"
    if "methylkit" in lower and name.endswith((".tsv", ".csv", ".txt")):
        return "methylkit_dmr_table"
    if "dss" in lower and name.endswith((".tsv", ".csv", ".txt")):
        return "dss_dmr_table"
    if any(token in lower for token in known_external_callers) and name.endswith((".tsv", ".csv", ".txt")):
        return "dmr_input_table"
    if "regional_evidence" in lower and name.endswith((".tsv", ".csv", ".txt")):
        return "regional_evidence_dmr_table"
    if "dmr_regions_canonical" in name or "dmr_regions_annotated" in name:
        return "canonical_dmr_table"
    if "high_confidence_dmr" in name or "significant_dmr" in name:
        return "dmr_input_table"
    if name.endswith((".bed", ".gff", ".gff3", ".gff3.gz")) and (
        "gene" in lower or "annotation" in lower or "gff" in lower
    ):
        return "gene_annotation"
    if "dmr" in lower and name.endswith((".tsv", ".csv", ".txt")):
        return "dmr_input_table"
    return "other"

REQUIRED_COLUMNS: dict[str, set[str]] = {
    "canonical_dmr_table": {"chrom", "start", "end", "context"},
    "regional_evidence_dmr_table": {"chrom", "start", "end", "context"},
    "dss_dmr_table": {"chrom", "start", "end", "context"},
    "methylkit_dmr_table": {"chrom", "start", "end", "context"},
    "dmr_input_table": {"chrom", "start", "end", "context"},
    "glm_results": {"region_id"},
    "glmm_results": {"region_id"},
    "glm_vs_glmm_comparison": {"region_id"},
    "occupancy_matrix": {"gene_id"},
    "occupancy_count_matrix": {"gene_id"},
    "mapped_dmr_centers": {"chrom", "start", "end", "gene_id", "metagene_bin"},
    "unmapped_dmr_centers": {"chrom", "start", "end", "gene_id", "reason"},
    "density_profile": set(),
    "method_support_matrix": {"chrom", "start", "end", "context"},
    "gene_annotation": set(),
}

def write_tsv(path: Path, rows: Iterable[dict], fieldnames: Sequence[str] | None = None) -> None:
    rows = list(rows)
    if fieldnames is None:
        keys: list[str] = []
        for row in rows:
            for key in row:
                if key not in keys:
                    keys.append(key)
        fieldnames = keys or ["status", "notes"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fieldnames), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def write_markdown_table(path: Path, title: str, rows: list[dict], max_rows: int = 50) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# {title}\n\n")
        fh.write("Lightweight audit output. Heavy raw FASTQ/Bismark pipelines are not run.\n\n")
        if not rows:
            fh.write("No rows.\n")
            return
        cols = list(rows[0].keys())
        fh.write("| " + " | ".join(cols) + " |\n")
        fh.write("| " + " | ".join(["---"] * len(cols)) + " |\n")
        for row in rows[:max_rows]:
            fh.write("| " + " | ".join(str(row.get(col, "")).replace("|", "\\|") for col in cols) + " |\n")
        if len(rows) > max_rows:
            fh.write(f"\nShowing {max_rows} of {len(rows)} rows.\n")

@dataclass
class CallerTable:
    name: str
    path: Path
    data: pd.DataFrame
    note: str
    source_kind: str = "native"

def load_caller_tables(roots: Sequence[Path] | None = None, *, include_support_matrix: bool = True) -> list[CallerTable]:
    files = find_candidate_files(roots)
    known_caller_tokens = [
        "regional_evidence",
        "bsx2",
        "dss",
        "methylkit",
        "dmrseq",
        "radmeth",
        "moabs",
        "methylsig",
        "bsmooth",
        "bsseq",
        "biseq",
        "dmrcate",
        "methylasso",
        "hmm-dm",
        "hmmdm",
        "comb-p",
        "combp",
        "metilene",
    ]
    known_caller_names = {
        "Regional Evidence",
        "DSS",
        "methylKit",
        "dmrseq",
        "RADMeth",
        "MOABS",
        "methylSig",
        "BSmooth",
        "BiSeq",
        "DMRcate",
        "MethyLasso",
        "HMM-DM",
        "comb-p",
        "metilene",
    }
    preferred: list[Path] = []
    for path in files:
        lower = str(path).lower()
        if path.name == "dmr_regions_canonical.tsv" and "harmonized" in lower:
            preferred.insert(0, path)
        elif path.name in {"dmr_method_support_matrix.tsv", "dmr_regions_annotated.tsv"}:
            preferred.append(path)
        elif (
            any(token in lower for token in known_caller_tokens if token != "bsx2") or "bsx2" in path.name.lower()
        ) and path.suffix in {
            ".tsv",
            ".csv",
            ".txt",
        }:
            preferred.append(path)

    tables: dict[str, CallerTable] = {}
    source_priority = {"support_matrix": 0, "native": 1}

    def store_table(cname: str, path: Path, filtered: pd.DataFrame, note: str, source_kind: str) -> None:
        current = tables.get(cname)
        priority = source_priority.get(source_kind, 0)
        current_priority = source_priority.get(current.source_kind, 0) if current is not None else -1
        if current is None or priority > current_priority or (
            priority == current_priority and len(filtered) > len(current.data)
        ):
            tables[cname] = CallerTable(cname, path, filtered, note, source_kind)

    for path in preferred[:30]:
        try:
            raw = read_table(path)
        except Exception:
            continue
        support_specs = [
            ("Regional Evidence", "regional_evidence_significant", ["region_q_value", "q_value"]),
            ("DSS", "dss_significant", ["dss_q_value", "dss_q"]),
            ("methylKit", "methylkit_significant", ["methylkit_q_value", "methylkit_q"]),
            ("dmrseq", "dmrseq_significant", ["dmrseq_q_value", "dmrseq_q", "dmrseq_qval"]),
            ("RADMeth", "radmeth_significant", ["radmeth_q_value", "radmeth_q", "radmeth_fdr"]),
            ("MOABS", "moabs_significant", ["moabs_q_value", "moabs_q", "moabs_fdr"]),
            ("methylSig", "methylsig_significant", ["methylsig_q_value", "methylsig_q", "methylsig_fdr"]),
            ("BSmooth", "bsmooth_significant", ["bsmooth_q_value", "bsmooth_q", "bsmooth_fdr"]),
            ("BiSeq", "biseq_significant", ["biseq_q_value", "biseq_q", "biseq_fdr"]),
            ("DMRcate", "dmrcate_significant", ["dmrcate_q_value", "dmrcate_q", "dmrcate_fdr", "dmrcate_minfdr"]),
            ("MethyLasso", "methylasso_significant", ["methylasso_q_value", "methylasso_q", "methylasso_fdr"]),
            ("HMM-DM", "hmmdm_significant", ["hmmdm_q_value", "hmmdm_q", "hmmdm_fdr"]),
            ("comb-p", "combp_significant", ["combp_q_value", "combp_q", "combp_fdr", "combp_sidak"]),
        ]
        if include_support_matrix and {"chrom", "start", "end", "context"}.issubset(raw.columns) and any(
            sig_col in raw.columns for _, sig_col, _ in support_specs
        ):
            base_cols = ["chrom", "start", "end", "context"]
            rid_col = first_existing(raw.columns, ["region_id", "dmr_id", "harmonized_region_id"])
            delta_col = first_existing(raw.columns, ALIASES["delta"])
            for cname, sig_col, q_aliases in support_specs:
                if sig_col not in raw.columns:
                    continue
                mask = raw[sig_col].astype(str).str.lower().isin(["true", "1", "yes"])
                sub = raw[mask].copy()
                if sub.empty:
                    continue
                q_col = first_existing(raw.columns, q_aliases)
                norm = pd.DataFrame(
                    {
                        "region_id": sub[rid_col] if rid_col else [f"{cname}:{i}" for i in range(len(sub))],
                        "caller": cname,
                        "chrom": sub["chrom"],
                        "start": pd.to_numeric(sub["start"], errors="coerce"),
                        "end": pd.to_numeric(sub["end"], errors="coerce"),
                        "context": sub["context"],
                    }
                )
                if delta_col:
                    norm["delta"] = pd.to_numeric(sub[delta_col], errors="coerce")
                if q_col:
                    norm["q_value"] = pd.to_numeric(sub[q_col], errors="coerce")
                norm = norm.dropna(subset=["start", "end"]).copy()
                norm = norm[norm["start"] < norm["end"]].reset_index(drop=True)
                filtered, sig_note = significant_or_all(norm)
                store_table(
                    cname,
                    path,
                    filtered,
                    f"loaded from method support matrix column {sig_col}; {sig_note}",
                    "support_matrix",
                )
        norm, notes = normalize_region_table(raw, source=path.stem)
        if norm.empty:
            continue
        if "caller" in norm.columns and norm["caller"].nunique(dropna=True) > 1:
            for caller_value, sub in norm.groupby("caller", dropna=True):
                cname = canonical_caller_name(str(caller_value))
                filtered, sig_note = significant_or_all(sub)
                store_table(cname, path, filtered, "; ".join(notes + [sig_note]), "native")
        else:
            lower = str(path).lower()
            name_lower = path.name.lower()
            cname = None
            caller_values = norm["caller"].dropna().astype(str).unique() if "caller" in norm.columns else []
            if len(caller_values) == 1:
                inferred = canonical_caller_name(str(caller_values[0]))
                if inferred in known_caller_names:
                    cname = inferred
            if cname is None and ("regional_evidence" in lower or "bsx2" in name_lower):
                cname = "Regional Evidence"
            elif cname is None and "dss" in lower:
                cname = "DSS"
            elif cname is None and "methylkit" in lower:
                cname = "methylKit"
            elif cname is None and "dmrseq" in lower:
                cname = "dmrseq"
            elif cname is None and "radmeth" in lower:
                cname = "RADMeth"
            elif cname is None and "moabs" in lower:
                cname = "MOABS"
            elif cname is None and "methylsig" in lower:
                cname = "methylSig"
            elif cname is None and ("bsmooth" in lower or "bsseq" in lower):
                cname = "BSmooth"
            elif cname is None and "biseq" in lower:
                cname = "BiSeq"
            elif cname is None and "dmrcate" in lower:
                cname = "DMRcate"
            elif cname is None and "methylasso" in lower:
                cname = "MethyLasso"
            elif cname is None and ("hmm-dm" in lower or "hmmdm" in lower):
                cname = "HMM-DM"
            elif cname is None and ("comb-p" in lower or "combp" in lower):
                cname = "comb-p"
            elif cname is None and "metilene" in lower:
                cname = "metilene"
            if cname:
                filtered, sig_note = significant_or_all(norm)
                store_table(cname, path, filtered, "; ".join(notes + [sig_note]), "native")
    priority = [
        "Regional Evidence",
        "DSS",
        "methylKit",
        "dmrseq",
        "RADMeth",
        "MOABS",
        "methylSig",
        "BSmooth",
        "BiSeq",
        "DMRcate",
        "MethyLasso",
        "HMM-DM",
        "comb-p",
        "metilene",
    ]
    ordered = [tables[k] for k in priority if k in tables]
    ordered.extend(tables[k] for k in sorted(tables) if k not in set(priority))
    return ordered

def find_first(role: str, roots: Sequence[Path] | None = None) -> Path | None:
    for path in find_candidate_files(roots):
        if classify_file_role(path) == role:
            return path
    return None

def find_preferred_file(names: Sequence[str], roots: Sequence[Path] | None = None) -> Path | None:
    files = find_candidate_files(roots)
    for wanted in names:
        for path in files:
            if path.name == wanted:
                return path
    for wanted in names:
        for path in files:
            if wanted.lower() in str(path).lower():
                return path
    return None


def write_html_table(path: str | Path, title: str, rows: Iterable[dict], *, max_rows: int = 200) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    shown = rows[:max_rows]
    header = "".join(f"<th>{html.escape(str(col))}</th>" for col in columns)
    body = "\n".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(col, '')))}</td>" for col in columns) + "</tr>"
        for row in shown
    )
    overflow = "" if len(rows) <= max_rows else f"<p>Showing {max_rows} of {len(rows)} rows.</p>"
    path.write_text(
        "\n".join(
            [
                "<!doctype html>",
                "<html><head><meta charset=\"utf-8\">",
                f"<title>{html.escape(title)}</title>",
                "<style>body{font-family:Arial,sans-serif;margin:24px}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:5px 7px;font-size:12px}th{background:#f4f4f4}</style>",
                "</head><body>",
                f"<h1>{html.escape(title)}</h1>",
                overflow,
                f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>",
                "</body></html>",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
