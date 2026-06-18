#!/usr/bin/env python3
"""Inventory lightweight validation inputs for the DMR/metagene workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from dmr_validation_framework.core.columns import has_columns
from dmr_validation_framework.core.io import (
    REQUIRED_COLUMNS,
    classify_file_role,
    count_rows,
    default_search_roots,
    ensure_out_dir,
    find_candidate_files,
    read_table,
    write_markdown_table,
    write_tsv,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--search-root", action="append", type=Path, default=None)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    roots = args.search_root or default_search_roots()

    rows: list[dict] = []
    for path in find_candidate_files(roots):
        role = classify_file_role(path)
        if role == "other":
            continue
        required = REQUIRED_COLUMNS.get(role, set())
        required_ok = "NA"
        missing = ""
        status = "PASS"
        notes = ""
        n_rows = "NA"
        if path.exists() and path.suffix not in {".gff", ".gff3"} and not path.name.endswith(".gff3.gz"):
            try:
                header_df = read_table(path, nrows=5)
                n_rows = count_rows(path)
                if required:
                    ok, missing = has_columns(header_df, required)
                    required_ok = str(ok)
                    if not ok:
                        status = "WARN"
                else:
                    required_ok = "NA"
            except Exception as exc:
                status = "WARN"
                notes = f"could not read table header: {exc}"
        elif path.exists():
            n_rows = "NA"
            required_ok = "NA"
            notes = "annotation-like file; header validation not applied"
        else:
            status = "SKIPPED"
            notes = "file does not exist"

        rows.append(
            {
                "file_role": role,
                "path": str(path),
                "exists": str(path.exists()),
                "n_rows": n_rows,
                "required_columns_present": required_ok,
                "missing_columns": missing,
                "status": status,
                "notes": notes,
            }
        )

    rows.sort(key=lambda row: (row["file_role"], row["path"]))
    write_tsv(out_dir / "input_inventory.tsv", rows)
    write_markdown_table(out_dir / "input_inventory.md", "DMR validation input inventory", rows)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
