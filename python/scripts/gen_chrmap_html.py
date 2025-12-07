#!/usr/bin/env python
"""Generate chromosome methylation line/box plots from a window-level table.

Input must contain columns: chr, window, sum, count and optional lower/upper.
Output HTML files land in the chosen directory (default: python/assets).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import holoviews as hv
import polars as pl

from bsx2.plots.chrmap import prepare_chr_line_data, prepare_chr_box_data
from bsx2.plots.chrmap_vis import chr_line_hv, chr_box_hv

hv.extension("bokeh", logo=False)


def _read_table(path: Path) -> pl.DataFrame:
    suf = path.suffix.lower()
    if suf in (".csv", ".tsv"):
        return pl.read_csv(path, separator="\t" if suf == ".tsv" else ",")
    if suf in (".parquet", ".pq"):
        return pl.read_parquet(path)
    raise SystemExit(f"Unsupported input format: {path}")


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path, help="CSV/TSV/Parquet with chr/window/sum/count columns")
    p.add_argument("--out-dir", default=Path("python/assets"), type=Path, help="Where to save HTML plots")
    p.add_argument("--prefix", default="chrmap", help="Prefix for output filenames")
    p.add_argument(
        "--smooth",
        type=int,
        default=0,
        help="Savitzky-Golay window length for line plot (0 to disable)",
    )
    p.add_argument("--box-kind", choices=["box", "violin"], default="box", help="Box or violin plot for per-chr densities")
    p.add_argument("--label", default=None, help="Optional label for the line plot legend")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    df = _read_table(args.input)

    line_data = prepare_chr_line_data(df, smooth=args.smooth)
    box_data = prepare_chr_box_data(df)

    line_plot = chr_line_hv(line_data, label=args.label)
    box_plot = chr_box_hv(box_data, kind=args.box_kind)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    line_path = args.out_dir / f"{args.prefix}_line.html"
    box_path = args.out_dir / f"{args.prefix}_{args.box_kind}.html"

    hv.save(line_plot, line_path, backend="bokeh", resources="cdn")
    hv.save(box_plot, box_path, backend="bokeh", resources="cdn")

    print(f"Saved {line_path}")
    print(f"Saved {box_path}")


if __name__ == "__main__":
    main()
