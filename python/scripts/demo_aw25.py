#!/usr/bin/env python
"""Generate AW25 demo artifacts (metagene, clustering, chrmap) in one shot."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import holoviews as hv
import numpy as np
import polars as pl

from bsx2.io import RegionReader
from bsx2._bsx2 import AggMethod, HcAnnotStore, BsxFileReader
from bsx2.plots import Segment, collect_contigs_from_hcannot, compute_discrete_regions
from bsx2.plots.cluster import prepare_matrix, run_kmeans, run_linkage, run_pca
from bsx2.plots.cluster_vis import (
    dendrogram_plot,
    heatmap_ordered,
    kmeans_centroids_heatmap,
    pca_scatter,
)
from bsx2.plots.metagene import segments_total_bins
from bsx2.plots.polars_html import (
    box_html_from_annot,
    heatmap_html_from_annot,
    line_html_from_annot,
    violin_html_from_annot,
)
from bsx2.plots.chrmap import prepare_chr_box_data, prepare_chr_line_data
from bsx2.plots.chrmap_vis import chr_box_hv, chr_line_hv

# Используем plotly backend по умолчанию для единого стека
hv.extension("plotly", logo=False)


def _sorted_contigs(
    annot: HcAnnotStore,
    *,
    limit: int | None,
) -> Tuple[List[object], List[str]]:
    contigs, labels = collect_contigs_from_hcannot(annot, feature_type=None, limit=limit)
    pairs = sorted(zip(labels, contigs), key=lambda x: x[0])
    if limit:
        pairs = pairs[:limit]
    lbls, cts = zip(*pairs) if pairs else ([], [])
    return list(cts), list(lbls)


def _save_html(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    print(f"Saved {path}")


def _metagene_set(
    name: str,
    segments: Sequence[Segment],
    bsx_path: Path,
    annot: HcAnnotStore,
    *,
    out_dir: Path,
    limit: int | None,
) -> None:
    rr = RegionReader(str(bsx_path))
    line = line_html_from_annot(
        rr,
        annot,
        segments=list(segments),
        agg="mean",
        agg_method=AggMethod.Mean,
        feature_type=None,
        limit=limit,
        full_html=False,
    )
    heat = heatmap_html_from_annot(
        rr,
        annot,
        segments=list(segments),
        agg_method=AggMethod.Mean,
        feature_type=None,
        limit=limit,
        full_html=False,
    )
    box = box_html_from_annot(
        rr,
        annot,
        segments=list(segments),
        agg_method=AggMethod.Mean,
        feature_type=None,
        limit=limit,
        full_html=False,
    )
    violin = violin_html_from_annot(
        rr,
        annot,
        segments=list(segments),
        agg_method=AggMethod.Mean,
        feature_type=None,
        limit=limit,
        full_html=False,
    )

    _save_html(out_dir / f"metagene_{name}_line.html", line)
    _save_html(out_dir / f"metagene_{name}_heatmap.html", heat)
    _save_html(out_dir / f"metagene_{name}_box.html", box)
    _save_html(out_dir / f"metagene_{name}_violin.html", violin)


def _cluster_artifacts(
    bsx_path: Path,
    contigs: Sequence[object],
    labels: Sequence[str],
    segments: Sequence[Segment],
    *,
    out_dir: Path,
    seed: int,
) -> None:
    rr = RegionReader(str(bsx_path))
    contig_to_label = {id(c): lbl for c, lbl in zip(contigs, labels)}
    contigs_sorted = rr.index().sort(list(contigs))
    labels_sorted = [contig_to_label.get(id(c), "") for c in contigs_sorted]
    drd = compute_discrete_regions(
        rr,
        contigs_sorted,
        segments=segments,
        agg_method=AggMethod.Mean,
        reverse_negative=True,
        labels=labels_sorted,
    )
    mat = prepare_matrix(drd, norm="zscore")
    km = run_kmeans(mat, n_clusters=4, random_state=seed)
    pca = run_pca(mat, n_components=3)
    _save_html(out_dir / "cluster_pca_scatter.html", pca_scatter(pca, labels=mat.region_ids, clusters=km.labels))
    _save_html(out_dir / "cluster_centroids_heatmap.html", kmeans_centroids_heatmap(km, bins=mat.bins))
    order = None
    try:
        from bsx2.plots.cluster import filter_constant_rows

        mat_link = filter_constant_rows(mat, drop_nan_rows=False, nan_fill=0.0)
        link = run_linkage(mat_link, drop_nan_rows=False, nan_fill=0.0)
        order = link.order
        _save_html(out_dir / "cluster_dendrogram.html", dendrogram_plot(link, labels=mat_link.region_ids))
        _save_html(out_dir / "cluster_ordered_heatmap.html", heatmap_ordered(mat_link, order=order))
    except ImportError:
        print("scipy not available; skipping dendrogram")
    if order is None:
        _save_html(out_dir / "cluster_ordered_heatmap.html", heatmap_ordered(mat, order=None))


def _windows_from_bsx(
    bsx_path: Path,
    *,
    window_size: int,
    max_batches: int | None,
    max_chroms: int = 3,
) -> pl.DataFrame:
    window_size = max(1, min(int(window_size), 1000))
    rf = BsxFileReader(str(bsx_path))
    acc: List[pl.DataFrame] = []
    seen = []
    for idx, batch in enumerate(rf):
        if max_batches is not None and idx >= max_batches:
            break
        if batch is None:
            continue
        df = batch.data()
        # filter by first N chromosomes seen
        chr_col = df["chr"].cast(pl.Utf8)
        if len(seen) < max_chroms:
            new = [c for c in chr_col.unique().to_list() if c not in seen]
            seen.extend(new)
        df = df.filter(pl.col("chr").cast(pl.Utf8).is_in(seen[:max_chroms]))
        if df.is_empty():
            continue
        df = df.with_columns((pl.col("position") // window_size).alias("window"))
        acc.append(
            df.group_by(["chr", "window"]).agg(
                [
                    pl.col("count_m").sum().alias("sum"),
                    pl.col("count_total").sum().alias("count"),
                ]
            )
        )
    if not acc:
        return pl.DataFrame({"chr": [], "window": [], "sum": [], "count": []})
    return (
        pl.concat(acc, how="vertical_relaxed")
        .group_by(["chr", "window"])
        .agg([pl.col("sum").sum(), pl.col("count").sum()])
        .sort(["chr", "window"])
    )


def _chrmap_artifacts(
    df: pl.DataFrame,
    *,
    out_dir: Path,
) -> None:
    line_data = prepare_chr_line_data(df, smooth=5)
    box_data = prepare_chr_box_data(df)

    line_plot = chr_line_hv(line_data, label="chrmap")
    box_plot = chr_box_hv(box_data, kind="violin")

    line_path = out_dir / "chrmap_line.html"
    box_path = out_dir / "chrmap_violin.html"
    hv.save(line_plot, line_path, backend="plotly", resources="cdn")
    hv.save(box_plot, box_path, backend="plotly", resources="cdn")
    print(f"Saved {line_path}")
    print(f"Saved {box_path}")


def _write_index(out_dir: Path) -> None:
    lines = [
        "# AW25 demo artifacts",
        "",
        "## Metagene classic",
        "- metagene_classic_line.html",
        "- metagene_classic_heatmap.html",
        "- metagene_classic_box.html",
        "- metagene_classic_violin.html",
        "",
        "## Metagene arbitrary",
        "- metagene_arbitrary_line.html",
        "- metagene_arbitrary_heatmap.html",
        "- metagene_arbitrary_box.html",
        "- metagene_arbitrary_violin.html",
        "",
        "## Clustering",
        "- cluster_pca_scatter.html",
        "- cluster_centroids_heatmap.html",
        "- cluster_ordered_heatmap.html",
        "- cluster_dendrogram.html",
        "",
        "## Chrmap",
        "- chrmap_line.html",
        "- chrmap_violin.html",
    ]
    (out_dir / "INDEX.md").write_text("\n".join(lines), encoding="utf-8")


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bsx", default=Path("bsxplorer2/tests/data/report.bsx"), type=Path, help="Path to BSX file")
    p.add_argument("--annot", default=Path("bsxplorer2/tests/data/annot.gff"), type=Path, help="Path to GFF/BED")
    p.add_argument("--out-dir", default=Path("artifacts"), type=Path, help="Output directory for HTML files")
    p.add_argument("--limit", type=int, default=250, help="Limit number of regions for speed")
    p.add_argument("--full", action="store_true", help="Ignore limit and process all regions")
    p.add_argument("--windows", type=Path, default=None, help="Optional precomputed windows table for chrmap")
    p.add_argument("--window-size", type=int, default=1_000, help="Window size for chrmap aggregation (bp)")
    p.add_argument("--max-batches", type=int, default=400, help="Max batches to read from BSX for chrmap fallback")
    p.add_argument("--max-chroms", type=int, default=3, help="Max chromosomes to include in chrmap fallback")
    p.add_argument("--skip-chrmap", action="store_true", help="Skip chrmap generation")
    p.add_argument("--seed", type=int, default=42, help="Random seed for clustering/determinism")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    annot = HcAnnotStore.from_gff(str(args.annot))

    limit = None if args.full else args.limit
    contigs, labels = _sorted_contigs(annot, limit=limit)

    # Сегменты по умолчанию ближе к статье: фланги короче, тело длиннее
    classic = [Segment("up", 50), Segment("body", 200), Segment("down", 50)]
    arbitrary = [Segment("up", 10), Segment("exon_like", 30), Segment("gap", 5), Segment("down", 55)]

    _metagene_set("classic", classic, args.bsx, annot, out_dir=out_dir, limit=limit)
    _metagene_set("arbitrary", arbitrary, args.bsx, annot, out_dir=out_dir, limit=limit)

    # clustering on classic segments with deterministic labels
    _cluster_artifacts(args.bsx, contigs, labels, classic, out_dir=out_dir, seed=args.seed)

    if not args.skip_chrmap:
        if args.windows is not None:
            df = pl.read_csv(args.windows) if args.windows.suffix.lower() in {".csv", ".tsv"} else pl.read_parquet(args.windows)
        else:
            df = _windows_from_bsx(
                args.bsx,
                window_size=args.window_size,
                max_batches=None if args.full else args.max_batches,
                max_chroms=args.max_chroms,
            )
        if df.is_empty():
            print("No windows available; skipping chrmap")
        else:
            _chrmap_artifacts(df, out_dir=out_dir)
    else:
        print("Chrmap generation skipped by flag")

    _write_index(out_dir)


if __name__ == "__main__":
    main()
