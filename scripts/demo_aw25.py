#!/usr/bin/env python3
"""
AW25 one-shot demo on real data (BSX + GFF/BED).

Generates:
- Metagene (classic/arbitrary) line/heatmap/box/violin
- Clustering/PCA/dendrogram/ordered heatmap
- Chrmap line/violin (if windows TSV provided)
- index.html with links

Defaults use bundled test data (bsxplorer2/tests/data/report.bsx, annot.gff),
override via CLI flags or env vars: BSX_FILE, ANNOT_FILE, WINDOWS_FILE.
"""

from __future__ import annotations

import argparse
import os
import pathlib
from typing import List, Sequence, Tuple

import holoviews as hv
import numpy as np
import pandas as pd
import plotly.figure_factory as ff
import polars as pl

from bsx2.io import RegionReader
from bsx2._bsx2 import HcAnnotStore, AggMethod
from bsx2.plots import Segment
from bsx2.plots.metagene import compute_discrete_regions, collect_contigs_from_hcannot
from bsx2.plots.polars_html import line_html, heatmap_html, box_html, violin_html
from bsx2.plots.cluster import prepare_matrix, run_kmeans, run_pca, run_linkage, reorder_matrix
from bsx2.plots.chrmap import prepare_chr_line_data, prepare_chr_box_data
from bsx2.plots.chrmap_vis import chr_line_hv, chr_box_hv

hv.extension("bokeh")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def resolve_default_paths() -> Tuple[pathlib.Path, pathlib.Path]:
    """Return bundled test bsx/gff paths relative to repo root."""
    here = pathlib.Path(__file__).resolve()
    repo_root = here.parent.parent
    default_bsx = repo_root / "bsxplorer2" / "tests" / "data" / "report.bsx"
    default_annot = repo_root / "bsxplorer2" / "tests" / "data" / "annot.gff"
    return default_bsx, default_annot


def ensure_dir(path: pathlib.Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def save_html(text: str, path: pathlib.Path) -> None:
    path.write_text(text, encoding="utf-8")


def save_hv(plot, path: pathlib.Path) -> None:
    hv.save(plot, str(path))


def write_index(out_dir: pathlib.Path, entries: List[Tuple[str, str]], meta: dict) -> None:
    lines = ["<!DOCTYPE html>", "<html><head><meta charset='utf-8'><title>AW25 artifacts</title></head><body>"]
    lines.append("<h1>AW25 Demo Artifacts</h1>")
    lines.append("<p><strong>BSX:</strong> {bsx}<br><strong>Annot:</strong> {annot}<br><strong>Windows:</strong> {windows}<br><strong>Regions:</strong> {regions}<br><strong>Segments classic:</strong> {classic}<br><strong>Segments arbitrary:</strong> {arbitrary}</p>".format(
        bsx=meta.get("bsx", ""),
        annot=meta.get("annot", ""),
        windows=meta.get("windows", "—"),
        regions=meta.get("regions", ""),
        classic=meta.get("classic", ""),
        arbitrary=meta.get("arbitrary", ""),
    ))
    lines.append("<ul>")
    for label, fname in entries:
        lines.append(f"<li><a href='{fname}' target='_blank'>{label}</a></li>")
    lines.append("</ul></body></html>")
    (out_dir / "index.html").write_text("\n".join(lines), encoding="utf-8")


# -----------------------------------------------------------------------------
# Data prep
# -----------------------------------------------------------------------------
def load_sources(bsx_path: pathlib.Path, annot_path: pathlib.Path) -> Tuple[RegionReader, HcAnnotStore]:
    if not bsx_path.exists():
        raise FileNotFoundError(f"BSX file not found: {bsx_path}")
    if not annot_path.exists():
        raise FileNotFoundError(f"Annotation file not found: {annot_path}")
    reader = RegionReader(str(bsx_path))
    annot = HcAnnotStore.from_gff(str(annot_path))
    return reader, annot


def get_contigs(annot: HcAnnotStore, feature: str | None, limit: int | None) -> Tuple[Sequence, Sequence[str]]:
    contigs, labels = collect_contigs_from_hcannot(annot, feature_type=feature, limit=limit)
    return contigs, labels


def _sanitize_drd(drd):
    """Clean infinities; keep NaN (no data)."""
    from bsx2.plots.data import DiscreteRegionData

    clean = DiscreteRegionData()
    for x, y, lbl in zip(drd.positions, drd.densities, drd.labels):
        y_clean = y.astype(float, copy=True)
        bad = ~np.isfinite(y_clean)
        if np.any(bad):
            y_clean[bad] = np.nan
        mask = np.isfinite(y_clean)
        if np.any(mask):
            y_clean[mask] = np.clip(y_clean[mask], 0.0, 1.0)
        clean.insert(x, y_clean, lbl)
    return clean


def drd_from_contigs(
    reader: RegionReader,
    contigs: Sequence,
    labels: Sequence[str] | None,
    segments: List[Segment],
):
    raw = compute_discrete_regions(
        reader,
        contigs,
        segments=segments,
        agg_method=AggMethod.Mean,
        reverse_negative=True,
        labels=labels,
    )
    return _sanitize_drd(raw)


def _auto_windows(reader: RegionReader, contigs: Sequence, window_size: int = 100_000, max_contigs: int = 3) -> pl.DataFrame | None:
    """Generate a simple windows table from BSX if external TSV is not provided."""
    try:
        from itertools import islice
    except ImportError:
        pass
    rows = []
    contigs_sel = list(contigs)[:max_contigs]
    for batch in reader.iter_contigs(contigs_sel):
        df = batch.data() if hasattr(batch, "data") else batch.into_dataframe()
        cols = set(df.columns)
        pos_col = "pos" if "pos" in cols else "position" if "position" in cols else None
        m_col = "count_m" if "count_m" in cols else "sum_m" if "sum_m" in cols else None
        t_col = "count_total" if "count_total" in cols else "sum_t" if "sum_t" in cols else None
        if not pos_col or not m_col or not t_col:
            continue
        chr_name = getattr(batch, "chr", None) if hasattr(batch, "chr") else None
        if callable(chr_name):
            try:
                chr_name = chr_name()
            except Exception:
                chr_name = None
        if chr_name is None:
            for attr in ("chrom", "name"):
                val = getattr(batch, attr, None)
                if val is not None:
                    chr_name = val() if callable(val) else val
                    break
        if chr_name is None and hasattr(batch, "contig"):
            val = getattr(batch.contig, "seqname", None)
            chr_name = val() if callable(val) else val
        # Normalize chr_name to plain string
        if chr_name is None:
            chr_lit = "chr"
        elif isinstance(chr_name, str):
            chr_lit = chr_name
        else:
            try:
                import polars as pl  # local import to avoid hard dep here
                if isinstance(chr_name, pl.Series):
                    chr_lit = str(chr_name.item()) if chr_name.len() > 0 else "chr"
                else:
                    chr_lit = str(chr_name)
            except Exception:
                chr_lit = str(chr_name)
        df = df.select([pos_col, m_col, t_col]).rename({pos_col: "pos", m_col: "sum_m", t_col: "sum_t"})
        df = df.with_columns((pl.col("pos") // window_size).alias("window"))
        df = df.with_columns(pl.lit(chr_lit, allow_object=True).alias("chr"))
        rows.append(df.select(["chr", "window", "sum_m", "sum_t"]))
    if not rows:
        return None
    return (
        pl.concat(rows, how="vertical_relaxed")
        .group_by(["chr", "window"], maintain_order=True)
        .agg([pl.col("sum_m").sum().alias("sum"), pl.col("sum_t").sum().alias("count")])
        .sort(["chr", "window"])
    )


# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------
def metagene_htmls(drd, name: str, out_dir: pathlib.Path, files: list[Tuple[str, str]], segments: Sequence[Segment]) -> None:
    line = line_html(drd, agg="mean", full_html=True, segments=list(segments), as_percent=True)
    heat = heatmap_html(drd, full_html=True, segments=list(segments), as_percent=True)
    box = box_html(drd, full_html=True, segments=list(segments), as_percent=True)
    violin = violin_html(drd, full_html=True, segments=list(segments), as_percent=True)
    outputs = [
        (f"metagene_{name}_line.html", line, f"Metagene {name} line"),
        (f"metagene_{name}_heatmap.html", heat, f"Metagene {name} heatmap"),
        (f"metagene_{name}_box.html", box, f"Metagene {name} box"),
        (f"metagene_{name}_violin.html", violin, f"Metagene {name} violin"),
    ]
    for fname, content, label in outputs:
        save_html(content, out_dir / fname)
        files.append((label, fname))


def clustering_plots(drd, seed: int, out_dir: pathlib.Path, files: list[Tuple[str, str]]) -> None:
    mat = prepare_matrix(drd, norm="zscore")
    if mat.matrix.size == 0:
        print("Clustering: empty matrix, skipping")
        return

    try:
        kmeans = run_kmeans(mat, n_clusters=4, random_state=seed)
    except ValueError as e:
        print(f"Clustering skipped: {e}")
        return
    pca = run_pca(mat, n_components=2)
    scatter_df = pd.DataFrame(
        {"pc1": pca.scores[:, 0], "pc2": pca.scores[:, 1], "cluster": kmeans.labels}
    )
    scatter = hv.Points(scatter_df, kdims=["pc1", "pc2"], vdims="cluster").opts(
        color="cluster", cmap="Category10", size=7, height=400, width=550, title="PCA scatter"
    )
    save_hv(scatter, out_dir / "cluster_pca_scatter.html")
    files.append(("PCA scatter", "cluster_pca_scatter.html"))

    # Centroids heatmap
    cent_df = (
        pd.DataFrame(kmeans.centroids)
        .reset_index()
        .melt(id_vars="index", var_name="bin", value_name="value")
    )
    cent = hv.HeatMap(cent_df, kdims=["index", "bin"], vdims="value").opts(
        title="Cluster centroids heatmap", cmap="Viridis", colorbar=True, height=400, width=600
    )
    save_hv(cent, out_dir / "cluster_centroids_heatmap.html")
    files.append(("Centroid heatmap", "cluster_centroids_heatmap.html"))

    # Linkage + ordered heatmap
    try:
        linkage = run_linkage(mat)
        ordered = reorder_matrix(mat, linkage.order)
        order_df = (
            pd.DataFrame(ordered.matrix)
            .reset_index()
            .melt(id_vars="index", var_name="bin", value_name="value")
        )
        ordered_hm = hv.HeatMap(order_df, kdims=["index", "bin"], vdims="value").opts(
            title="Ordered heatmap", cmap="Viridis", colorbar=True, height=400, width=600
        )
        save_hv(ordered_hm, out_dir / "cluster_ordered_heatmap.html")
        files.append(("Ordered heatmap", "cluster_ordered_heatmap.html"))

        fig = ff.create_dendrogram(mat.matrix, orientation="left", labels=mat.region_ids)
        fig.update_layout(width=600, height=500, title="Dendrogram")
        fig.write_html(out_dir / "cluster_dendrogram.html", include_plotlyjs="cdn")
        files.append(("Dendrogram", "cluster_dendrogram.html"))

        # Clustermap (dendrogram + heatmap)
        import plotly.graph_objects as go
        dendro = ff.create_dendrogram(mat.matrix, orientation="left", labels=mat.region_ids)
        for i in range(len(dendro['data'])):
            dendro['data'][i]['xaxis'] = 'x2'
            dendro['data'][i]['yaxis'] = 'y2'
        heat = go.Heatmap(
            x=list(range(ordered.matrix.shape[1])),
            y=[mat.region_ids[i] for i in linkage.order],
            z=ordered.matrix,
            colorscale="Viridis",
            colorbar=dict(title="zscore"),
            showscale=True,
        )
        heat['xaxis'] = 'x1'
        heat['yaxis'] = 'y1'
        from plotly.subplots import make_subplots
        fig_cm = make_subplots(
            rows=1, cols=2,
            column_widths=[0.2, 0.8],
            shared_yaxes=True,
            specs=[[{'type': 'heatmap'}, {'type': 'xy'}]],
            horizontal_spacing=0.02,
        )
        fig_cm.add_trace(heat, 1, 1)
        for tr in dendro['data']:
            fig_cm.add_trace(tr, 1, 2)
        fig_cm.update_layout(
            title="Clustermap (dendrogram + ordered heatmap)",
            xaxis_title="bin",
            yaxis_title="region",
            xaxis2_title="distance",
        )
        fig_cm.update_yaxes(autorange="reversed", showticklabels=True, row=1, col=1)
        fig_cm.update_yaxes(showticklabels=False, row=1, col=2)
        fig_cm.write_html(out_dir / "cluster_clustermap.html", include_plotlyjs="cdn")
        files.append(("Clustermap", "cluster_clustermap.html"))
    except ImportError as exc:
        print(f"Linkage skipped (scipy missing): {exc}")


def chrmap_plots(windows_path: pathlib.Path | None, df: pl.DataFrame | None, out_dir: pathlib.Path, files: list[Tuple[str, str]]) -> None:
    if df is None and (not windows_path or not windows_path.exists()):
        print("Chrmap: windows file not provided, skipping")
        return
    if df is None:
        df = pl.read_csv(windows_path, separator="\t")
    line = prepare_chr_line_data(df, smooth=5)
    box = prepare_chr_box_data(df)
    line_plot = chr_line_hv(line, label="density")
    violin_plot = chr_box_hv(box, kind="violin")
    save_hv(line_plot, out_dir / "chrmap_line.html")
    save_hv(violin_plot, out_dir / "chrmap_violin.html")
    files.append(("Chrmap line", "chrmap_line.html"))
    files.append(("Chrmap violin", "chrmap_violin.html"))


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    default_bsx, default_annot = resolve_default_paths()
    parser = argparse.ArgumentParser(description="Generate AW25 artifacts from real BSX/GFF data.")
    parser.add_argument("--bsx", type=pathlib.Path, default=pathlib.Path(os.environ.get("BSX_FILE", default_bsx)), help="Path to report.bsx")
    parser.add_argument("--annot", type=pathlib.Path, default=pathlib.Path(os.environ.get("ANNOT_FILE", default_annot)), help="Path to annotation GFF/BED")
    parser.add_argument("--windows", type=pathlib.Path, default=os.environ.get("WINDOWS_FILE"), help="Optional windows TSV for chrmap")
    parser.add_argument("--out-dir", type=pathlib.Path, default=pathlib.Path("artifacts"), help="Output directory")
    parser.add_argument("--limit", type=int, default=300, help="Limit number of regions/contigs")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for clustering")
    args = parser.parse_args()

    ensure_dir(args.out_dir)

    reader, annot = load_sources(args.bsx, args.annot)
    contigs, labels = get_contigs(annot, feature="gene", limit=args.limit)
    if not contigs:
        raise RuntimeError("No contigs found in annotation for feature 'gene'.")
    print(f"Using {len(contigs)} regions from {args.annot}")
    for l, c in list(zip(labels, contigs))[:5]:
        strand = c.strand_str if hasattr(c, "strand_str") else getattr(c, "strand", "")
        print(f"  {l}: {c.seqname}:{c.start}-{c.end} {strand}")

    files: List[Tuple[str, str]] = []

    classic_segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
    arbitrary_segments = [Segment("up", 10), Segment("exon_like", 30), Segment("gap", 5), Segment("down", 55)]

    drd_classic = drd_from_contigs(reader, contigs, labels, classic_segments)
    drd_arbitrary = drd_from_contigs(reader, contigs, labels, arbitrary_segments)

    metagene_htmls(drd_classic, "classic", args.out_dir, files, classic_segments)
    metagene_htmls(drd_arbitrary, "arbitrary", args.out_dir, files, arbitrary_segments)

    clustering_plots(drd_classic, args.seed, args.out_dir, files)

    windows_df = None
    if args.windows and args.windows.exists():
        windows_df = None  # will be loaded in chrmap_plots
    else:
        windows_df = _auto_windows(reader, contigs, window_size=100_000, max_contigs=3)
        if windows_df is None:
            print("Chrmap: auto-generation failed, skipping")
    chrmap_plots(args.windows if args.windows else None, windows_df, args.out_dir, files)

    meta = {
        "bsx": args.bsx,
        "annot": args.annot,
        "windows": args.windows if args.windows else ("auto-generated" if windows_df is not None else "not provided"),
        "regions": len(contigs),
        "classic": f"{classic_segments}",
        "arbitrary": f"{arbitrary_segments}",
    }
    write_index(args.out_dir, files, meta)
    print(f"Artifacts saved to {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
