from __future__ import annotations

from dataclasses import dataclass, field

import holoviews as hv
import numpy as np
from beartype import beartype
from beartype.typing import Optional, Sequence
from scipy.cluster.hierarchy import dendrogram

from bsx2.clustering.models import GeneClusterResult
from bsx2.validation import validate_positive_int

from ..compute.clustering import (
    ClusterMetageneData,
    GeneDendrogramData,
    GeneEmbeddingData,
    build_cluster_metagene_data,
)
from ..compute.data import DiscreteRegionData
from ._common import _hv_init
from .line import LinePlotComposer


def _cluster_label(cluster: int, *, cluster_source: str) -> str:
    return f"{cluster_source} cluster {cluster}"


@beartype
@dataclass
class GeneEmbeddingPlotComposer:
    title: Optional[str] = None
    width: int | None = None
    height: int | None = None
    x_label: str = "PC1"
    y_label: str = "PC2"
    datasets: list[GeneEmbeddingData] = field(default_factory=list)
    labels: list[str] = field(default_factory=list)

    @beartype
    def set_width(self, width: int | None) -> GeneEmbeddingPlotComposer:
        self.width = None if width is None else validate_positive_int(width, name="width")
        return self

    @beartype
    def set_height(self, height: int | None) -> GeneEmbeddingPlotComposer:
        self.height = None if height is None else validate_positive_int(height, name="height")
        return self

    @beartype
    def add_data(
        self,
        data: GeneEmbeddingData,
        *,
        name: str = "genes",
    ) -> GeneEmbeddingPlotComposer:
        self.datasets.append(data)
        self.labels.append(name)
        return self

    def finish(self):
        _hv_init()
        if not self.datasets:
            return hv.Points([]).opts(
                title=self.title or "Gene-level PCA",
                xlabel=self.x_label,
                ylabel=self.y_label,
                show_legend=False,
                **({} if self.width is None else {"width": int(self.width)}),
                **({} if self.height is None else {"height": int(self.height)}),
            )

        overlays = []
        multi_dataset = len(self.datasets) > 1
        for dataset_name, data in zip(self.labels, self.datasets, strict=True):
            y_vals = (
                data.embedding[:, 1]
                if data.embedding.shape[1] > 1
                else np.zeros(data.embedding.shape[0], dtype=np.float64)
            )
            for cluster in sorted(np.unique(data.labels).tolist()):
                idx = np.flatnonzero(data.labels == cluster)
                label = _cluster_label(int(cluster), cluster_source=data.cluster_source)
                if multi_dataset:
                    label = f"{dataset_name}: {label}"
                points = hv.Points(
                    {
                        "PC1": data.embedding[idx, 0],
                        "PC2": y_vals[idx],
                        "gene_id": [data.gene_ids[i] for i in idx],
                        "gene_name": [data.gene_names[i] or "" for i in idx],
                        "chrom": [data.chromosomes[i] for i in idx],
                        "start": data.starts[idx],
                        "end": data.ends[idx],
                        "strand": [data.strands[i] for i in idx],
                    },
                    kdims=["PC1", "PC2"],
                    vdims=["gene_id", "gene_name", "chrom", "start", "end", "strand"],
                ).relabel(label).opts(size=7)
                overlays.append(points)

        plot = overlays[0]
        for points in overlays[1:]:
            plot *= points

        opts_kwargs = dict(
            title=self.title or "Gene-level PCA",
            xlabel=self.x_label,
            ylabel=self.y_label,
            show_legend=len(overlays) > 1,
        )
        if self.width is not None:
            opts_kwargs["width"] = int(self.width)
        if self.height is not None:
            opts_kwargs["height"] = int(self.height)
        return plot.opts(**opts_kwargs)


@beartype
@dataclass
class GeneDendrogramPlotComposer:
    title: Optional[str] = None
    width: int | None = None
    height: int | None = None
    x_label: str = "Genes"
    y_label: str = "Distance"
    datasets: list[GeneDendrogramData] = field(default_factory=list)
    labels: list[str] = field(default_factory=list)

    @beartype
    def set_width(self, width: int | None) -> GeneDendrogramPlotComposer:
        self.width = None if width is None else validate_positive_int(width, name="width")
        return self

    @beartype
    def set_height(self, height: int | None) -> GeneDendrogramPlotComposer:
        self.height = None if height is None else validate_positive_int(height, name="height")
        return self

    @beartype
    def add_data(
        self,
        data: GeneDendrogramData,
        *,
        name: str = "genes",
    ) -> GeneDendrogramPlotComposer:
        self.datasets.append(data)
        self.labels.append(name)
        return self

    def finish(self):
        _hv_init()
        if not self.datasets:
            return hv.Curve([]).opts(
                title=self.title or "Gene-level dendrogram",
                xlabel=self.x_label,
                ylabel=self.y_label,
                show_legend=False,
                **({} if self.width is None else {"width": int(self.width)}),
                **({} if self.height is None else {"height": int(self.height)}),
            )

        overlays = []
        xticks: list[tuple[float, str]] = []
        multi_dataset = len(self.datasets) > 1

        for dataset_name, data in zip(self.labels, self.datasets, strict=True):
            dendro = dendrogram(data.linkage_matrix, labels=data.leaf_labels, no_plot=True)
            if not xticks:
                xticks = [(5.0 + 10.0 * idx, label) for idx, label in enumerate(dendro["ivl"])]
            curve_label = dataset_name if multi_dataset else "dendrogram"
            for branch_index, (xs, ys) in enumerate(zip(dendro["icoord"], dendro["dcoord"], strict=True)):
                curve = hv.Curve((xs, ys), kdims=["leaf"], vdims=["distance"])
                if branch_index == 0:
                    curve = curve.relabel(curve_label)
                overlays.append(curve)

        plot = overlays[0]
        for curve in overlays[1:]:
            plot *= curve

        opts_kwargs = dict(
            title=self.title or "Gene-level dendrogram",
            xlabel=self.x_label,
            ylabel=self.y_label,
            show_legend=multi_dataset,
            xticks=xticks,
            xrotation=90,
        )
        if self.width is not None:
            opts_kwargs["width"] = int(self.width)
        if self.height is not None:
            opts_kwargs["height"] = int(self.height)
        return plot.opts(**opts_kwargs)


@beartype
def build_cluster_metagene_plot(
    data: ClusterMetageneData,
    *,
    title: str | None = None,
    width: int | None = None,
    height: int | None = None,
):
    """
    Build a HoloViews line plot for cluster metagene profiles.
    """
    composer = LinePlotComposer(
        segments=data.segments,
        n_windows=sum(segment.n_bins for segment in data.segments),
        smooth=None,
        title=title or ("Cluster metagene profile" if len(data.groups) == 1 else "Cluster metagene profiles"),
        width=width,
        height=height,
    )
    for positions, densities, label in zip(
        data.profiles.positions,
        data.profiles.densities,
        data.profiles.labels,
        strict=True,
    ):
        profile = DiscreteRegionData()
        profile.insert_unchecked(positions.copy(), densities.copy(), label)
        composer.add_data(profile, name=label or "cluster")
    return composer.finish()


@beartype
def cluster_metagene_plot(
    result: GeneClusterResult,
    *,
    cluster_ids: Sequence[int] | None = None,
    gene_ids: Sequence[str] | None = None,
    collapse: bool = False,
    label: str | None = None,
    title: str | None = None,
    width: int | None = None,
    height: int | None = None,
):
    """
    Build a cluster metagene plot directly from a clustering result.
    """
    data = build_cluster_metagene_data(
        result,
        cluster_ids=cluster_ids,
        gene_ids=gene_ids,
        collapse=collapse,
        label=label,
    )
    return build_cluster_metagene_plot(
        data,
        title=title,
        width=width,
        height=height,
    )


__all__ = [
    "GeneDendrogramPlotComposer",
    "GeneEmbeddingPlotComposer",
    "build_cluster_metagene_plot",
    "cluster_metagene_plot",
]
