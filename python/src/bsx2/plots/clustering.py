from __future__ import annotations

from dataclasses import dataclass, field
from beartype.typing import Optional

import holoviews as hv
import numpy as np
from beartype import beartype
from scipy.cluster.hierarchy import dendrogram

from bsx2.validation import validate_positive_int
from bsx2.clustering.models import FeatureBin, GeneClusterResult

from ._common import _ensure_plotly, _hv_init
from .data import DiscreteRegionData
from .line import LinePlotComposer
from .metagene import MetageneProfileSegment


@beartype
@dataclass(frozen=True)
class GeneEmbeddingData:
    gene_ids: list[str]
    gene_names: list[str | None]
    chromosomes: list[str]
    starts: np.ndarray
    ends: np.ndarray
    strands: list[str]
    labels: np.ndarray
    embedding: np.ndarray
    cluster_source: str

    def __post_init__(self) -> None:
        n_genes = len(self.gene_ids)
        if self.embedding.ndim != 2:
            raise ValueError("GeneEmbeddingData.embedding must be 2D")
        if self.embedding.shape[0] != n_genes:
            raise ValueError("GeneEmbeddingData.embedding row count must match gene_ids")
        if self.labels.shape != (n_genes,):
            raise ValueError("GeneEmbeddingData.labels shape must match gene_ids")
        if len(self.gene_names) != n_genes:
            raise ValueError("GeneEmbeddingData.gene_names length must match gene_ids")
        if len(self.chromosomes) != n_genes:
            raise ValueError("GeneEmbeddingData.chromosomes length must match gene_ids")
        if self.starts.shape != (n_genes,):
            raise ValueError("GeneEmbeddingData.starts shape must match gene_ids")
        if self.ends.shape != (n_genes,):
            raise ValueError("GeneEmbeddingData.ends shape must match gene_ids")
        if len(self.strands) != n_genes:
            raise ValueError("GeneEmbeddingData.strands length must match gene_ids")


@beartype
@dataclass(frozen=True)
class GeneDendrogramData:
    linkage_matrix: np.ndarray
    leaf_labels: list[str]
    leaf_order: np.ndarray

    def __post_init__(self) -> None:
        if self.linkage_matrix.ndim != 2 or self.linkage_matrix.shape[1] != 4:
            raise ValueError("GeneDendrogramData.linkage_matrix must have shape [n-1, 4]")
        if self.leaf_order.ndim != 1:
            raise ValueError("GeneDendrogramData.leaf_order must be 1D")
        if len(self.leaf_labels) == 0:
            raise ValueError("GeneDendrogramData.leaf_labels must not be empty")


@beartype
@dataclass(frozen=True)
class ClusterMetageneData:
    profiles: DiscreteRegionData
    segments: list[MetageneProfileSegment]


def _cluster_label(cluster: int, *, cluster_source: str) -> str:
    return f"{cluster_source} cluster {cluster}"


@beartype
def cluster_profile_segments(feature_bins: list[FeatureBin]) -> list[MetageneProfileSegment]:
    segments: list[MetageneProfileSegment] = []
    current_name: str | None = None
    current_count = 0

    for feature_bin in feature_bins:
        segment = str(getattr(feature_bin, "segment"))
        if current_name is None:
            current_name = segment
            current_count = 1
            continue
        if segment == current_name:
            current_count += 1
            continue
        segments.append(MetageneProfileSegment(current_name, current_count))
        current_name = segment
        current_count = 1

    if current_name is not None:
        segments.append(MetageneProfileSegment(current_name, current_count))
    return segments


@beartype
def build_cluster_metagene_data(result: GeneClusterResult) -> ClusterMetageneData:
    profiles = DiscreteRegionData()
    values = np.asarray(result.feature_matrix.values, dtype=float)
    labels = np.asarray(result.labels, dtype=np.int64)
    total_bins = result.feature_matrix.n_features
    x_vals = (np.arange(total_bins, dtype=np.float64) + 0.5) / float(total_bins)

    for cluster in sorted(np.unique(labels).tolist()):
        mask = labels == cluster
        cluster_values = values[mask]
        finite = np.isfinite(cluster_values)
        sums = np.where(finite, cluster_values, 0.0).sum(axis=0)
        counts = finite.sum(axis=0, dtype=np.int64)
        means = np.divide(
            sums,
            counts,
            out=np.full(total_bins, np.nan, dtype=np.float64),
            where=counts > 0,
        )
        profiles.insert_unchecked(
            x_vals.copy(),
            means.astype(np.float64, copy=False),
            _cluster_label(int(cluster), cluster_source=result.cluster_source),
        )

    return ClusterMetageneData(
        profiles=profiles,
        segments=cluster_profile_segments(result.feature_matrix.feature_bins),
    )


@beartype
def build_gene_embedding_data(result: GeneClusterResult) -> GeneEmbeddingData:
    genes = result.feature_matrix.genes
    return GeneEmbeddingData(
        gene_ids=[gene.gene_id for gene in genes],
        gene_names=[gene.gene_name for gene in genes],
        chromosomes=[gene.chrom for gene in genes],
        starts=np.asarray([gene.start for gene in genes], dtype=np.int64),
        ends=np.asarray([gene.end for gene in genes], dtype=np.int64),
        strands=[gene.strand for gene in genes],
        labels=np.asarray(result.labels, dtype=np.int64),
        embedding=np.asarray(result.embedding, dtype=np.float64),
        cluster_source=str(result.cluster_source),
    )


@beartype
def build_gene_dendrogram_data(result: GeneClusterResult) -> GeneDendrogramData | None:
    if result.linkage_matrix is None or result.leaf_order is None:
        return None
    return GeneDendrogramData(
        linkage_matrix=np.asarray(result.linkage_matrix, dtype=np.float64),
        leaf_labels=result.feature_matrix.gene_ids,
        leaf_order=np.asarray(result.leaf_order, dtype=np.int64),
    )


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
    def set_width(self, width: int | None) -> "GeneEmbeddingPlotComposer":
        self.width = None if width is None else validate_positive_int(width, name="width")
        return self

    @beartype
    def set_height(self, height: int | None) -> "GeneEmbeddingPlotComposer":
        self.height = None if height is None else validate_positive_int(height, name="height")
        return self

    @beartype
    def add_data(
        self,
        data: GeneEmbeddingData,
        *,
        name: str = "genes",
    ) -> "GeneEmbeddingPlotComposer":
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
    def to_html(
        self,
        *,
        full_html: bool = False,
        include_js: str = "cdn",
    ) -> str:
        fig = _ensure_plotly(hv.render(self.finish(), backend="plotly"))
        fig.update_layout(margin=dict(l=70, r=30, t=60, b=70))
        return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


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
    def set_width(self, width: int | None) -> "GeneDendrogramPlotComposer":
        self.width = None if width is None else validate_positive_int(width, name="width")
        return self

    @beartype
    def set_height(self, height: int | None) -> "GeneDendrogramPlotComposer":
        self.height = None if height is None else validate_positive_int(height, name="height")
        return self

    @beartype
    def add_data(
        self,
        data: GeneDendrogramData,
        *,
        name: str = "genes",
    ) -> "GeneDendrogramPlotComposer":
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
            dendro = dendrogram(
                data.linkage_matrix,
                labels=data.leaf_labels,
                no_plot=True,
            )
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
    def to_html(
        self,
        *,
        full_html: bool = False,
        include_js: str = "cdn",
    ) -> str:
        fig = _ensure_plotly(hv.render(self.finish(), backend="plotly"))
        fig.update_layout(margin=dict(l=70, r=30, t=60, b=140))
        return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


@beartype
def build_cluster_metagene_plot(data: ClusterMetageneData):
    composer = LinePlotComposer(
        segments=data.segments,
        n_windows=sum(segment.n_bins for segment in data.segments),
        smooth=None,
        title="Cluster metagene profiles",
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
