from __future__ import annotations

import numpy as np

from bsx2.clustering.models import FeatureBin, GeneAnnotation, GeneClusterResult, GeneProfileMatrix
from bsx2.viz.compute.clustering import (
    ClusterMetageneData,
    ClusterMetageneGroup,
    GeneDendrogramData,
    GeneEmbeddingData,
)
from bsx2.viz.compute.data import DiscreteRegionData
from bsx2.viz.compute.metagene import AnnotProfileLayout, AnnotProfilePart
from bsx2.viz.compute.studio_workspaces import (
    build_studio_cluster_config,
    cluster_matrix_family_key,
    cluster_plot_family_key,
    metagene_family_key,
    prepare_cluster_dendrogram_family,
    prepare_cluster_family,
    prepare_cluster_matrix_workspace,
    prepare_metagene_family,
)


def _make_layout() -> AnnotProfileLayout:
    return AnnotProfileLayout(
        (
            AnnotProfilePart("upstream", 25, source="flank5", flank_bp=2_000),
            AnnotProfilePart("body", 50, source="gene"),
            AnnotProfilePart("downstream", 25, source="flank3", flank_bp=2_000),
        )
    )


def _make_drd() -> DiscreteRegionData:
    drd = DiscreteRegionData()
    drd.insert_unchecked(
        np.asarray([0.1, 0.5, 0.9], dtype=float),
        np.asarray([0.2, 0.4, 0.3], dtype=float),
        "gene_1",
    )
    drd.insert_unchecked(
        np.asarray([0.1, 0.5, 0.9], dtype=float),
        np.asarray([0.1, 0.6, 0.2], dtype=float),
        "gene_2",
    )
    return drd


def _make_matrix() -> GeneProfileMatrix:
    genes = [
        GeneAnnotation("gene_1", "chr1", 10, 50, "+", "G1"),
        GeneAnnotation("gene_2", "chr1", 60, 100, "-", "G2"),
    ]
    feature_bins = [
        FeatureBin("metagene", "upstream", 0, 0),
        FeatureBin("metagene", "body", 0, 1),
        FeatureBin("metagene", "downstream", 0, 2),
    ]
    values = np.asarray([[0.1, 0.4, 0.2], [0.3, 0.6, 0.5]], dtype=float)
    gene_missing_rate = np.zeros(2, dtype=float)
    gene_variance = np.asarray([0.01, 0.02], dtype=float)
    feature_missing_rate = np.zeros(3, dtype=float)
    feature_variance = np.asarray([0.02, 0.03, 0.04], dtype=float)
    return GeneProfileMatrix(
        genes=genes,
        feature_bins=feature_bins,
        values=values,
        gene_missing_rate=gene_missing_rate,
        gene_variance=gene_variance,
        feature_missing_rate=feature_missing_rate,
        feature_variance=feature_variance,
    )


def _make_cluster_result(matrix: GeneProfileMatrix) -> GeneClusterResult:
    return GeneClusterResult(
        feature_matrix=matrix,
        labels=np.asarray([0, 1], dtype=np.int64),
        cluster_source="kmeans",
        kmeans_labels=np.asarray([0, 1], dtype=np.int64),
        hierarchical_labels=None,
        embedding=np.asarray([[0.1, 0.2], [0.3, 0.4]], dtype=float),
        components=np.asarray([[1.0, 0.0], [0.0, 1.0]], dtype=float),
        centroids=np.asarray([[0.1, 0.2], [0.3, 0.4]], dtype=float),
        explained_variance_ratio=np.asarray([0.6, 0.4], dtype=float),
        inertia=1.23,
        silhouette_score=0.42,
        linkage_matrix=np.asarray([[0.0, 1.0, 0.5, 2.0]], dtype=float),
        leaf_order=np.asarray([0, 1], dtype=np.int64),
        metadata={"hierarchical": {"leaf_gene_ids": ["gene_1", "gene_2"]}},
    )


def _make_embedding_data() -> GeneEmbeddingData:
    return GeneEmbeddingData(
        gene_ids=["gene_1", "gene_2"],
        gene_names=["G1", "G2"],
        chromosomes=["chr1", "chr1"],
        starts=np.asarray([10, 60], dtype=np.int64),
        ends=np.asarray([50, 100], dtype=np.int64),
        strands=["+", "-"],
        labels=np.asarray([0, 1], dtype=np.int64),
        embedding=np.asarray([[0.1, 0.2], [0.3, 0.4]], dtype=float),
        cluster_source="kmeans",
    )


def _make_cluster_metagene_data() -> ClusterMetageneData:
    drd = DiscreteRegionData()
    drd.insert_unchecked(
        np.asarray([0.1, 0.5, 0.9], dtype=float),
        np.asarray([0.2, 0.4, 0.3], dtype=float),
        "kmeans cluster 0",
    )
    drd.insert_unchecked(
        np.asarray([0.1, 0.5, 0.9], dtype=float),
        np.asarray([0.3, 0.6, 0.5], dtype=float),
        "kmeans cluster 1",
    )
    return ClusterMetageneData(
        profiles=drd,
        segments=list(_make_layout().segments),
        groups=[
            ClusterMetageneGroup("kmeans cluster 0", ["gene_1"], 1, cluster_id=0),
            ClusterMetageneGroup("kmeans cluster 1", ["gene_2"], 1, cluster_id=1),
        ],
    )


def test_metagene_family_key_changes_with_assembly_mode() -> None:
    layout = _make_layout()
    key_annotation = metagene_family_key(
        bsx_path="report.bsx",
        annot_path="annot.gff",
        context="CG",
        assembly_mode="annotation-driven",
        layout=layout,
        limit_regions=None,
    )
    key_manual = metagene_family_key(
        bsx_path="report.bsx",
        annot_path="annot.gff",
        context="CG",
        assembly_mode="manual-composed",
        layout=layout,
        limit_regions=None,
    )
    assert key_annotation != key_manual


def test_cluster_family_keys_split_matrix_and_plot_layers() -> None:
    config_a = build_studio_cluster_config(
        bsx_path="report.bsx",
        annot_path="annot.gff",
        context="CG",
        limit_genes=100,
        min_coverage=5,
        query_block_merge_gap_bp=250,
        n_clusters=4,
        seed=0,
        hierarchical_enabled=False,
        hierarchical_max_genes=5_000,
    )
    config_b = build_studio_cluster_config(
        bsx_path="report.bsx",
        annot_path="annot.gff",
        context="CG",
        limit_genes=100,
        min_coverage=5,
        query_block_merge_gap_bp=250,
        n_clusters=6,
        seed=99,
        hierarchical_enabled=False,
        hierarchical_max_genes=5_000,
    )
    matrix_key_a = cluster_matrix_family_key(config_a)
    matrix_key_b = cluster_matrix_family_key(config_b)
    assert matrix_key_a == matrix_key_b
    assert cluster_plot_family_key(matrix_key_a, config_a) != cluster_plot_family_key(matrix_key_a, config_b)


def test_prepare_metagene_family_manual_uses_composed_builder(monkeypatch) -> None:
    drd = _make_drd()
    observed: dict[str, object] = {}

    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces._new_reader",
        lambda *args, **kwargs: object(),
    )
    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces.HcAnnotStore.from_gff",
        lambda path: object(),
    )

    def fake_collect_parts(*args, **kwargs):
        observed["parts_called"] = True
        return {"gene": ([], [])}

    def fake_build_manual(reader, *, part_map, layout):
        observed["manual_called"] = True
        observed["part_map"] = part_map
        observed["layout"] = layout
        return drd

    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces.collect_layout_parts_from_hcannot",
        fake_collect_parts,
    )
    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces.build_manual_metagene",
        fake_build_manual,
    )

    family = prepare_metagene_family(
        bsx_path="report.bsx",
        annot_path="annot.gff",
        context="CG",
        assembly_mode="manual-composed",
        layout=_make_layout(),
        limit_regions=25,
    )

    assert observed["parts_called"] is True
    assert observed["manual_called"] is True
    assert family.n_regions == 2
    assert len(family.box_segments.rows) > 0
    assert len(family.violin_segments.rows) > 0


def test_prepare_metagene_family_accepts_explicit_manual_region_specs(monkeypatch) -> None:
    drd = _make_drd()
    observed: dict[str, object] = {}

    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces._new_reader",
        lambda *args, **kwargs: object(),
    )

    def fail_collect_parts(*args, **kwargs):
        raise AssertionError("annotation collection should not run for explicit manual specs")

    def fake_manual_part_map(specs, *, layout):
        observed["manual_specs"] = specs
        return {"gene": ([object()], ["gene_1"])}

    def fake_build_manual(reader, *, part_map, layout):
        observed["manual_called"] = True
        observed["part_map"] = part_map
        return drd

    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces.collect_layout_parts_from_hcannot",
        fail_collect_parts,
    )
    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces._manual_part_map_from_specs",
        fake_manual_part_map,
    )
    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces.build_manual_metagene",
        fake_build_manual,
    )

    family = prepare_metagene_family(
        bsx_path="report.bsx",
        annot_path="annot.gff",
        context="CG",
        assembly_mode="manual-composed",
        layout=_make_layout(),
        limit_regions=None,
        manual_part_specs={
            "upstream": [("gene_1", "chr1", 10, 30, "+")],
            "body": [("gene_1", "chr1", 30, 60, "+")],
            "downstream": [("gene_1", "chr1", 60, 90, "+")],
        },
    )

    assert observed["manual_called"] is True
    assert family.n_regions == 2


def test_prepare_cluster_workspaces_reuse_matrix_and_build_plot_data(monkeypatch) -> None:
    matrix = _make_matrix()
    result = _make_cluster_result(matrix)
    embedding_data = _make_embedding_data()
    cluster_metagene_data = _make_cluster_metagene_data()
    dendrogram_data = GeneDendrogramData(
        linkage_matrix=np.asarray([[0.0, 1.0, 0.5, 2.0]], dtype=float),
        leaf_labels=["gene_1", "gene_2"],
        leaf_order=np.asarray([0, 1], dtype=np.int64),
    )

    config = build_studio_cluster_config(
        bsx_path="report.bsx",
        annot_path="annot.gff",
        context="CG",
        limit_genes=100,
        min_coverage=5,
        query_block_merge_gap_bp=250,
        n_clusters=4,
        seed=0,
        hierarchical_enabled=False,
        hierarchical_max_genes=5_000,
    )
    dendrogram_config = build_studio_cluster_config(
        bsx_path="report.bsx",
        annot_path="annot.gff",
        context="CG",
        limit_genes=100,
        min_coverage=5,
        query_block_merge_gap_bp=250,
        n_clusters=4,
        seed=0,
        hierarchical_enabled=True,
        hierarchical_max_genes=5_000,
    )

    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces.build_gene_profile_matrix",
        lambda current_config: matrix,
    )
    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces.cluster_gene_profiles",
        lambda current_matrix, current_config: result,
    )
    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces.build_gene_embedding_data",
        lambda current_result: embedding_data,
    )
    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces.build_cluster_metagene_data",
        lambda current_result: cluster_metagene_data,
    )
    monkeypatch.setattr(
        "bsx2.viz.compute.studio_workspaces.build_gene_dendrogram_data",
        lambda current_result: dendrogram_data,
    )

    matrix_workspace = prepare_cluster_matrix_workspace(config)
    cluster_family = prepare_cluster_family(matrix_workspace, config)
    dendrogram_family = prepare_cluster_dendrogram_family(matrix_workspace, dendrogram_config)

    assert matrix_workspace.feature_matrix is matrix
    assert cluster_family.feature_matrix is matrix
    assert cluster_family.embedding_data is embedding_data
    assert cluster_family.cluster_metagene_data is cluster_metagene_data
    assert dendrogram_family.dendrogram_data is dendrogram_data
