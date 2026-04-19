from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np

import bsx2.clustering.gene_profile as gene_profile
from bsx2.clustering.config import (
    BlockCacheConfig,
    BlockCacheMode,
    ClusterConfig,
    GeneProfileConfig,
    OutputConfig,
    ReadConfig,
)
from bsx2.clustering.gene_profile import (
    build_gene_profile_matrix,
    build_metagene_bins,
    clear_region_reader_cache,
)
from bsx2.clustering.models import GeneAnnotation


def _config() -> ClusterConfig:
    return ClusterConfig(
        bsx_path=Path("sample.bsx"),
        annotation_path=Path("genes.gff"),
        gene_profile=GeneProfileConfig(
            upstream_bp=10,
            downstream_bp=10,
            upstream_bins=1,
            body_bins=2,
            downstream_bins=1,
            max_gene_missing_rate=1.0,
            max_feature_missing_rate=1.0,
        ),
        output=OutputConfig(output_dir=Path(".")),
    )


def test_build_metagene_bins_is_strand_aware() -> None:
    config = _config()
    plus_bins = build_metagene_bins(
        GeneAnnotation("gene_plus", "chr1", 10, 50, "+"),
        config,
    )
    minus_bins = build_metagene_bins(
        GeneAnnotation("gene_minus", "chr1", 10, 50, "-"),
        config,
    )

    assert [(bin.segment, bin.start, bin.end) for bin in plus_bins] == [
        ("up", 0, 10),
        ("body", 10, 30),
        ("body", 30, 50),
        ("down", 50, 60),
    ]
    assert [(bin.segment, bin.start, bin.end) for bin in minus_bins] == [
        ("up", 50, 60),
        ("body", 30, 50),
        ("body", 10, 30),
        ("down", 0, 10),
    ]


def test_build_metagene_interval_arrays_match_public_bins() -> None:
    config = _config()
    for gene in (
        GeneAnnotation("gene_plus", "chr1", 10, 50, "+"),
        GeneAnnotation("gene_minus", "chr1", 10, 50, "-"),
    ):
        starts, ends = gene_profile._build_metagene_interval_arrays(gene, config)
        public_bins = build_metagene_bins(gene, config)
        assert list(zip(starts.tolist(), ends.tolist(), strict=True)) == [
            (profile_bin.start, profile_bin.end) for profile_bin in public_bins
        ]


def test_build_metagene_interval_matrices_match_per_gene_arrays() -> None:
    config = _config()
    genes = [
        GeneAnnotation("gene_plus", "chr1", 10, 50, "+"),
        GeneAnnotation("gene_minus", "chr1", 10, 50, "-"),
        GeneAnnotation("gene_plus_2", "chr1", 25, 75, "+"),
    ]

    starts_matrix, ends_matrix = gene_profile._build_metagene_interval_matrices(genes, config)

    for row_idx, gene in enumerate(genes):
        starts, ends = gene_profile._build_metagene_interval_arrays(gene, config)
        assert np.array_equal(starts_matrix[row_idx], starts)
        assert np.array_equal(ends_matrix[row_idx], ends)


def test_build_query_blocks_merges_overlapping_gene_spans() -> None:
    config = _config()
    blocks = gene_profile._build_query_blocks(
        [
            GeneAnnotation("gene_a", "chr1", 10, 50, "+"),
            GeneAnnotation("gene_b", "chr1", 45, 80, "+"),
            GeneAnnotation("gene_c", "chr1", 200, 230, "+"),
        ],
        config,
    )

    assert [(block.chrom, block.start, block.end) for block in blocks] == [
        ("chr1", 0, 90),
        ("chr1", 190, 240),
    ]


def test_build_query_blocks_can_merge_small_gaps_when_configured() -> None:
    config = ClusterConfig(
        bsx_path=Path("sample.bsx"),
        annotation_path=Path("genes.gff"),
        read=ReadConfig(query_block_merge_gap_bp=35),
        gene_profile=_config().gene_profile,
        output=OutputConfig(output_dir=Path(".")),
    )
    blocks = gene_profile._build_query_blocks(
        [
            GeneAnnotation("gene_a", "chr1", 10, 50, "+"),
            GeneAnnotation("gene_b", "chr1", 105, 120, "+"),
            GeneAnnotation("gene_c", "chr1", 200, 230, "+"),
        ],
        config,
    )

    assert [(block.chrom, block.start, block.end) for block in blocks] == [
        ("chr1", 0, 130),
        ("chr1", 190, 240),
    ]


def test_assign_genes_to_query_blocks_matches_expected_blocks() -> None:
    config = _config()
    genes = [
        GeneAnnotation("gene_a", "chr1", 10, 50, "+"),
        GeneAnnotation("gene_b", "chr1", 45, 80, "+"),
        GeneAnnotation("gene_c", "chr1", 200, 230, "+"),
        GeneAnnotation("gene_d", "chr2", 25, 60, "-"),
    ]
    spans = gene_profile._collect_gene_profile_spans(genes, config)
    blocks = gene_profile._build_query_blocks(genes, config, gene_spans=spans)
    assigned = gene_profile._assign_genes_to_query_blocks(
        genes,
        blocks,
        config,
        gene_spans=spans,
    )

    assert [(item.block.chrom, item.block.start, item.block.end) for item in assigned] == [
        ("chr1", 0, 90),
        ("chr1", 190, 240),
        ("chr2", 15, 70),
    ]
    assert [item.gene_indices for item in assigned] == [
        (0, 1),
        (2,),
        (3,),
    ]


def test_build_gene_profile_matrix_uses_coverage_weighted_bins(monkeypatch) -> None:
    config = _config()

    class _FakeRegionReader:
        def chr_order(self):
            return ["chr1"]

        def reset(self):
            return None

    monkeypatch.setattr(
        gene_profile,
        "load_gene_annotations",
        lambda *args, **kwargs: [
            GeneAnnotation("gene_plus", "chr1", 10, 50, "+"),
            GeneAnnotation("gene_minus", "chr1", 10, 50, "-"),
        ],
    )
    monkeypatch.setattr(
        gene_profile,
        "_make_region_reader",
        lambda *args, **kwargs: gene_profile._PreparedRegionReader(
            reader=_FakeRegionReader(),
            chr_order=("chr1",),
            cache_hit=False,
        ),
    )
    monkeypatch.setattr(
        gene_profile,
        "_query_block_arrays",
        lambda *args, **kwargs: (
            np.array([5, 15, 25, 35, 45, 55], dtype=np.int64),
            np.array([1, 2, 3, 8, 9, 10], dtype=np.int64),
            np.array([2, 4, 6, 10, 12, 20], dtype=np.int64),
        ),
    )

    matrix = build_gene_profile_matrix(config)

    assert matrix.gene_ids == ["gene_plus", "gene_minus"]
    assert [feature_bin.feature_name for feature_bin in matrix.feature_bins] == [
        "up_1",
        "body_1",
        "body_2",
        "down_1",
    ]
    assert np.allclose(
        matrix.values,
        np.array(
            [
                [0.5, 0.5, 17 / 22, 0.5],
                [0.5, 17 / 22, 0.5, 0.5],
            ],
            dtype=float,
        ),
    )


def test_profile_genes_in_block_matches_per_gene_path() -> None:
    config = _config()
    genes = [
        GeneAnnotation("gene_plus", "chr1", 10, 50, "+"),
        GeneAnnotation("gene_minus", "chr1", 10, 50, "-"),
        GeneAnnotation("gene_plus_2", "chr1", 25, 75, "+"),
    ]
    block_cumsums = gene_profile._prepare_block_cumsums(
        np.array([5, 15, 25, 35, 45, 55, 65, 75, 85], dtype=np.int64),
        np.array([1, 2, 3, 8, 9, 10, 6, 4, 2], dtype=np.int64),
        np.array([2, 4, 6, 10, 12, 20, 12, 8, 4], dtype=np.int64),
    )

    batch_profiles = gene_profile._profile_genes_in_block(
        genes,
        (0, 1, 2),
        config,
        block_cumsums,
    )
    stacked_profiles = np.vstack(
        [gene_profile._profile_gene(gene, config, block_cumsums) for gene in genes]
    )

    assert np.allclose(batch_profiles, stacked_profiles, equal_nan=True)


def test_make_region_reader_reuses_python_cache(monkeypatch) -> None:
    clear_region_reader_cache()
    bsx_path = Path("python/tests/clustering/_cache_test_sample.bsx")
    bsx_path.write_bytes(b"bsx")
    try:
        config = ClusterConfig(
            bsx_path=bsx_path,
            annotation_path=Path("genes.gff"),
            output=OutputConfig(output_dir=Path(".")),
        )
        events: list[str] = []

        class _FakeRegionReader:
            def __init__(self, path: str):
                assert path == str(bsx_path)
                events.append("init")

            def filter_coverage_gt(self, value: int):
                events.append(f"coverage:{value}")

            def filter_context(self, value):
                events.append(f"context:{value}")

            def chr_order(self):
                events.append("chr_order")
                return ["chr1", "chr2"]

            def reset(self):
                events.append("reset")

            def clear_filters(self):
                events.append("clear_filters")

        monkeypatch.setattr(gene_profile, "RegionReader", _FakeRegionReader)

        first = gene_profile._make_region_reader(config)
        second = gene_profile._make_region_reader(config)

        assert first.cache_hit is False
        assert second.cache_hit is True
        assert first.reader is second.reader
        assert first.chr_order == ("chr1", "chr2")
        assert second.chr_order == ("chr1", "chr2")
        assert events.count("init") == 1
        assert events.count("chr_order") == 1
        assert events.count("reset") == 1
        assert events.count("clear_filters") == 1
    finally:
        clear_region_reader_cache()
        bsx_path.unlink(missing_ok=True)


def test_build_gene_profile_matrix_reuses_persistent_block_cache(monkeypatch) -> None:
    clear_region_reader_cache()
    cache_dir = Path("python/tests/clustering/_query_block_cache")
    shutil.rmtree(cache_dir, ignore_errors=True)
    config = ClusterConfig(
        bsx_path=Path("sample.bsx"),
        annotation_path=Path("genes.gff"),
        block_cache=BlockCacheConfig(enabled=True, cache_dir=cache_dir),
        gene_profile=GeneProfileConfig(
            upstream_bp=10,
            downstream_bp=10,
            upstream_bins=1,
            body_bins=2,
            downstream_bins=1,
            max_gene_missing_rate=1.0,
            max_feature_missing_rate=1.0,
        ),
        output=OutputConfig(output_dir=Path(".")),
    )
    query_calls = 0

    class _FakeRegionReader:
        def reset(self):
            return None

    monkeypatch.setattr(
        gene_profile,
        "load_gene_annotations",
        lambda *args, **kwargs: [GeneAnnotation("gene_plus", "chr1", 10, 50, "+")],
    )
    monkeypatch.setattr(
        gene_profile,
        "_make_region_reader",
        lambda *args, **kwargs: gene_profile._PreparedRegionReader(
            reader=_FakeRegionReader(),
            chr_order=("chr1",),
            cache_hit=False,
        ),
    )

    def _fake_query(*args, **kwargs):
        nonlocal query_calls
        query_calls += 1
        return (
            np.array([5, 15, 25, 35, 45, 55], dtype=np.int64),
            np.array([1, 2, 3, 8, 9, 10], dtype=np.int64),
            np.array([2, 4, 6, 10, 12, 20], dtype=np.int64),
        )

    monkeypatch.setattr(gene_profile, "_query_block_arrays", _fake_query)

    try:
        first = build_gene_profile_matrix(config)
        clear_region_reader_cache()
        second = build_gene_profile_matrix(config)

        assert query_calls == 1
        assert first.metadata["query_block_cache"]["hits"] == 0
        assert first.metadata["query_block_cache"]["misses"] == 1
        assert first.metadata["query_block_cache"]["writes"] == 1
        assert second.metadata["query_block_cache"]["hits"] == 1
        assert second.metadata["query_block_cache"]["misses"] == 0
        assert second.metadata["query_block_cache"]["writes"] == 0
        assert np.allclose(first.values, second.values)
    finally:
        clear_region_reader_cache()
        shutil.rmtree(cache_dir, ignore_errors=True)


def test_build_gene_profile_matrix_supports_uncompressed_block_cache(monkeypatch) -> None:
    clear_region_reader_cache()
    cache_dir = Path("python/tests/clustering/_query_block_cache_uncompressed")
    shutil.rmtree(cache_dir, ignore_errors=True)
    config = ClusterConfig(
        bsx_path=Path("sample.bsx"),
        annotation_path=Path("genes.gff"),
        block_cache=BlockCacheConfig(
            enabled=True,
            cache_dir=cache_dir,
            mode=BlockCacheMode.UNCOMPRESSED,
        ),
        gene_profile=GeneProfileConfig(
            upstream_bp=10,
            downstream_bp=10,
            upstream_bins=1,
            body_bins=2,
            downstream_bins=1,
            max_gene_missing_rate=1.0,
            max_feature_missing_rate=1.0,
        ),
        output=OutputConfig(output_dir=Path(".")),
    )
    query_calls = 0

    class _FakeRegionReader:
        def reset(self):
            return None

    monkeypatch.setattr(
        gene_profile,
        "load_gene_annotations",
        lambda *args, **kwargs: [GeneAnnotation("gene_plus", "chr1", 10, 50, "+")],
    )
    monkeypatch.setattr(
        gene_profile,
        "_make_region_reader",
        lambda *args, **kwargs: gene_profile._PreparedRegionReader(
            reader=_FakeRegionReader(),
            chr_order=("chr1",),
            cache_hit=False,
        ),
    )

    def _fake_query(*args, **kwargs):
        nonlocal query_calls
        query_calls += 1
        return (
            np.array([5, 15, 25, 35, 45, 55], dtype=np.int64),
            np.array([1, 2, 3, 8, 9, 10], dtype=np.int64),
            np.array([2, 4, 6, 10, 12, 20], dtype=np.int64),
        )

    monkeypatch.setattr(gene_profile, "_query_block_arrays", _fake_query)

    try:
        first = build_gene_profile_matrix(config)
        clear_region_reader_cache()
        second = build_gene_profile_matrix(config)

        cache_paths = list(cache_dir.glob("*.npz"))
        assert len(cache_paths) == 1
        assert query_calls == 1
        assert first.metadata["query_block_cache"]["mode"] == "uncompressed"
        assert second.metadata["query_block_cache"]["hits"] == 1
        assert np.allclose(first.values, second.values)
    finally:
        clear_region_reader_cache()
        shutil.rmtree(cache_dir, ignore_errors=True)
