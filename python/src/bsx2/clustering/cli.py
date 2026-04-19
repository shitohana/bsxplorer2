from __future__ import annotations

import argparse
from pathlib import Path
from time import perf_counter

from bsx2 import Context

from .backend import cluster_gene_profiles
from .config import (
    AnnotationFormat,
    BackendConfig,
    BlockCacheConfig,
    BlockCacheMode,
    ClusterConfig,
    ClusterSource,
    GeneProfileConfig,
    HierarchicalConfig,
    HierarchicalDistance,
    HierarchicalLinkage,
    HierarchicalMode,
    NormalizationMode,
    OutputConfig,
    PcaSolver,
    ReadConfig,
    SilhouetteConfig,
    SilhouetteMode,
    TableFormat,
)
from .gene_profile import build_gene_profile_matrix
from .io import write_cluster_outputs


def _parse_context(value: str | None) -> Context | None:
    if value is None:
        return None
    try:
        return getattr(Context, value.upper())
    except AttributeError as exc:
        raise ValueError(f"Unsupported methylation context: {value}") from exc


def _parse_annotation_format(value: str | None) -> AnnotationFormat | None:
    return None if value is None else AnnotationFormat(value)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cluster-bsx",
        description=(
            "Cluster genes from a single BSX sample using metagene profiles, PCA, "
            "KMeans, and hierarchical clustering."
        ),
    )
    parser.add_argument("--bsx", required=True, help="Input .bsx path")
    parser.add_argument("--genes", required=True, help="Gene annotation path (GFF/GTF/BED)")
    parser.add_argument(
        "--annotation-format",
        choices=[fmt.value for fmt in AnnotationFormat],
        default=None,
        help="Optional explicit annotation format",
    )
    parser.add_argument(
        "--context",
        choices=["CG", "CHG", "CHH"],
        default=None,
        help="Optional methylation context filter",
    )
    parser.add_argument(
        "--min-coverage",
        type=int,
        default=5,
        help="Minimum per-cytosine coverage filter",
    )
    parser.add_argument(
        "--query-block-cache",
        action="store_true",
        help="Persist queried BSX blocks on disk and reuse them across runs",
    )
    parser.add_argument(
        "--query-block-cache-dir",
        default=None,
        help="Optional directory for persistent queried-block cache",
    )
    parser.add_argument(
        "--query-block-cache-mode",
        choices=[mode.value for mode in BlockCacheMode],
        default=BlockCacheMode.COMPRESSED.value,
        help="Persistent queried-block cache mode: compressed or uncompressed",
    )
    parser.add_argument(
        "--query-block-merge-gap-bp",
        type=int,
        default=0,
        help="Merge nearby gene profile spans into a shared query block when gaps are below this size",
    )
    parser.add_argument("--upstream-bp", type=int, default=2000)
    parser.add_argument("--downstream-bp", type=int, default=2000)
    parser.add_argument("--upstream-bins", type=int, default=20)
    parser.add_argument("--body-bins", type=int, default=50)
    parser.add_argument("--downstream-bins", type=int, default=20)
    parser.add_argument("--min-gene-length-bp", type=int, default=0)
    parser.add_argument("--min-bin-total-coverage", type=int, default=1)
    parser.add_argument("--max-gene-missing-rate", type=float, default=0.4)
    parser.add_argument("--max-feature-missing-rate", type=float, default=0.95)
    parser.add_argument("--min-gene-profile-variance", type=float, default=0.0)
    parser.add_argument("--min-feature-variance", type=float, default=0.0)
    parser.add_argument(
        "--normalization",
        choices=[mode.value for mode in NormalizationMode],
        default=NormalizationMode.ROW_ZSCORE.value,
    )
    parser.add_argument("--limit-genes", type=int, default=None)
    parser.add_argument("--n-components", type=int, default=2)
    parser.add_argument("--n-clusters", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--pca-solver",
        choices=[solver.value for solver in PcaSolver],
        default=PcaSolver.AUTO.value,
        help="PCA backend: exact, truncated, or auto",
    )
    parser.add_argument(
        "--pca-exact-max-matrix-size",
        type=int,
        default=1_000_000,
        help=(
            "Use exact PCA at or below this matrix size (n_genes * n_bins) in auto mode; "
            "use 0 to disable the threshold"
        ),
    )
    parser.add_argument(
        "--no-silhouette",
        action="store_true",
        help="Disable silhouette scoring entirely",
    )
    parser.add_argument(
        "--silhouette-mode",
        choices=[mode.value for mode in SilhouetteMode],
        default=SilhouetteMode.AUTO.value,
        help="Silhouette execution policy: exact, sampled, or auto",
    )
    parser.add_argument(
        "--silhouette-max-samples",
        type=int,
        default=2_000,
        help="Maximum genes to sample for silhouette in sampled mode; use 0 to disable",
    )
    parser.add_argument(
        "--silhouette-exact-max-genes",
        type=int,
        default=2_000,
        help="Use exact silhouette at or below this retained gene count; use 0 to disable",
    )
    parser.add_argument(
        "--distance",
        choices=[distance.value for distance in HierarchicalDistance],
        default=HierarchicalDistance.CORRELATION.value,
    )
    parser.add_argument(
        "--linkage",
        choices=[linkage.value for linkage in HierarchicalLinkage],
        default=HierarchicalLinkage.AVERAGE.value,
    )
    parser.add_argument(
        "--no-hierarchical",
        action="store_true",
        help="Disable hierarchical clustering and skip dendrogram/linkage outputs",
    )
    parser.add_argument(
        "--hierarchical-mode",
        choices=[mode.value for mode in HierarchicalMode],
        default=HierarchicalMode.AUTO.value,
        help="Hierarchical execution policy: exact, subsample, skip, or auto",
    )
    parser.add_argument(
        "--hierarchical-max-genes",
        type=int,
        default=5_000,
        help=(
            "Skip hierarchical clustering when retained gene count exceeds this threshold; "
            "use 0 to disable the threshold"
        ),
    )
    parser.add_argument(
        "--hierarchical-subsample-genes",
        type=int,
        default=2_000,
        help=(
            "Maximum genes to retain in subsampled hierarchical mode; use 0 to disable "
            "subsampling fallback"
        ),
    )
    parser.add_argument(
        "--cluster-source",
        choices=[source.value for source in ClusterSource],
        default=ClusterSource.KMEANS.value,
    )
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--prefix", default="", help="Prefix for output files")
    parser.add_argument(
        "--table-format",
        action="append",
        choices=[fmt.value for fmt in TableFormat],
        default=None,
        help="Repeat to write multiple table formats; default is tsv",
    )
    parser.add_argument(
        "--no-table-files",
        action="store_true",
        help="Do not write tabular outputs to disk",
    )
    parser.add_argument(
        "--no-metrics-file",
        action="store_true",
        help="Do not write metrics.json to disk",
    )
    return parser


def parse_args(argv: list[str] | None = None) -> ClusterConfig:
    args = build_parser().parse_args(argv)
    return ClusterConfig(
        bsx_path=Path(args.bsx),
        annotation_path=Path(args.genes),
        annotation_format=_parse_annotation_format(args.annotation_format),
        read=ReadConfig(
            context=_parse_context(args.context),
            min_coverage=args.min_coverage,
            query_block_merge_gap_bp=max(args.query_block_merge_gap_bp, 0),
        ),
        block_cache=BlockCacheConfig(
            enabled=bool(args.query_block_cache or args.query_block_cache_dir),
            cache_dir=None if args.query_block_cache_dir is None else Path(args.query_block_cache_dir),
            mode=BlockCacheMode(args.query_block_cache_mode),
        ),
        gene_profile=GeneProfileConfig(
            upstream_bp=args.upstream_bp,
            downstream_bp=args.downstream_bp,
            upstream_bins=args.upstream_bins,
            body_bins=args.body_bins,
            downstream_bins=args.downstream_bins,
            min_gene_length_bp=args.min_gene_length_bp,
            min_bin_total_coverage=args.min_bin_total_coverage,
            max_gene_missing_rate=args.max_gene_missing_rate,
            max_feature_missing_rate=args.max_feature_missing_rate,
            min_gene_profile_variance=args.min_gene_profile_variance,
            min_feature_variance=args.min_feature_variance,
            normalization=NormalizationMode(args.normalization),
            limit_genes=args.limit_genes,
        ),
        backend=BackendConfig(
            n_components=args.n_components,
            n_clusters=args.n_clusters,
            seed=args.seed,
            pca_solver=PcaSolver(args.pca_solver),
            pca_exact_max_matrix_size=(
                None
                if args.pca_exact_max_matrix_size is None or args.pca_exact_max_matrix_size <= 0
                else args.pca_exact_max_matrix_size
            ),
        ),
        silhouette=SilhouetteConfig(
            enabled=not args.no_silhouette,
            mode=SilhouetteMode(args.silhouette_mode),
            max_samples=(
                None
                if args.silhouette_max_samples is None or args.silhouette_max_samples <= 0
                else args.silhouette_max_samples
            ),
            exact_max_genes=(
                None
                if args.silhouette_exact_max_genes is None
                or args.silhouette_exact_max_genes <= 0
                else args.silhouette_exact_max_genes
            ),
            seed=args.seed,
        ),
        hierarchical=HierarchicalConfig(
            enabled=not args.no_hierarchical,
            max_genes=(
                None
                if args.hierarchical_max_genes is None or args.hierarchical_max_genes <= 0
                else args.hierarchical_max_genes
            ),
            mode=HierarchicalMode(args.hierarchical_mode),
            subsample_genes=(
                None
                if args.hierarchical_subsample_genes is None
                or args.hierarchical_subsample_genes <= 0
                else args.hierarchical_subsample_genes
            ),
            distance=HierarchicalDistance(args.distance),
            linkage=HierarchicalLinkage(args.linkage),
        ),
        cluster_source=ClusterSource(args.cluster_source),
        output=OutputConfig(
            output_dir=Path(args.output),
            prefix=args.prefix,
            table_formats=tuple(
                TableFormat(fmt) for fmt in (args.table_format or [TableFormat.TSV.value])
            ),
            write_table_files=not args.no_table_files,
            write_metrics_file=not args.no_metrics_file,
        ),
    )


def main(argv: list[str] | None = None) -> int:
    t0 = perf_counter()
    config = parse_args(argv)
    t_parse = perf_counter()
    matrix = build_gene_profile_matrix(config)
    t_build = perf_counter()
    result = cluster_gene_profiles(matrix, config)
    t_cluster = perf_counter()
    result.metadata["pipeline_timings_s"] = {
        "parse_args": round(t_parse - t0, 6),
        "build_gene_profile_matrix": round(t_build - t_parse, 6),
        "cluster_gene_profiles": round(t_cluster - t_build, 6),
        "pre_write_total": round(t_cluster - t0, 6),
    }
    write_cluster_outputs(result, config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
