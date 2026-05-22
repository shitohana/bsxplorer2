"""Candidate gene prioritization for DMR-linked functional evidence.

This module is intentionally additive: it consumes existing DMR evidence and
external annotation/evidence tables, then produces a rule-based prioritization
table. It does not call DMRs, prove gene function, infer causality, or run
RNA-seq/chromatin processing pipelines.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable

import pandas as pd


FUNCTIONAL_PRIORITIZATION_COLUMNS = [
    "dmr_id",
    "chrom",
    "start",
    "end",
    "context",
    "region_delta",
    "region_q_value",
    "evidence_class",
    "linked_gene_id",
    "link_type",
    "distance_to_gene",
    "te_overlap",
    "te_family",
    "chromatin_overlap",
    "peak_type",
    "rna_seq_log2fc",
    "rna_seq_q_value",
    "direction_consistency",
    "functional_support_score",
    "functional_support_class",
    "interpretation",
    "missing_evidence_flags",
]


def _open_text(path: str | Path):
    table_path = Path(path)
    if table_path.suffix == ".gz":
        return gzip.open(table_path, "rt", encoding="utf-8", errors="replace")
    return table_path.open("r", encoding="utf-8", errors="replace")


def _is_number(value: object) -> bool:
    try:
        float(str(value))
        return True
    except (TypeError, ValueError):
        return False


def _find_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    normalized = {str(c).lower().replace("-", "_").replace(" ", "_"): c for c in df.columns}
    for candidate in candidates:
        key = candidate.lower()
        if key in lower:
            return lower[key]
        norm_key = key.replace("-", "_").replace(" ", "_")
        if norm_key in normalized:
            return normalized[norm_key]
    return None


def _first_data_line(path: str | Path) -> str | None:
    with _open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                return stripped
    return None


def _read_table(path: str | Path) -> pd.DataFrame:
    table_path = Path(path)
    first = _first_data_line(table_path)
    if first is None:
        return pd.DataFrame()
    delimiter = "\t" if "\t" in first else "," if "," in first else r"\s+"
    fields = first.split("\t") if delimiter == "\t" else first.split(",") if delimiter == "," else first.split()
    header = not (len(fields) >= 3 and _is_number(fields[1]) and _is_number(fields[2]))
    if header:
        return pd.read_csv(table_path, sep=None if delimiter != r"\s+" else delimiter, engine="python", comment="#")
    df = pd.read_csv(table_path, sep=delimiter, engine="python", comment="#", header=None)
    df.columns = [f"col_{i}" for i in range(df.shape[1])]
    return df


def _parse_attributes(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in str(value).strip().strip(";").split(";"):
        item = part.strip()
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
        elif " " in item:
            key, val = item.split(" ", 1)
            val = val.strip().strip('"')
        else:
            continue
        attrs[key.strip()] = val.strip().strip('"')
    return attrs


def _read_gff_like(path: str | Path, id_prefix: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, feature_type, start, end, _score, strand, _phase, attributes = fields[:9]
            if not _is_number(start) or not _is_number(end):
                continue
            attrs = _parse_attributes(attributes)
            feature = str(feature_type).lower()
            record_id = (
                attrs.get("gene_id")
                or attrs.get("ID")
                or attrs.get("Parent")
                or attrs.get("Name")
                or f"{id_prefix}_{len(rows) + 1}"
            )
            rows.append({
                f"{id_prefix}_id": record_id,
                "gene_id": record_id,
                "chrom": chrom,
                "start": max(0, int(float(start)) - 1),
                "end": int(float(end)),
                "strand": strand if strand in {"+", "-"} else ".",
                "feature_type": feature,
                "te_family": attrs.get("family") or attrs.get("Family") or attrs.get("te_family") or attrs.get("Name", "NA"),
                "te_class": attrs.get("class") or attrs.get("Class") or attrs.get("te_class") or feature_type,
            })
    return pd.DataFrame(rows)


def _normalize_interval_table(
    df: pd.DataFrame,
    *,
    id_name: str,
    default_feature_type: str,
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[id_name, "chrom", "start", "end", "strand", "feature_type"])
    chrom_col = _find_column(df, ("chrom", "chr", "chromosome", "seqname", "seqid", "col_0"))
    start_col = _find_column(df, ("start", "begin", "start_bp", "col_1"))
    end_col = _find_column(df, ("end", "stop", "end_bp", "col_2"))
    if chrom_col is None or start_col is None or end_col is None:
        raise ValueError("interval table must contain chrom/start/end columns")
    id_col = _find_column(df, (id_name, "gene_id", "id", "name", "ID", "col_3"))
    strand_col = _find_column(df, ("strand", "col_5"))
    feature_col = _find_column(df, ("feature_type", "type", "feature"))
    out = pd.DataFrame({
        id_name: df[id_col].astype(str) if id_col is not None else [f"{id_name}_{i + 1}" for i in range(len(df))],
        "chrom": df[chrom_col].astype(str),
        "start": pd.to_numeric(df[start_col], errors="coerce"),
        "end": pd.to_numeric(df[end_col], errors="coerce"),
        "strand": df[strand_col].astype(str) if strand_col is not None else ".",
        "feature_type": df[feature_col].astype(str).str.lower() if feature_col is not None else default_feature_type,
    })
    out = out.dropna(subset=["start", "end"]).copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    out.loc[out["start"] < 0, "start"] = 0
    if id_name != "gene_id" and "gene_id" in df.columns:
        out["gene_id"] = df["gene_id"].astype(str)
    return out


def read_dmr_evidence_table(path: str | Path) -> pd.DataFrame:
    """Read and normalize an existing DMR evidence table."""
    df = _read_table(path)
    if df.empty:
        return pd.DataFrame(columns=["dmr_id", "chrom", "start", "end", "context", "region_delta", "region_q_value", "evidence_class"])
    chrom_col = _find_column(df, ("chrom", "chr", "chromosome", "seqname", "seqid"))
    start_col = _find_column(df, ("start", "begin", "start_bp"))
    end_col = _find_column(df, ("end", "stop", "end_bp"))
    if chrom_col is None or start_col is None or end_col is None:
        raise ValueError("DMR evidence table must contain chrom/start/end columns")
    dmr_col = _find_column(df, ("dmr_id", "region_id", "id", "name"))
    context_col = _find_column(df, ("context", "methylation_context"))
    delta_col = _find_column(df, ("region_delta", "mean_delta", "delta", "methylation_delta"))
    q_col = _find_column(df, ("region_q_value", "q_value", "qvalue", "padj", "fdr"))
    class_col = _find_column(df, ("evidence_class", "class"))
    out = pd.DataFrame({
        "dmr_id": df[dmr_col].astype(str) if dmr_col is not None else [f"dmr_{i + 1}" for i in range(len(df))],
        "chrom": df[chrom_col].astype(str),
        "start": pd.to_numeric(df[start_col], errors="coerce"),
        "end": pd.to_numeric(df[end_col], errors="coerce"),
        "context": df[context_col].astype(str).str.upper() if context_col is not None else "NA",
        "region_delta": pd.to_numeric(df[delta_col], errors="coerce") if delta_col is not None else pd.NA,
        "region_q_value": pd.to_numeric(df[q_col], errors="coerce") if q_col is not None else pd.NA,
        "evidence_class": df[class_col].astype(str) if class_col is not None else "NA",
    })
    out = out.dropna(subset=["start", "end"]).copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    return out


def read_gene_annotation(path: str | Path) -> pd.DataFrame:
    """Read GFF/GTF/BED/TSV gene annotation into a stable interval schema."""
    suffixes = "".join(Path(path).suffixes).lower()
    if any(ext in suffixes for ext in (".gff", ".gff3", ".gtf")):
        df = _read_gff_like(path, "gene")
        if df.empty:
            return pd.DataFrame(columns=["gene_id", "chrom", "start", "end", "strand", "feature_type"])
        return df[["gene_id", "chrom", "start", "end", "strand", "feature_type"]].copy()
    return _normalize_interval_table(_read_table(path), id_name="gene_id", default_feature_type="gene")


def read_te_annotation(path: str | Path) -> pd.DataFrame:
    """Read TE annotation, if provided."""
    suffixes = "".join(Path(path).suffixes).lower()
    if any(ext in suffixes for ext in (".gff", ".gff3", ".gtf")):
        df = _read_gff_like(path, "te")
        if df.empty:
            return pd.DataFrame(columns=["te_id", "chrom", "start", "end", "te_family", "te_class"])
        return df.rename(columns={"te_id": "te_id"})[["te_id", "chrom", "start", "end", "te_family", "te_class"]].copy()
    raw = _read_table(path)
    out = _normalize_interval_table(raw, id_name="te_id", default_feature_type="TE")
    family_col = _find_column(raw, ("te_family", "family", "name", "col_3"))
    class_col = _find_column(raw, ("te_class", "class", "feature_type", "type"))
    out["te_family"] = raw[family_col].astype(str).values if family_col is not None and len(raw) == len(out) else "NA"
    out["te_class"] = raw[class_col].astype(str).values if class_col is not None and len(raw) == len(out) else "NA"
    return out[["te_id", "chrom", "start", "end", "te_family", "te_class"]]


def read_expression_table(path: str | Path) -> pd.DataFrame:
    """Read a ready-made RNA-seq differential expression table."""
    df = _read_table(path)
    if df.empty:
        return pd.DataFrame(columns=["gene_id", "log2FC", "p_value", "q_value"])
    gene_col = _find_column(df, ("gene_id", "gene", "id", "target_id"))
    logfc_col = _find_column(df, ("log2FC", "log2_fold_change", "log2FoldChange", "logFC"))
    p_col = _find_column(df, ("p_value", "pvalue", "p.val", "pval"))
    q_col = _find_column(df, ("q_value", "padj", "fdr", "adj_p_value", "adj.P.Val"))
    if gene_col is None:
        raise ValueError("expression table must contain a gene_id column")
    return pd.DataFrame({
        "gene_id": df[gene_col].astype(str),
        "log2FC": pd.to_numeric(df[logfc_col], errors="coerce") if logfc_col is not None else pd.NA,
        "p_value": pd.to_numeric(df[p_col], errors="coerce") if p_col is not None else pd.NA,
        "q_value": pd.to_numeric(df[q_col], errors="coerce") if q_col is not None else pd.NA,
    })


def read_chromatin_peaks(path: str | Path) -> pd.DataFrame:
    """Read BED/narrowPeak-like chromatin peak evidence."""
    raw = _read_table(path)
    out = _normalize_interval_table(raw, id_name="peak_id", default_feature_type="peak")
    signal_col = _find_column(raw, ("signal", "signalValue", "score", "col_4", "col_6"))
    peak_type_col = _find_column(raw, ("peak_type", "feature_type", "type"))
    out["signal"] = pd.to_numeric(raw[signal_col], errors="coerce").values if signal_col is not None and len(raw) == len(out) else pd.NA
    out["peak_type"] = raw[peak_type_col].astype(str).values if peak_type_col is not None and len(raw) == len(out) else "peak"
    return out[["peak_id", "chrom", "start", "end", "signal", "peak_type"]]


def _overlap_length(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(int(a_end), int(b_end)) - max(int(a_start), int(b_start)))


def _distance_to_interval(start: int, end: int, other_start: int, other_end: int) -> int:
    if _overlap_length(start, end, other_start, other_end) > 0:
        return 0
    if end < other_start:
        return int(other_start - end)
    return int(start - other_end)


def _gene_body_rows(gene_df: pd.DataFrame) -> pd.DataFrame:
    if gene_df.empty:
        return gene_df
    preferred = {"gene", "mrna", "transcript", "gene_body", "mrna_as_gene_body"}
    body = gene_df[gene_df["feature_type"].astype(str).str.lower().isin(preferred)].copy()
    if body.empty:
        body = gene_df.copy()
    return body


def _promoter_interval(gene: pd.Series, upstream: int, downstream: int) -> tuple[int, int]:
    strand = str(gene.get("strand", "."))
    start = int(gene["start"])
    end = int(gene["end"])
    if strand == "-":
        return max(0, end - downstream), max(0, end + upstream)
    return max(0, start - upstream), max(0, start + downstream)


def link_dmrs_to_genes(
    dmr_df: pd.DataFrame,
    gene_df: pd.DataFrame,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
    max_distance: int = 10000,
) -> pd.DataFrame:
    """Link DMRs to promoter/gene-body/proximal/intergenic gene candidates."""
    dmrs = dmr_df.copy()
    genes = _gene_body_rows(gene_df.copy())
    rows: list[dict[str, object]] = []
    genes_by_chrom = {chrom: group.reset_index(drop=True) for chrom, group in genes.groupby("chrom")} if not genes.empty else {}
    feature_by_chrom = {chrom: group.reset_index(drop=True) for chrom, group in gene_df.groupby("chrom")} if not gene_df.empty else {}

    for _, dmr in dmrs.iterrows():
        chrom = str(dmr["chrom"])
        start = int(dmr["start"])
        end = int(dmr["end"])
        base = dmr.to_dict()
        base.update({
            "linked_gene_id": "NA",
            "nearest_gene": "NA",
            "distance_to_gene": pd.NA,
            "promoter_overlap": False,
            "gene_body_overlap": False,
            "exon_overlap": False,
            "intron_overlap": False,
            "link_type": "intergenic",
        })
        chrom_genes = genes_by_chrom.get(chrom)
        if chrom_genes is None or chrom_genes.empty:
            rows.append(base)
            continue

        promoter_hits = []
        body_hits = []
        distances = []
        for _, gene in chrom_genes.iterrows():
            gene_id = str(gene["gene_id"])
            p_start, p_end = _promoter_interval(gene, promoter_upstream, promoter_downstream)
            promoter_ol = _overlap_length(start, end, p_start, p_end)
            body_ol = _overlap_length(start, end, int(gene["start"]), int(gene["end"]))
            distance = _distance_to_interval(start, end, int(gene["start"]), int(gene["end"]))
            distances.append((distance, gene_id))
            if promoter_ol > 0:
                promoter_hits.append((promoter_ol, distance, gene_id))
            if body_ol > 0:
                body_hits.append((body_ol, distance, gene_id))

        if promoter_hits:
            promoter_hits.sort(key=lambda item: (-item[0], item[1], item[2]))
            _ol, distance, gene_id = promoter_hits[0]
            base.update({"linked_gene_id": gene_id, "nearest_gene": gene_id, "distance_to_gene": distance, "promoter_overlap": True, "link_type": "promoter"})
        elif body_hits:
            body_hits.sort(key=lambda item: (-item[0], item[1], item[2]))
            _ol, distance, gene_id = body_hits[0]
            base.update({"linked_gene_id": gene_id, "nearest_gene": gene_id, "distance_to_gene": distance, "gene_body_overlap": True, "link_type": "gene_body"})
        elif distances:
            distances.sort(key=lambda item: (item[0], item[1]))
            distance, gene_id = distances[0]
            base.update({"nearest_gene": gene_id, "distance_to_gene": distance})
            if distance <= max_distance:
                base.update({"linked_gene_id": gene_id, "link_type": "proximal"})

        chrom_features = feature_by_chrom.get(chrom)
        linked_gene_id = str(base.get("linked_gene_id", "NA"))
        if chrom_features is not None and linked_gene_id != "NA":
            linked_features = chrom_features[chrom_features["gene_id"].astype(str) == linked_gene_id]
            for feature_name, column_name in (("exon", "exon_overlap"), ("intron", "intron_overlap")):
                feature_rows = linked_features[linked_features["feature_type"].astype(str).str.lower() == feature_name]
                base[column_name] = any(_overlap_length(start, end, int(row["start"]), int(row["end"])) > 0 for _, row in feature_rows.iterrows())
        rows.append(base)
    return pd.DataFrame(rows)


def _best_overlap(interval: pd.Series, candidates: pd.DataFrame, id_col: str) -> tuple[pd.Series | None, int, float]:
    best_row = None
    best_overlap = 0
    start = int(interval["start"])
    end = int(interval["end"])
    length = max(1, end - start)
    for _, candidate in candidates.iterrows():
        overlap = _overlap_length(start, end, int(candidate["start"]), int(candidate["end"]))
        if overlap > best_overlap:
            best_overlap = overlap
            best_row = candidate
    return best_row, best_overlap, best_overlap / length


def link_dmrs_to_te(dmr_df: pd.DataFrame, te_df: pd.DataFrame | None) -> pd.DataFrame:
    """Annotate DMR links to TE intervals, if TE annotation is available."""
    out = dmr_df.copy()
    out["te_input_provided"] = te_df is not None
    out["te_overlap"] = False
    out["te_family"] = "NA"
    out["te_class"] = "NA"
    out["te_overlap_fraction"] = 0.0
    if te_df is None or te_df.empty:
        return out
    te_by_chrom = {chrom: group.reset_index(drop=True) for chrom, group in te_df.groupby("chrom")}
    for index, row in out.iterrows():
        candidates = te_by_chrom.get(str(row["chrom"]))
        if candidates is None or candidates.empty:
            continue
        best, overlap, fraction = _best_overlap(row, candidates, "te_id")
        if best is not None and overlap > 0:
            out.at[index, "te_overlap"] = True
            out.at[index, "te_family"] = best.get("te_family", "NA")
            out.at[index, "te_class"] = best.get("te_class", "NA")
            out.at[index, "te_overlap_fraction"] = fraction
    return out


def link_dmrs_to_chromatin(dmr_df: pd.DataFrame, peaks_df: pd.DataFrame | None) -> pd.DataFrame:
    """Annotate DMR links to ready-made chromatin peak intervals."""
    out = dmr_df.copy()
    out["chromatin_input_provided"] = peaks_df is not None
    out["chromatin_overlap"] = False
    out["peak_type"] = "NA"
    out["peak_signal"] = pd.NA
    out["peak_overlap_fraction"] = 0.0
    if peaks_df is None or peaks_df.empty:
        return out
    peaks_by_chrom = {chrom: group.reset_index(drop=True) for chrom, group in peaks_df.groupby("chrom")}
    for index, row in out.iterrows():
        candidates = peaks_by_chrom.get(str(row["chrom"]))
        if candidates is None or candidates.empty:
            continue
        best, overlap, fraction = _best_overlap(row, candidates, "peak_id")
        if best is not None and overlap > 0:
            out.at[index, "chromatin_overlap"] = True
            out.at[index, "peak_type"] = best.get("peak_type", "peak")
            out.at[index, "peak_signal"] = best.get("signal", pd.NA)
            out.at[index, "peak_overlap_fraction"] = fraction
    return out


def _expression_direction(log2fc: object, q_value: object) -> str:
    if pd.isna(log2fc):
        return "no_expression_data"
    q_ok = pd.isna(q_value) or float(q_value) <= 0.1
    if not q_ok:
        return "no_change"
    if float(log2fc) > 0:
        return "up"
    if float(log2fc) < 0:
        return "down"
    return "no_change"


def _expression_contrast_status(
    expression_df: pd.DataFrame | None,
    dmr_contrast_label: str | None,
    expression_contrast_label: str | None,
    require_matched_contrast: bool,
) -> str:
    if expression_df is None or expression_df.empty:
        return "no_expression_data"
    if not dmr_contrast_label or not expression_contrast_label:
        return "not_validated"
    if str(dmr_contrast_label) == str(expression_contrast_label):
        return "matched"
    if require_matched_contrast:
        raise ValueError(
            "DMR and expression contrast labels do not match: "
            f"{dmr_contrast_label!r} != {expression_contrast_label!r}"
        )
    return "mismatch"


def join_expression_evidence(
    linked_df: pd.DataFrame,
    expression_df: pd.DataFrame | None,
    *,
    dmr_contrast_label: str | None = None,
    expression_contrast_label: str | None = None,
    require_matched_contrast: bool = False,
) -> pd.DataFrame:
    """Join ready-made RNA-seq differential expression evidence by linked gene."""
    out = linked_df.copy()
    contrast_status = _expression_contrast_status(
        expression_df,
        dmr_contrast_label,
        expression_contrast_label,
        require_matched_contrast,
    )
    out["rna_seq_log2fc"] = pd.NA
    out["rna_seq_q_value"] = pd.NA
    out["expression_direction"] = "no_expression_data"
    out["expression_supported"] = False
    out["expression_contrast_status"] = contrast_status
    if expression_df is None or expression_df.empty:
        return out
    expr = expression_df.rename(columns={"q_value": "rna_seq_q_value", "log2FC": "rna_seq_log2fc"})
    merged = out.merge(
        expr[["gene_id", "rna_seq_log2fc", "rna_seq_q_value"]],
        how="left",
        left_on="linked_gene_id",
        right_on="gene_id",
        suffixes=("", "_expr"),
    )
    merged["rna_seq_log2fc"] = merged["rna_seq_log2fc_expr"]
    merged["rna_seq_q_value"] = merged["rna_seq_q_value_expr"]
    merged = merged.drop(columns=[c for c in ("gene_id", "rna_seq_log2fc_expr", "rna_seq_q_value_expr") if c in merged.columns])
    merged["expression_direction"] = [
        _expression_direction(log2fc, qval)
        for log2fc, qval in zip(merged["rna_seq_log2fc"], merged["rna_seq_q_value"])
    ]
    merged["expression_supported"] = (
        pd.to_numeric(merged["rna_seq_q_value"], errors="coerce").le(0.1)
        & pd.to_numeric(merged["rna_seq_log2fc"], errors="coerce").abs().gt(0)
    ).fillna(False)
    merged["expression_contrast_status"] = contrast_status
    return merged


def _dmr_direction(delta: object) -> str:
    if pd.isna(delta):
        return "unknown"
    return "hyper" if float(delta) > 0 else "hypo" if float(delta) < 0 else "neutral"


def _direction_consistency(row: pd.Series) -> str:
    link_type = str(row.get("link_type", "intergenic"))
    direction = _dmr_direction(row.get("region_delta"))
    expr_direction = str(row.get("expression_direction", "no_expression_data"))
    contrast_status = str(row.get("expression_contrast_status", "matched"))
    if contrast_status == "not_validated":
        return "not_evaluated_contrast_not_validated"
    if contrast_status == "mismatch":
        return "not_evaluated_contrast_mismatch"
    if expr_direction == "no_expression_data":
        return "no_expression_data"
    if link_type == "promoter" and direction == "hyper" and expr_direction == "down":
        return "consistent"
    if link_type == "promoter" and direction == "hypo" and expr_direction == "up":
        return "consistent"
    if link_type == "promoter" and expr_direction in {"up", "down"}:
        return "inconsistent"
    if link_type == "gene_body":
        return "context_dependent"
    if bool(row.get("te_overlap", False)):
        return "te_context"
    return "not_applicable"


def _score_and_classify(row: pd.Series) -> tuple[int, str, str, str]:
    evidence_class = str(row.get("evidence_class", "")).lower()
    dmr_component = {"strong": 30, "moderate": 22, "weak": 12, "candidate_only": 6, "qc_limited": 0}.get(evidence_class, 8)
    link_type = str(row.get("link_type", "intergenic"))
    annotation_component = {"promoter": 20, "gene_body": 12, "proximal": 8, "intergenic": 0}.get(link_type, 0)
    expression_supported = bool(row.get("expression_supported", False))
    expr_q = row.get("rna_seq_q_value")
    expression_component = 0
    if expression_supported:
        expression_component = 15 if pd.isna(expr_q) or float(expr_q) <= 0.05 else 8
    consistency = _direction_consistency(row)
    consistency_component = {"consistent": 15, "context_dependent": 5, "te_context": 4, "inconsistent": -5}.get(consistency, 0)
    te_component = 10 if bool(row.get("te_overlap", False)) else 0
    chromatin_component = 10 if bool(row.get("chromatin_overlap", False)) else 0
    qc_penalty = 20 if evidence_class == "qc_limited" else 0
    score = int(max(0, min(100, dmr_component + annotation_component + expression_component + consistency_component + te_component + chromatin_component - qc_penalty)))

    if bool(row.get("te_overlap", False)) and score >= 35:
        support_class = "te_associated"
    elif expression_supported and dmr_component < 10:
        support_class = "expression_only"
    elif score >= 70 and consistency in {"consistent", "context_dependent"}:
        support_class = "high_confidence_candidate"
    elif score >= 50:
        support_class = "moderate_candidate"
    elif dmr_component >= 10 and not expression_supported:
        support_class = "methylation_only"
    else:
        support_class = "insufficient_support"

    direction = _dmr_direction(row.get("region_delta"))
    if link_type == "promoter" and direction == "hyper" and consistency == "consistent":
        interpretation = "Promoter hypermethylation with reduced expression candidate"
    elif link_type == "promoter" and direction == "hypo" and consistency == "consistent":
        interpretation = "Promoter hypomethylation with increased expression candidate"
    elif bool(row.get("te_overlap", False)):
        interpretation = "TE-associated methylation candidate"
    elif link_type == "gene_body":
        interpretation = "Gene-body methylation candidate; expression direction is context-dependent"
    elif bool(row.get("chromatin_overlap", False)):
        interpretation = "DMR candidate overlapping chromatin evidence"
    else:
        interpretation = "DMR-linked candidate gene prioritization result"
    return score, support_class, consistency, interpretation


def add_functional_prioritization_scores(linked_df: pd.DataFrame) -> pd.DataFrame:
    """Add rule-based functional support scores and conservative classes."""
    out = linked_df.copy()
    scores = [_score_and_classify(row) for _, row in out.iterrows()]
    if scores:
        out["functional_support_score"] = [item[0] for item in scores]
        out["functional_support_class"] = [item[1] for item in scores]
        out["direction_consistency"] = [item[2] for item in scores]
        out["interpretation"] = [item[3] for item in scores]
    else:
        out["functional_support_score"] = pd.Series(dtype="int64")
        out["functional_support_class"] = pd.Series(dtype="object")
        out["direction_consistency"] = pd.Series(dtype="object")
        out["interpretation"] = pd.Series(dtype="object")

    flags: list[str] = []
    for _, row in out.iterrows():
        row_flags: list[str] = []
        if str(row.get("linked_gene_id", "NA")) == "NA":
            row_flags.append("no_gene_link")
        if str(row.get("expression_direction", "no_expression_data")) == "no_expression_data":
            row_flags.append("no_expression_data")
        contrast_status = str(row.get("expression_contrast_status", "matched"))
        if contrast_status == "not_validated":
            row_flags.append("expression_contrast_not_validated")
        elif contrast_status == "mismatch":
            row_flags.append("expression_contrast_mismatch")
        if not bool(row.get("te_input_provided", False)):
            row_flags.append("te_input_missing")
        elif not bool(row.get("te_overlap", False)):
            row_flags.append("te_no_overlap")
        if not bool(row.get("chromatin_input_provided", False)):
            row_flags.append("chromatin_input_missing")
        elif not bool(row.get("chromatin_overlap", False)):
            row_flags.append("chromatin_no_overlap")
        flags.append(";".join(row_flags) if row_flags else "none")
    out["missing_evidence_flags"] = flags

    for column in FUNCTIONAL_PRIORITIZATION_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    return out[FUNCTIONAL_PRIORITIZATION_COLUMNS]


def prioritize_dmr_linked_genes(
    dmr_df: pd.DataFrame,
    gene_df: pd.DataFrame,
    te_df: pd.DataFrame | None = None,
    expression_df: pd.DataFrame | None = None,
    chromatin_df: pd.DataFrame | None = None,
    *,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
    max_distance: int = 10000,
    dmr_contrast_label: str | None = None,
    expression_contrast_label: str | None = None,
    require_matched_contrast: bool = False,
) -> pd.DataFrame:
    """Run the full candidate gene prioritization workflow from normalized tables."""
    linked = link_dmrs_to_genes(dmr_df, gene_df, promoter_upstream, promoter_downstream, max_distance)
    linked = link_dmrs_to_te(linked, te_df)
    linked = link_dmrs_to_chromatin(linked, chromatin_df)
    linked = join_expression_evidence(
        linked,
        expression_df,
        dmr_contrast_label=dmr_contrast_label,
        expression_contrast_label=expression_contrast_label,
        require_matched_contrast=require_matched_contrast,
    )
    return add_functional_prioritization_scores(linked)


def write_dmr_functional_prioritization_table(df: pd.DataFrame, out_path: str | Path) -> None:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
