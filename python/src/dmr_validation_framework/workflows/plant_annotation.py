#!/usr/bin/env python
"""Download and normalize plant annotation inputs for DMR downstream analysis."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import urllib.request
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd


SGN_ITAG32_BASE_URL = "https://solgenomics.net/ftp/tomato_genome/annotation/annotation/ITAG3.2_release/"


@dataclass(frozen=True)
class AnnotationProfile:
    name: str
    species: str
    assembly: str
    annotation_version: str
    base_url: str
    gene_gff: str | None
    te_gff: str | None = None
    go_table: str | None = None
    source_note: str = ""
    output_prefix: str = "plant_annotation"


PROFILES = {
    "custom": AnnotationProfile(
        name="custom",
        species="custom",
        assembly="custom",
        annotation_version="custom",
        base_url="",
        gene_gff=None,
        te_gff=None,
        go_table=None,
        source_note="User-provided local annotation inputs.",
        output_prefix="custom",
    ),
    "tomato_itag3_2": AnnotationProfile(
        name="tomato_itag3_2",
        species="Solanum lycopersicum",
        assembly="SL3.0",
        annotation_version="ITAG3.2",
        base_url=SGN_ITAG32_BASE_URL,
        gene_gff="ITAG3.2_gene_models.gff",
        te_gff="ITAG3.2_RepeatModeler_repeats_light.gff",
        go_table="ITAG3.2_protein_go.tsv",
        source_note="SGN ITAG3.2 annotation release.",
        output_prefix="tomato_itag3_2",
    ),
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare a common annotation package for DMR validation: normalized gene GFF3, "
            "metagene BED, optional TE BED, optional GO table, QC summaries, and manifests."
        )
    )
    parser.add_argument("--root", help="Dataset root. Defaults out-dir to <root>/annotation_<profile> when out-dir is omitted.")
    parser.add_argument("--out-dir", help="Output annotation directory.")
    parser.add_argument("--profile", choices=sorted(PROFILES), default="tomato_itag3_2")
    parser.add_argument("--species", help="Override species metadata in the output manifest.")
    parser.add_argument("--assembly", help="Override assembly metadata in the output manifest.")
    parser.add_argument("--annotation-version", help="Override annotation version metadata in the output manifest.")
    parser.add_argument("--output-prefix", help="Override normalized gene GFF3 filename prefix.")
    parser.add_argument("--gene-gff", help="Local gene GFF/GFF3 path. Overrides profile download for genes.")
    parser.add_argument("--te-gff", help="Local TE/repeat GFF/GFF3 path. Overrides profile download for repeats.")
    parser.add_argument("--go-table", help="Local GO TSV path. Overrides profile download for GO.")
    parser.add_argument("--base-url", help="Override profile base URL for downloading profile filenames.")
    parser.add_argument("--download-dir", help="Defaults to <out-dir>/downloads.")
    parser.add_argument("--skip-download", action="store_true", help="Use existing local/downloaded files only.")
    parser.add_argument("--force-download", action="store_true", help="Re-download files even if they already exist.")
    parser.add_argument("--include-te", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--include-go", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--gene-feature-types", default="gene", help="Comma-separated gene feature types to export.")
    parser.add_argument("--te-feature-types", default="", help="Optional comma-separated repeat feature types. Empty keeps all GFF rows.")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def _default_out_dir(args: argparse.Namespace, profile: AnnotationProfile) -> Path:
    if args.out_dir:
        return Path(args.out_dir)
    if args.root:
        return Path(args.root) / f"annotation_{profile.annotation_version.lower().replace('.', '_')}"
    return Path("outputs") / f"annotation_{profile.name}"


def _split_csv(value: str) -> set[str]:
    return {item.strip().lower() for item in value.split(",") if item.strip()}


def _parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in str(raw).strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def _clean_feature_id(value: str | None, *, prefixes: tuple[str, ...] = ()) -> str:
    if value is None:
        return ""
    cleaned = str(value).strip()
    for prefix in prefixes:
        if cleaned.startswith(prefix):
            cleaned = cleaned[len(prefix) :]
    return cleaned


def _download(url: str, target: Path, *, force: bool = False) -> dict[str, object]:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not force:
        return {
            "source": url,
            "path": str(target),
            "status": "existing",
            "bytes": target.stat().st_size,
        }
    tmp = target.with_suffix(target.suffix + ".tmp")
    with urllib.request.urlopen(url, timeout=120) as response, tmp.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    tmp.replace(target)
    return {
        "source": url,
        "path": str(target),
        "status": "downloaded",
        "bytes": target.stat().st_size,
    }


def _resolve_input(
    *,
    local_path: str | None,
    filename: str | None,
    base_url: str,
    download_dir: Path,
    skip_download: bool,
    force_download: bool,
) -> tuple[Path | None, dict[str, object] | None]:
    if local_path:
        path = Path(local_path)
        if not path.exists():
            raise FileNotFoundError(path)
        return path, {"source": str(path), "path": str(path), "status": "local", "bytes": path.stat().st_size}
    if not filename:
        return None, None
    target = download_dir / filename
    if skip_download:
        if not target.exists():
            raise FileNotFoundError(f"{target} does not exist and --skip-download was used")
        return target, {"source": str(target), "path": str(target), "status": "existing", "bytes": target.stat().st_size}
    row = _download(urljoin(base_url, filename), target, force=force_download)
    return target, row


def _iter_gff_rows(path: Path):
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature_type, start, end, score, strand, phase, attrs = fields[:9]
            try:
                start_i = int(float(start))
                end_i = int(float(end))
            except ValueError:
                continue
            yield {
                "chrom": chrom,
                "source": source,
                "feature_type": feature_type,
                "start": start_i,
                "end": end_i,
                "score": score,
                "strand": strand if strand in {"+", "-"} else ".",
                "phase": phase,
                "attributes": attrs,
                "attrs": _parse_attributes(attrs),
            }


def normalize_gene_annotation(
    source: Path,
    *,
    out_dir: Path,
    profile: AnnotationProfile,
    feature_types: set[str],
) -> tuple[pd.DataFrame, dict[str, Path]]:
    rows: list[dict[str, object]] = []
    feature_counts: dict[str, int] = {}
    seq_counts: dict[str, int] = {}
    for row in _iter_gff_rows(source):
        feature = str(row["feature_type"]).lower()
        feature_counts[feature] = feature_counts.get(feature, 0) + 1
        seq_counts[str(row["chrom"])] = seq_counts.get(str(row["chrom"]), 0) + 1
        if feature not in feature_types:
            continue
        attrs = row["attrs"]
        gene_id = (
            _clean_feature_id(attrs.get("ID"), prefixes=("gene:", "Gene:"))
            or _clean_feature_id(attrs.get("gene_id"), prefixes=("gene:", "Gene:"))
            or _clean_feature_id(attrs.get("Name"), prefixes=("gene:", "Gene:"))
            or f"gene_{len(rows) + 1}"
        )
        alias = _clean_feature_id(attrs.get("Alias"))
        rows.append(
            {
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "gene_id": gene_id,
                "gene_alias": alias,
                "strand": row["strand"],
                "feature_type": feature,
                "source": row["source"],
                "length_bp": max(0, int(row["end"]) - int(row["start"]) + 1),
            }
        )
    genes = pd.DataFrame(rows)
    if genes.empty:
        raise ValueError(f"No gene rows found in {source} for feature types {sorted(feature_types)}")
    genes = genes.sort_values(["chrom", "start", "end", "gene_id"], kind="mergesort").reset_index(drop=True)
    genes["gene_index"] = range(len(genes))

    normalized_gff = out_dir / f"{profile.output_prefix}_genes.normalized.gff3"
    with normalized_gff.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for rec in genes.itertuples(index=False):
            attrs = f"ID={rec.gene_id};Name={rec.gene_id}"
            if rec.gene_alias:
                attrs += f";Alias={rec.gene_alias}"
            handle.write(
                "\t".join(
                    [
                        str(rec.chrom),
                        "dmr_validation_framework",
                        "gene",
                        str(int(rec.start)),
                        str(int(rec.end)),
                        ".",
                        str(rec.strand),
                        ".",
                        attrs,
                    ]
                )
                + "\n"
            )

    bed = genes[["chrom", "start", "end", "gene_id", "strand"]].copy()
    bed["bed_start"] = (bed["start"].astype(int) - 1).clip(lower=0)
    bed_out = bed[["chrom", "bed_start", "end", "gene_id"]].copy()
    bed_out["score"] = 0
    bed_out["strand"] = bed["strand"]
    bed_path = out_dir / "metagene_gene_regions.bed"
    bed_out.to_csv(bed_path, sep="\t", index=False, header=False)

    qc = genes[["chrom", "start", "end", "gene_id", "strand", "length_bp", "gene_index"]]
    qc_path = out_dir / "metagene_gene_regions_qc.tsv"
    qc.to_csv(qc_path, sep="\t", index=False)

    seq_map = pd.DataFrame(
        [{"raw_seqname": chrom, "normalized_seqname": chrom, "n_rows": count} for chrom, count in sorted(seq_counts.items())]
    )
    seq_map_path = out_dir / "seqname_mapping.tsv"
    seq_map.to_csv(seq_map_path, sep="\t", index=False)

    feature_summary = pd.DataFrame(
        [{"feature_type": feature, "n_rows": count} for feature, count in sorted(feature_counts.items())]
    )
    feature_summary_path = out_dir / "gff3_feature_summary.tsv"
    feature_summary.to_csv(feature_summary_path, sep="\t", index=False)

    seq_summary_path = out_dir / "gff3_seqname_summary.tsv"
    seq_map.to_csv(seq_summary_path, sep="\t", index=False)

    return genes, {
        "normalized_gene_gff3": normalized_gff,
        "metagene_gene_regions_bed": bed_path,
        "metagene_gene_regions_qc": qc_path,
        "seqname_mapping": seq_map_path,
        "gff3_feature_summary": feature_summary_path,
        "gff3_seqname_summary": seq_summary_path,
    }


def normalize_te_annotation(
    source: Path,
    *,
    out_dir: Path,
    feature_types: set[str],
) -> tuple[pd.DataFrame, dict[str, Path]]:
    rows: list[dict[str, object]] = []
    for row in _iter_gff_rows(source):
        feature = str(row["feature_type"]).lower()
        if feature_types and feature not in feature_types:
            continue
        attrs = row["attrs"]
        te_id = (
            _clean_feature_id(attrs.get("ID"), prefixes=("repeat:", "te:", "transposable_element:"))
            or _clean_feature_id(attrs.get("Name"))
            or f"te_{len(rows) + 1}"
        )
        te_family = attrs.get("Name") or attrs.get("Target") or attrs.get("family") or attrs.get("Family") or feature
        te_class = attrs.get("Class") or attrs.get("class") or attrs.get("Classification") or attrs.get("classification") or feature
        rows.append(
            {
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "te_id": te_id,
                "strand": row["strand"],
                "feature_type": feature,
                "te_family": te_family,
                "te_class": te_class,
                "length_bp": max(0, int(row["end"]) - int(row["start"]) + 1),
            }
        )
    tes = pd.DataFrame(rows)
    if tes.empty:
        return tes, {}
    tes = tes.sort_values(["chrom", "start", "end", "te_id"], kind="mergesort").reset_index(drop=True)
    bed = tes[["chrom", "start", "end", "te_id", "strand", "te_family", "te_class"]].copy()
    bed["bed_start"] = (bed["start"].astype(int) - 1).clip(lower=0)
    bed_out = bed[["chrom", "bed_start", "end", "te_id"]].copy()
    bed_out["score"] = 0
    bed_out["strand"] = bed["strand"]
    bed_out["te_family"] = bed["te_family"]
    bed_out["te_class"] = bed["te_class"]
    bed_path = out_dir / "te_regions.bed"
    bed_out.to_csv(bed_path, sep="\t", index=False, header=False)

    qc = tes[["chrom", "start", "end", "te_id", "strand", "te_family", "te_class", "length_bp"]]
    qc_path = out_dir / "te_regions_qc.tsv"
    qc.to_csv(qc_path, sep="\t", index=False)

    class_summary = (
        tes.groupby(["te_class"], dropna=False)
        .size()
        .reset_index(name="n_regions")
        .sort_values("n_regions", ascending=False)
    )
    class_summary_path = out_dir / "te_class_summary.tsv"
    class_summary.to_csv(class_summary_path, sep="\t", index=False)
    return tes, {
        "te_regions_bed": bed_path,
        "te_regions_qc": qc_path,
        "te_class_summary": class_summary_path,
    }


def _protein_to_gene_id(protein_id: str) -> str:
    value = str(protein_id).strip()
    match = re.match(r"^(Solyc\d+g\d+\.\d+)\.\d+$", value)
    if match:
        return match.group(1)
    parts = value.rsplit(".", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return value


def normalize_go_table(source: Path, *, out_dir: Path) -> tuple[pd.DataFrame, dict[str, Path]]:
    raw = pd.read_csv(source, sep="\t", header=None, comment="#", names=["protein_id", "go_terms"], dtype=str)
    rows: list[dict[str, str]] = []
    for rec in raw.dropna(subset=["protein_id", "go_terms"]).itertuples(index=False):
        protein_id = str(rec.protein_id).strip()
        gene_id = _protein_to_gene_id(protein_id)
        for go_id in re.findall(r"GO:\d+", str(rec.go_terms)):
            rows.append(
                {
                    "gene_id": gene_id,
                    "protein_id": protein_id,
                    "go_id": go_id,
                    "source": "profile_go_table",
                }
            )
    go = pd.DataFrame(rows).drop_duplicates() if rows else pd.DataFrame(columns=["gene_id", "protein_id", "go_id", "source"])
    go = go.sort_values(["gene_id", "protein_id", "go_id"], kind="mergesort").reset_index(drop=True)
    go_path = out_dir / "go_annotations.tsv"
    go.to_csv(go_path, sep="\t", index=False)
    summary = (
        go.groupby("gene_id", dropna=False)
        .size()
        .reset_index(name="n_go_terms")
        .sort_values(["n_go_terms", "gene_id"], ascending=[False, True])
    )
    summary_path = out_dir / "go_gene_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    return go, {"go_annotations": go_path, "go_gene_summary": summary_path}


def write_report(
    *,
    out_dir: Path,
    profile: AnnotationProfile,
    gene_count: int,
    te_count: int | None,
    go_count: int | None,
    outputs: dict[str, Path],
    downloads: list[dict[str, object]],
) -> Path:
    report = out_dir / "annotation_preflight_report.md"
    lines = [
        "# Plant Annotation Preflight Report",
        "",
        f"- profile: `{profile.name}`",
        f"- species: `{profile.species}`",
        f"- assembly: `{profile.assembly}`",
        f"- annotation_version: `{profile.annotation_version}`",
        f"- genes: `{gene_count}`",
        f"- TE/repeat regions: `{te_count if te_count is not None else 'not_requested'}`",
        f"- GO rows: `{go_count if go_count is not None else 'not_requested'}`",
        "",
        "## Outputs",
        "",
    ]
    for name, path in sorted(outputs.items()):
        lines.append(f"- `{name}`: `{path}`")
    lines.extend(["", "## Download Sources", ""])
    for row in downloads:
        lines.append(f"- `{row['status']}` `{row['source']}` -> `{row['path']}` ({row['bytes']} bytes)")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def run(args: argparse.Namespace) -> int:
    base_profile = PROFILES[args.profile]
    profile = replace(
        base_profile,
        species=args.species or base_profile.species,
        assembly=args.assembly or base_profile.assembly,
        annotation_version=args.annotation_version or base_profile.annotation_version,
        output_prefix=args.output_prefix or base_profile.output_prefix,
    )
    out_dir = _default_out_dir(args, profile)
    out_dir.mkdir(parents=True, exist_ok=True)
    download_dir = Path(args.download_dir) if args.download_dir else out_dir / "downloads"
    base_url = args.base_url or profile.base_url

    downloads: list[dict[str, object]] = []
    gene_path, gene_download = _resolve_input(
        local_path=args.gene_gff,
        filename=profile.gene_gff,
        base_url=base_url,
        download_dir=download_dir,
        skip_download=args.skip_download,
        force_download=args.force_download,
    )
    if gene_path is None:
        raise ValueError("gene annotation is required")
    if gene_download:
        downloads.append(gene_download)

    te_path: Path | None = None
    if args.include_te:
        te_path, te_download = _resolve_input(
            local_path=args.te_gff,
            filename=profile.te_gff,
            base_url=base_url,
            download_dir=download_dir,
            skip_download=args.skip_download,
            force_download=args.force_download,
        )
        if te_download:
            downloads.append(te_download)

    go_path: Path | None = None
    if args.include_go:
        go_path, go_download = _resolve_input(
            local_path=args.go_table,
            filename=profile.go_table,
            base_url=base_url,
            download_dir=download_dir,
            skip_download=args.skip_download,
            force_download=args.force_download,
        )
        if go_download:
            downloads.append(go_download)

    outputs: dict[str, Path] = {}
    genes, gene_outputs = normalize_gene_annotation(
        gene_path,
        out_dir=out_dir,
        profile=profile,
        feature_types=_split_csv(args.gene_feature_types),
    )
    outputs.update(gene_outputs)

    te_count: int | None = None
    if te_path is not None:
        tes, te_outputs = normalize_te_annotation(
            te_path,
            out_dir=out_dir,
            feature_types=_split_csv(args.te_feature_types),
        )
        te_count = int(len(tes))
        outputs.update(te_outputs)

    go_count: int | None = None
    if go_path is not None:
        go, go_outputs = normalize_go_table(go_path, out_dir=out_dir)
        go_count = int(len(go))
        outputs.update(go_outputs)

    download_manifest = out_dir / "download_manifest.tsv"
    pd.DataFrame(downloads).to_csv(download_manifest, sep="\t", index=False)
    outputs["download_manifest"] = download_manifest

    report = write_report(
        out_dir=out_dir,
        profile=profile,
        gene_count=int(len(genes)),
        te_count=te_count,
        go_count=go_count,
        outputs=outputs,
        downloads=downloads,
    )
    outputs["annotation_preflight_report"] = report

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "profile": profile.name,
        "species": profile.species,
        "assembly": profile.assembly,
        "annotation_version": profile.annotation_version,
        "source_note": profile.source_note,
        "inputs": {
            "gene_gff": str(gene_path),
            "te_gff": str(te_path) if te_path else None,
            "go_table": str(go_path) if go_path else None,
        },
        "counts": {
            "genes": int(len(genes)),
            "te_regions": te_count,
            "go_rows": go_count,
        },
        "outputs": {key: str(path) for key, path in sorted(outputs.items())},
    }
    manifest_path = out_dir / "annotation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"annotation profile: {profile.name}")
    print(f"genes: {len(genes)}")
    if te_count is not None:
        print(f"TE/repeat regions: {te_count}")
    if go_count is not None:
        print(f"GO rows: {go_count}")
    print(f"annotation manifest: {manifest_path}")
    print(f"preflight report: {report}")
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
