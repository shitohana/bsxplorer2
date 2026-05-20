"""On-disk cache for extracted DiscreteRegionData.

Purpose:
    Store precomputed regional point data so display-level binning/rendering can
    be repeated without rereading methylation files.

Limitations:
    The cache preserves extracted points only. Changing biological extraction
    layout, flanks, region definitions, context, or strand policy requires
    recomputation.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

try:
    from .data import DiscreteRegionData
except ModuleNotFoundError as exc:
    if exc.name not in {"beartype", "bsx2._bsx2"}:
        raise

    class DiscreteRegionData:  # type: ignore[no-redef]
        def __init__(self) -> None:
            self.positions: list[np.ndarray] = []
            self.densities: list[np.ndarray] = []
            self.labels: list[str | None] = []

        def insert(self, positions: np.ndarray, densities: np.ndarray, label: str | None = None) -> None:
            if len(positions) != len(densities):
                raise ValueError("length mismatch between positions and densities")
            self.insert_unchecked(positions, densities, label)

        def insert_unchecked(self, positions: np.ndarray, densities: np.ndarray, label: str | None = None) -> None:
            self.positions.append(np.asarray(positions, dtype=float))
            self.densities.append(np.asarray(densities, dtype=float))
            self.labels.append(label)


SCHEMA_VERSION = "1.0"


def _as_object_array(values: list[Any]) -> np.ndarray:
    arr = np.empty(len(values), dtype=object)
    for i, value in enumerate(values):
        arr[i] = np.asarray(value, dtype=float)
    return arr


def discrete_region_data_fingerprint(drd: DiscreteRegionData) -> str:
    h = hashlib.sha256()
    for positions, densities, label in zip(drd.positions, drd.densities, drd.labels):
        h.update(np.asarray(positions, dtype=float).tobytes())
        h.update(np.asarray(densities, dtype=float).tobytes())
        h.update(str(label).encode("utf-8"))
    return h.hexdigest()


def save_discrete_region_data(drd: DiscreteRegionData, out_dir: str | Path, metadata: dict[str, Any] | None = None) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "data.npz", positions=_as_object_array(drd.positions), densities=_as_object_array(drd.densities), labels=np.asarray(drd.labels, dtype=object))
    (out / "metadata.json").write_text(json.dumps(metadata or {}, indent=2, sort_keys=True), encoding="utf-8")
    write_discrete_region_cache_manifest(out, drd, metadata or {})


def load_discrete_region_data(out_dir: str | Path) -> DiscreteRegionData:
    out = Path(out_dir)
    data_path = out / "data.npz"
    if not data_path.exists():
        raise FileNotFoundError(f"Missing cache data file: {data_path}")
    npz = np.load(data_path, allow_pickle=True)
    drd = DiscreteRegionData()
    for positions, densities, label in zip(npz["positions"], npz["densities"], npz["labels"]):
        if hasattr(drd, "insert_unchecked"):
            drd.insert_unchecked(np.asarray(positions, dtype=float), np.asarray(densities, dtype=float), None if label is None else str(label))
        else:
            drd.insert(np.asarray(positions, dtype=float), np.asarray(densities, dtype=float), None if label is None else str(label))
    return drd


def write_discrete_region_cache_manifest(out_dir: str | Path, drd: DiscreteRegionData, metadata: dict[str, Any] | None = None) -> dict[str, Any]:
    out = Path(out_dir)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "hash": discrete_region_data_fingerprint(drd),
        "n_regions": len(drd.positions),
        "metadata": metadata or {},
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
