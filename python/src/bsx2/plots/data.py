from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Callable

try:  # Prefer beartype-aware typing to silence PEP585 warnings
    from beartype.typing import List, Optional, Tuple  # type: ignore
except Exception:  # pragma: no cover - fallback
    from typing import List, Optional, Tuple  # type: ignore
import numpy as np

# beartype: runtime type checking with value constraints
try:
    from beartype import beartype  # type: ignore
    from beartype.vale import Is  # type: ignore
    try:  # Python >=3.9 has typing.Annotated; fall back to typing_extensions
        from typing import Annotated  # type: ignore
    except Exception:  # pragma: no cover - platform variance
        from typing_extensions import Annotated  # type: ignore
except ModuleNotFoundError:  # graceful fallback if beartype is not installed
    try:
        from typing import Annotated  # type: ignore
    except Exception:  # pragma: no cover
        from typing_extensions import Annotated  # type: ignore

    def beartype(obj):  # type: ignore
        return obj

    class Is:  # type: ignore
        def __class_getitem__(cls, item):
            # Return the predicate itself as metadata for Annotated; no-op without beartype
            return item


# Validators for arrays in [0, 1], 1D, and ordering where required
def _is_1d_sorted_unit_positions(a: np.ndarray) -> bool:
    if not isinstance(a, np.ndarray) or a.ndim != 1:
        return False
    if a.size == 0:
        return True
    if not np.issubdtype(a.dtype, np.number) or not np.all(np.isfinite(a)):
        return False
    return (
        0.0 <= a.min() <= a.max() <= 1.0
        and (a.size == 1 or np.all(a[:-1] <= a[1:]))  # sorted non-decreasing
    )


def _is_1d_unit_density(a: np.ndarray) -> bool:
    if not isinstance(a, np.ndarray) or a.ndim != 1:
        return False
    if a.size == 0:
        return True
    if not np.issubdtype(a.dtype, np.number):
        return False
    if np.isinf(a).any():
        return False
    finite = np.isfinite(a)
    if not finite.any():
        return True
    return (a[finite] >= 0.0).all() and (a[finite] <= 1.0).all()


Pos1D = Annotated[np.ndarray, Is[_is_1d_sorted_unit_positions]]
Density1D = Annotated[np.ndarray, Is[_is_1d_unit_density]]

try:
    import holoviews as hv  # type: ignore
except ModuleNotFoundError:
    hv = None  # lazy-load in to_curve; raise on use if missing

@dataclass
class MetageneData:
    """Единый контейнер для метагенной матрицы и сегментов."""
    matrix: np.ndarray
    region_ids: List[str]
    bins: np.ndarray
    segments: List[Any]  # Segment, но избегаем циклического импорта

    @classmethod
    def from_discrete(cls, drd: "DiscreteRegionData", segments: List[Any], *, drop_nan_rows: bool = False) -> "MetageneData":
        """Строит метагенную матрицу из DRD; можно выкинуть строки с NaN."""
        mat, labels, grid = drd._stack_to_common_grid()
        if drop_nan_rows and mat.size > 0:
            finite_mask = np.isfinite(mat).all(axis=1)
            mat = mat[finite_mask]
            labels = [l for l, k in zip(labels, finite_mask) if k]
        return cls(matrix=mat, region_ids=labels, bins=grid, segments=segments)

    def line_profile(self, agg: str = "mean", *, as_percent: bool = False, drop_nan_rows: bool = False) -> tuple[np.ndarray, np.ndarray]:
        data = self.matrix * (100.0 if as_percent else 1.0)
        if drop_nan_rows and data.size > 0:
            data = data[np.isfinite(data).any(axis=1)]
        agg_map = {
            "mean": np.nanmean,
            "median": np.nanmedian,
            "max": np.nanmax,
            "min": np.nanmin,
        }
        fn = agg_map.get(agg)
        if fn is None:
            raise ValueError(f"Unsupported agg: {agg}")
        y = fn(data, axis=0)
        # если весь столбец NaN — заменяем на глобальное среднее
        if np.isnan(y).any():
            global_mean = np.nanmean(data) if np.isfinite(data).any() else 0.0
            y = np.where(np.isnan(y), global_mean, y)
        return self.bins, y

    def dist_data(self, *, as_percent: bool = False, per_region: bool = False, drop_nan_rows: bool = False):
        data = self.matrix * (100.0 if as_percent else 1.0)
        if drop_nan_rows and data.size > 0:
            data = data[np.isfinite(data).any(axis=1)]
        if per_region:
            out = []
            for rid, row in zip(self.region_ids, data):
                if np.isfinite(row).any():
                    out.append((rid, float(np.nanmean(row))))
            return out
        out = []
        for b_idx, col in enumerate(data.T):
            for val in col:
                if np.isfinite(val):
                    out.append((b_idx, float(val)))
        return out

    def heatmap_data(self, *, as_percent: bool = False, drop_nan_rows: bool = False):
        data = self.matrix * (100.0 if as_percent else 1.0)
        regions = self.region_ids
        if drop_nan_rows and data.size > 0:
            mask = np.isfinite(data).any(axis=1)
            data = data[mask]
            regions = [r for r, k in zip(self.region_ids, mask) if k]
        return data, self.region_ids, list(range(data.shape[1]))


@beartype
@dataclass
class DiscreteRegionData:
    positions: List[np.ndarray] = field(default_factory=list)   # each: (n_bins,)
    densities: List[np.ndarray] = field(default_factory=list)   # each: (n_bins,)
    labels: List[Optional[str]] = field(default_factory=list)

    @beartype
    def insert(self, positions: Pos1D, densities: Density1D, label: Optional[str] = None) -> None:
        if len(positions) != len(densities):
            raise ValueError("length mismatch between positions and densities")
        self.positions.append(positions.astype(np.float64, copy=False))
        self.densities.append(densities.astype(np.float64, copy=False))
        self.labels.append(label)

    def __len__(self) -> int:
        return len(self.positions)

    def _stack_to_common_grid(self) -> Tuple[np.ndarray, List[str], np.ndarray]:
        if len(self) == 0:
            return np.empty((0, 0), dtype=np.float64), [], np.empty(0, dtype=np.float64)

        # Choose a common grid by the maximum number of bins across regions (normalized 0..1)
        target_bins = max(len(d) for d in self.densities)
        if target_bins == 0:
            return np.empty((len(self), 0), dtype=np.float64), [], np.empty(0, dtype=np.float64)
        # Bin edges and centers
        edges = np.linspace(0.0, 1.0, target_bins + 1, dtype=np.float64)
        centers = (edges[:-1] + edges[1:]) / 2.0

        rows = []
        for pos, dens in zip(self.positions, self.densities):
            # Simple binning: points whose coordinate falls into bin edges
            binned = np.full(target_bins, np.nan, dtype=np.float64)
            if len(pos) > 0:
                idx = np.clip(np.digitize(pos, edges[1:-1], right=True), 0, target_bins - 1)
                for b in range(target_bins):
                    mask = idx == b
                    if np.any(mask):
                        binned[b] = float(np.mean(dens[mask]))
            rows.append(binned)

        mat = np.vstack(rows)
        row_labels = [lbl if lbl is not None else f"region_{i+1}" for i, lbl in enumerate(self.labels)]
        return mat, row_labels, centers

    @beartype
    def stack_matrix(self) -> Tuple[np.ndarray, List[str]]:
        mat, row_labels, _ = self._stack_to_common_grid()
        return mat, row_labels


@beartype
@dataclass
class LinePlotData:
    x: Pos1D
    y: Density1D
    x_ticks: List[float] = field(default_factory=list)
    x_labels: List[str] = field(default_factory=list)
    y_ticks: List[float] = field(default_factory=list)
    y_labels: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if len(self.x) != len(self.y):
            raise ValueError("x and y must have the same length")
        if self.x_labels and (len(self.x_ticks) != len(self.x_labels)):
            raise ValueError("x_ticks and x_labels must have the same length when labels are provided")
        if self.y_labels and (len(self.y_ticks) != len(self.y_labels)):
            raise ValueError("y_ticks and y_labels must have the same length when labels are provided")

    @classmethod
    @beartype
    def from_discrete(cls, drd: "DiscreteRegionData", agg_fn: Callable = np.nanmean) -> "LinePlotData":
        if len(drd) == 0:
            return cls(np.empty(0), np.empty(0))
        mat, _, grid = drd._stack_to_common_grid()
        y_raw = agg_fn(mat, axis=0)
        if not isinstance(y_raw, np.ndarray):
            y_raw = np.asarray(y_raw)
        if y_raw.ndim != 1 or len(y_raw) != len(grid):
            raise ValueError("aggregated values must be a 1D array with the same length as positions")
        y = y_raw.astype(np.float64, copy=False)
        return cls(x=grid, y=y)

    @beartype
    def to_curve(self, x_shift: float | int = 0.0, y_shift: float | int = 0.0) -> Any:
        if hv is None:
            import importlib
            try:
                globals()["hv"] = importlib.import_module("holoviews")
            except ModuleNotFoundError as e:
                raise ImportError("holoviews is required to create curves; install with 'pip install holoviews'") from e
        # Re-read possible imported hv from globals
        local_hv = globals().get("hv")
        if local_hv is None:
            raise RuntimeError("Failed to import holoviews")
        curve = local_hv.Curve((self.x + x_shift, self.y + y_shift))
        if self.x_ticks:
            if self.x_labels:
                curve = curve.opts(xticks=list(zip([t + x_shift for t in self.x_ticks], self.x_labels)))
            else:
                curve = curve.opts(xticks=[t + x_shift for t in self.x_ticks])
        if self.y_ticks:
            if self.y_labels:
                curve = curve.opts(yticks=list(zip([t + y_shift for t in self.y_ticks], self.y_labels)))
            else:
                curve = curve.opts(yticks=[t + y_shift for t in self.y_ticks])
        return curve
