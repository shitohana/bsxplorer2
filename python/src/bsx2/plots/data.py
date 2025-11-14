from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple
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
    if not np.issubdtype(a.dtype, np.number):
        return False
    if not np.all(np.isfinite(a)):
        return False
    if a.min() < 0.0 or a.max() > 1.0:
        return False
    if a.size > 1 and not np.all(a[:-1] <= a[1:]):
        return False
    return True


def _is_1d_unit_density(a: np.ndarray) -> bool:
    if not isinstance(a, np.ndarray) or a.ndim != 1:
        return False
    if a.size == 0:
        return True
    if not np.issubdtype(a.dtype, np.number):
        return False
    # Allow NaNs in densities but constrain finite values to [0, 1]
    mask = ~np.isnan(a)
    if mask.any():
        v = a[mask]
        if v.min() < 0.0 or v.max() > 1.0:
            return False
    return True


Pos1D = Annotated[np.ndarray, Is[_is_1d_sorted_unit_positions]]
Density1D = Annotated[np.ndarray, Is[_is_1d_unit_density]]

try:
    import holoviews as hv  # type: ignore
except ModuleNotFoundError:
    hv = None  # lazy-load in to_curve; raise on use if missing


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

    @beartype
    def stack_matrix(self) -> Tuple[np.ndarray, List[str]]:
        if len(self) == 0:
            return np.empty((0, 0), dtype=np.float64), []
        n_bins = len(self.densities[0])
        if not all(len(d) == n_bins for d in self.densities):
            raise ValueError("all regions must have the same number of bins (n_bins)")
        mat = np.vstack([d.astype(np.float64, copy=False) for d in self.densities])
        row_labels = [lbl if lbl is not None else f"region_{i+1}" for i, lbl in enumerate(self.labels)]
        return mat, row_labels

    @beartype
    def to_line_plot(self, agg_fn: Callable = np.nanmean) -> "LinePlotData":
        if len(self) == 0:
            return LinePlotData(np.empty(0), np.empty(0))
        x = self.positions[0].astype(np.float64, copy=False)
        mat, _ = self.stack_matrix()
        y_raw = agg_fn(mat, axis=0)
        if not isinstance(y_raw, np.ndarray):
            y_raw = np.asarray(y_raw)
        if y_raw.ndim != 1 or len(y_raw) != len(x):
            raise ValueError("aggregated values must be a 1D array with the same length as positions")
        y = y_raw.astype(np.float64, copy=False)
        return LinePlotData(x=x, y=y)


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

    @beartype
    def to_curve(self, x_shift: float | int = 0.0, y_shift: float | int = 0.0):
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
