from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple
import numpy as np

try:
    import holoviews as hv  # type: ignore
except ModuleNotFoundError:
    hv = None  # lazy-load in to_curve; raise on use if missing


@dataclass
class DiscreteRegionData:
    positions: List[np.ndarray] = field(default_factory=list)   # each: (n_bins,)
    densities: List[np.ndarray] = field(default_factory=list)   # each: (n_bins,)
    labels: List[Optional[str]] = field(default_factory=list)

    def insert(self, positions: np.ndarray, densities: np.ndarray, label: Optional[str] = None) -> None:
        if not isinstance(positions, np.ndarray) or not isinstance(densities, np.ndarray):
            raise TypeError("positions and densities must be numpy.ndarray")
        if positions.ndim != 1 or densities.ndim != 1:
            raise ValueError("positions/densities must be 1D arrays")
        if len(positions) != len(densities):
            raise ValueError("length mismatch between positions and densities")
        if len(positions) > 1 and not np.all(positions[:-1] <= positions[1:]):
            raise ValueError("positions must be sorted in non-decreasing order")
        self.positions.append(positions.astype(np.float64, copy=False))
        self.densities.append(densities.astype(np.float64, copy=False))
        self.labels.append(label)

    def __len__(self) -> int:
        return len(self.positions)

    def stack_matrix(self) -> Tuple[np.ndarray, List[str]]:
        if len(self) == 0:
            return np.empty((0, 0), dtype=np.float64), []
        n_bins = len(self.densities[0])
        if not all(len(d) == n_bins for d in self.densities):
            raise ValueError("all regions must have the same number of bins (n_bins)")
        mat = np.vstack([d.astype(np.float64, copy=False) for d in self.densities])
        row_labels = [lbl if lbl is not None else f"region_{i+1}" for i, lbl in enumerate(self.labels)]
        return mat, row_labels

    def to_line_plot(self, agg_fn: Callable = np.nanmean) -> "LinePlotData":
        if len(self) == 0:
            return LinePlotData(np.empty(0), np.empty(0))
        if not callable(agg_fn):
            raise TypeError("agg_fn must be callable")
        x = self.positions[0].astype(np.float64, copy=False)
        mat, _ = self.stack_matrix()
        y_raw = agg_fn(mat, axis=0)
        if not isinstance(y_raw, np.ndarray):
            y_raw = np.asarray(y_raw)
        if y_raw.ndim != 1 or len(y_raw) != len(x):
            raise ValueError("aggregated values must be a 1D array with the same length as positions")
        y = y_raw.astype(np.float64, copy=False)
        return LinePlotData(x=x, y=y)


@dataclass
class LinePlotData:
    x: np.ndarray
    y: np.ndarray
    x_ticks: List[float] = field(default_factory=list)
    x_labels: List[str] = field(default_factory=list)
    y_ticks: List[float] = field(default_factory=list)
    y_labels: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not isinstance(self.x, np.ndarray) or not isinstance(self.y, np.ndarray):
            raise TypeError("x and y must be numpy.ndarray")
        if self.x.ndim != 1 or self.y.ndim != 1:
            raise ValueError("x and y must be 1D arrays")
        if len(self.x) != len(self.y):
            raise ValueError("x and y must have the same length")
        if self.x_labels and (len(self.x_ticks) != len(self.x_labels)):
            raise ValueError("x_ticks and x_labels must have the same length when labels are provided")
        if self.y_labels and (len(self.y_ticks) != len(self.y_labels)):
            raise ValueError("y_ticks and y_labels must have the same length when labels are provided")

    def to_curve(self, x_shift: float = 0.0, y_shift: float = 0.0):
        if not isinstance(x_shift, (int, float)) or not isinstance(y_shift, (int, float)):
            raise TypeError("x_shift and y_shift must be numbers")
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
