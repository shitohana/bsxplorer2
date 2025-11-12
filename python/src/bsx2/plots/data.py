from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple
import numpy as np

try:
    import holoviews as hv  # type: ignore
except Exception:  # noqa: BLE001
    hv = None  # lazy-load in to_curve


@dataclass
class DiscreteRegionData:
    positions: List[np.ndarray] = field(default_factory=list)   # each: (n_bins,)
    densities: List[np.ndarray] = field(default_factory=list)   # each: (n_bins,)
    labels: List[Optional[str]] = field(default_factory=list)

    def insert(self, positions: np.ndarray, densities: np.ndarray, label: Optional[str] = None) -> None:
        assert positions.ndim == 1 and densities.ndim == 1, "positions/densities must be 1D"
        assert len(positions) == len(densities), "length mismatch positions/densities"
        if len(positions) > 1:
            assert np.all(positions[:-1] <= positions[1:]), "positions must be sorted"
        self.positions.append(positions.astype(np.float64, copy=False))
        self.densities.append(densities.astype(np.float64, copy=False))
        self.labels.append(label)

    def __len__(self) -> int:
        return len(self.positions)

    def stack_matrix(self) -> Tuple[np.ndarray, List[str]]:
        if len(self) == 0:
            return np.empty((0, 0), dtype=np.float64), []
        n_bins = len(self.densities[0])
        assert all(len(d) == n_bins for d in self.densities), "all regions must have same n_bins"
        mat = np.vstack([d.astype(np.float64, copy=False) for d in self.densities])
        row_labels = [lbl if lbl is not None else f"region_{i+1}" for i, lbl in enumerate(self.labels)]
        return mat, row_labels

    def to_line_plot(self, agg_fn: Callable = np.nanmean) -> "LinePlotData":
        if len(self) == 0:
            return LinePlotData(np.empty(0), np.empty(0))
        x = self.positions[0].astype(np.float64, copy=False)
        mat, _ = self.stack_matrix()
        y = agg_fn(mat, axis=0).astype(np.float64, copy=False)
        return LinePlotData(x=x, y=y)


@dataclass
class LinePlotData:
    x: np.ndarray
    y: np.ndarray
    x_ticks: List[float] = field(default_factory=list)
    x_labels: List[str] = field(default_factory=list)
    y_ticks: List[float] = field(default_factory=list)
    y_labels: List[str] = field(default_factory=list)

    def to_curve(self, x_shift: float = 0.0, y_shift: float = 0.0):
        if hv is None:
            import importlib
            globals()["hv"] = importlib.import_module("holoviews")
        curve = hv.Curve((self.x + x_shift, self.y + y_shift))
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

