import heapq
from collections.abc import Callable
from typing import Any
from typing import Sequence, List, Optional, Tuple

import holoviews as hv
import numpy as np
import numpy.typing as npt
import pandas as pd
from numpy import floating
from pydantic import BaseModel, Field, model_validator

from bsx2.utils import LengthMismatchErr, DimMismatchErr, ValueRangeErr, UnsortedErr


class DiscreteRegionData(BaseModel):
    model_config = {"arbitrary_types_allowed": True}

    positions: List[npt.NDArray[np.float64]] = Field(default_factory=list)
    densities: List[npt.NDArray[np.float64]] = Field(default_factory=list)
    labels: List[Optional[str]] = Field(default_factory=list)

    @model_validator(mode="after")
    def _validate_self(self) -> "DiscreteRegionData":
        assert len(self.positions) == len(self.densities), LengthMismatchErr(
            positions=self.positions, densities=self.densities
        )
        assert (
            len(
                neq := list(
                    (p, d)
                    for p, d in zip(self.positions, self.densities)
                    if len(p) != len(d)
                )
            )
            == 0
        ), LengthMismatchErr(positions_arr=neq[0][0], densities_arr=neq[0][1])
        if len(self.labels) != 0:
            assert len(self.labels) == len(self.positions), LengthMismatchErr(
                labels=self.labels, positions=self.positions
            )
        assert all(np.all(arr[:-1] <= arr[1:]) for arr in self.positions), UnsortedErr(
            "All arrays"
        )

        return self

    def insert(
        self,
        positions: npt.NDArray[np.float64],
        densities: npt.NDArray[np.float64],
        label: Optional[str] = None,
    ) -> None:
        assert len(positions) == len(densities), LengthMismatchErr(
            positions=positions, densities=densities
        )
        assert np.all(positions[:-1] <= positions[1:]), UnsortedErr("positions")

        self.positions.append(positions)
        self.densities.append(densities)

        self.labels.append(label)

    def find(self, label: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        if index := self.labels.index(label) is not None:
            return self.positions[index], self.densities[index]
        return None

    def to_line_plot(
        self,
        n_points: int = 100,
        agg_fn: Callable[[np.ndarray], np.floating[Any]] = lambda arr: np.mean(arr),
    ) -> "LinePlotData":
        if self.__len__() == 0:
            return LinePlotData()

        combined_points = [
            np.column_stack((p, d)) for p, d in zip(self.positions, self.densities)
        ]
        sorted_points = np.array(
            list(
                heapq.merge(*combined_points, key=lambda row: row[0]),
            )
        )

        positions, densities = sorted_points.T
        x_fixed = np.linspace(0, 1, n_points, endpoint=False)
        split_indices = np.searchsorted(positions, x_fixed, side="right")[1:]

        position_chunks = np.split(positions, split_indices)
        density_chunks = np.split(densities, split_indices)

        x = np.array([agg_fn(chunk) for chunk in position_chunks])
        y = np.array([agg_fn(chunk) for chunk in density_chunks])
        return LinePlotData(x, y)

    def n_points(self) -> int:
        return sum(len(a) for a in self.positions)

    def mean_points(self) -> floating[Any]:
        return np.mean([len(a) for a in self.positions])

    def to_pandas(self) -> pd.DataFrame:
        return pd.DataFrame.from_dict(
            {
                "positions": self.positions,
                "densities": self.densities,
                "labels": self.labels,
            },
            orient="columns",
        )

    def __len__(self) -> int:
        return len(self.positions)


class PlotData(BaseModel):
    model_config = {"arbitrary_types_allowed": True}

    name: Optional[str] = Field(default=None)
    x_ticks: List[float] = Field(default_factory=list)
    x_labels: List[str] = Field(default_factory=list)
    y_ticks: List[float] = Field(default_factory=list)
    y_labels: List[str] = Field(default_factory=list)

    def __init__(
        self,
        name: Optional[str] = None,
        *,
        x_ticks: Optional[Sequence[float]] = None,
        y_ticks: Optional[Sequence[float]] = None,
        x_labels: Optional[Sequence[str]] = None,
        y_labels: Optional[Sequence[str]] = None,
    ):
        kwargs = dict(name=name)

        if x_ticks is not None:
            kwargs |= dict(x_ticks=list(x_ticks))
        if y_ticks is not None:
            kwargs |= dict(y_ticks=list(y_ticks))
        if x_labels is not None:
            kwargs |= dict(x_labels=list(x_labels))
        if y_labels is not None:
            kwargs |= dict(y_labels=list(y_labels))

        super().__init__(**kwargs)

    @property
    def has_x_labels(self):
        return len(self.x_labels) > 0

    @property
    def has_y_labels(self):
        return len(self.y_labels) > 0

    def set_x_ticks(
        self, ticks: Sequence[float], labels: Optional[Sequence[str]] = None
    ):
        if labels is not None:
            assert len(ticks) == len(labels), LengthMismatchErr(
                ticks=ticks, labels=labels
            )

        ticks = list(ticks)
        labels = list(labels) if labels is not None else []

        self.x_ticks, self.x_labels = ticks, labels

    def set_y_ticks(
        self, ticks: Sequence[float], labels: Optional[Sequence[str]] = None
    ):
        if labels is not None:
            assert len(ticks) == len(labels), LengthMismatchErr(
                ticks=ticks, labels=labels
            )

        ticks = list(ticks)
        labels = list(labels) if labels is not None else []

        self.y_ticks, self.y_labels = ticks, labels

    @model_validator(mode="after")
    def _validate_self(self) -> "PlotData":
        # Validation for x and y removed from here
        if len(self.x_labels) != 0:
            assert len(self.x_labels) == len(self.x_ticks), LengthMismatchErr(
                x_labels=self.x_labels, x_ticks=self.x_ticks
            )
        if len(self.y_labels) != 0:
            assert len(self.y_labels) == len(self.y_ticks), LengthMismatchErr(
                y_labels=self.y_labels, y_ticks=self.y_ticks
            )

        return self


class LinePlotData(PlotData):
    x: npt.NDArray[np.float64] = Field(
        default_factory=lambda: np.empty(0, dtype=np.float64)
    )
    y: npt.NDArray[np.float64] = Field(
        default_factory=lambda: np.empty(0, dtype=np.float64)
    )

    def __init__(
        self,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        name: Optional[str] = None,
        *,
        x_ticks: Optional[Sequence[float]] = None,
        y_ticks: Optional[Sequence[float]] = None,
        x_labels: Optional[Sequence[str]] = None,
        y_labels: Optional[Sequence[str]] = None,
    ):
        # Pass parent-relevant args to parent init
        super().__init__(
            name=name,
            x_ticks=x_ticks,
            y_ticks=y_ticks,
            x_labels=x_labels,
            y_labels=y_labels,
        )
        self.x = x
        self.y = y

    @model_validator(mode="after")
    def _validate_values(self) -> "LinePlotData":
        assert self.x.ndim == 1, DimMismatchErr(1, x=self.x)
        assert self.y.ndim == 1, DimMismatchErr(1, y=self.y)
        assert len(self.x) == len(self.y), LengthMismatchErr(x=self.x, y=self.y)

        assert all(np.isnan(v) or 0 <= v <= 1 for v in self.x), ValueRangeErr("x")
        assert all(np.isnan(v) or 0 <= v <= 1 for v in self.y), ValueRangeErr("y")
        assert all(np.isnan(v) or 0 <= v <= 1 for v in self.x_ticks), ValueRangeErr(
            "x_ticks"
        )
        assert all(np.isnan(v) or 0 <= v <= 1 for v in self.y_ticks), ValueRangeErr(
            "y_ticks"
        )

        return self

    def to_curve(
        self,
        x_shift: float = 0,
        y_shift: float = 0,
    ) -> hv.Curve:
        curve = hv.Curve((self.x + x_shift, self.y + y_shift))

        if self.x_ticks:
            if self.has_x_labels:
                curve = curve.opts(
                    xticks=list(zip([t + x_shift for t in self.x_ticks], self.x_labels))
                )
            else:
                curve = curve.opts(xticks=[t + x_shift for t in self.x_ticks])

        if self.y_ticks:
            if self.has_y_labels:
                curve = curve.opts(
                    yticks=list(zip([t + y_shift for t in self.y_ticks], self.y_labels))
                )
            else:
                curve = curve.opts(yticks=[t + y_shift for t in self.y_ticks])
        return curve

    def __len__(self) -> int:
        return len(self.x)
