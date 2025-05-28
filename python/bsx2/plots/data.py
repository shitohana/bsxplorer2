from collections.abc import Callable
from typing import Sequence, List, Optional, Tuple

import numpy as np
from pydantic import BaseModel, Field


class DiscreteRegionData(BaseModel):
    positions: List[np.array] = Field(default_factory=list)
    densities: List[np.array] = Field(default_factory=list)
    labels: List[Optional[str]] = Field(default_factory=list)

    def insert(self, positions: Sequence[float], densities: Sequence[float], label: Optional[str] = None) -> None:
        pos_array = np.array(positions, dtype=np.float64)
        density_array = np.array(densities, dtype=np.float64)

        self.positions.append(pos_array)
        self.densities.append(density_array)

        self.labels.append(label)

    def find(self, label: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        if index := self.labels.index(label) is not None:
            return self.positions[index], self.densities[index]
        return None

    def sort(self, fn: Callable[[np.ndarray], np.float64]) -> None:
        indices = list(sorted(range(len(self.positions)), key=lambda i: fn(self.densities[i])))

        self.positions = [self.positions[i] for i in indices]
        self.densities = [self.densities[i] for i in indices]
        self.labels = [self.labels[i] for i in indices]

    def to_line_plot(self) -> "LinePlotData":
        raise NotImplementedError()

    def __len__(self) -> int:
        return len(self.positions)


class PlotData(BaseModel):
    name: Optional[str] = Field(default=None)
    x_ticks: List[float] = Field(default_factory=list)
    x_labels: List[str] = Field(default_factory=list)
    y_ticks: List[float] = Field(default_factory=list)
    y_labels: List[str] = Field(default_factory=list)


class LinePlotData(BaseModel):
    x: np.array = Field(default_factory=lambda: np.empty(dtype=np.float64))
    y: np.array = Field(default_factory=lambda: np.empty(dtype=np.float64))
    x_ticks: List[str] = Field(default_factory=list)
    x_labels: List[str] = Field(default_factory=list)
    name: Optional[str] = Field(default=None)

    def __len__(self) -> int:
        return len(self.x)
