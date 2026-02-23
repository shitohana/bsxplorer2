from __future__ import annotations
from dataclasses import dataclass, field

from typing import List, Optional

import numpy as np
from bsx2.guards import require_equal_length
from bsx2.validation import validate_matrix_shape

# beartype: runtime type checking with value constraints
from beartype import beartype
from beartype.vale import Is
from typing import Annotated

# Validators for arrays in [0, 1], 1D, and ordering where required
def _is_1d_sorted_unit_positions(a: np.ndarray) -> bool:
    try:
        a = validate_matrix_shape(a, 1, name="positions")
    except ValueError:
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
    try:
        a = validate_matrix_shape(a, 1, name="densities")
    except ValueError:
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


@beartype
@dataclass
class DiscreteRegionData:
    positions: List[np.ndarray] = field(default_factory=list)   # each: (n_bins,)
    densities: List[np.ndarray] = field(default_factory=list)   # each: (n_bins,)
    labels: List[Optional[str]] = field(default_factory=list)

    @beartype
    def insert(
        self,
        positions: Pos1D,
        densities: Density1D,
        label: Optional[str] = None,
    ) -> None:
        positions = validate_matrix_shape(positions, 1, name="positions")
        densities = validate_matrix_shape(densities, 1, name="densities")
        require_equal_length(
            positions,
            densities,
            left_name="positions",
            right_name="densities",
            message="length mismatch between positions and densities",
        )
        self.positions.append(positions.astype(np.float64, copy=False))
        self.densities.append(densities.astype(np.float64, copy=False))
        self.labels.append(label)
