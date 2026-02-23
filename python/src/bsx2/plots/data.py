from __future__ import annotations
from dataclasses import dataclass, field
from typing import  Annotated, List, Optional

import numpy as np
from beartype import beartype
from beartype.vale import Is

from bsx2.guards import require_equal_length
from bsx2.validation import validate_matrix_shape

# ---------- shared low-level predicates ----------

def _is_ndarray_1d(a: np.ndarray) -> bool:
    return isinstance(a, np.ndarray) and a.ndim == 1


def _is_real_numeric_dtype(a: np.ndarray) -> bool:
    """
    Accept only real numeric dtypes (int/float), reject complex/object/etc.
    """
    return (
        isinstance(a, np.ndarray)
        and (
            np.issubdtype(a.dtype, np.integer)
            or np.issubdtype(a.dtype, np.floating)
        )
    )


def _is_real_numeric_1d(a: np.ndarray) -> bool:
    return _is_ndarray_1d(a) and _is_real_numeric_dtype(a)


# ---------- shape predicate factory ----------

def _make_is_shape_1d(name: str) -> Callable[[np.ndarray], bool]:
    def _pred(a: np.ndarray) -> bool:
        try:
            validate_matrix_shape(a, 1, name=name)
            return True
        except ValueError:
            return False

    _pred.__name__ = f"_is_{name}_shape_1d"
    return _pred


# ---------- domain predicates ----------

def _is_unit_finite(a: np.ndarray) -> bool:
    """
    For positions:
    - empty array -> True
    - real numeric dtype
    - all values are finite
    - all values in [0, 1]
    """
    if not _is_ndarray_1d(a):
        return False

    if a.size == 0:
        return True

    if not _is_real_numeric_dtype(a):
        return False

    return bool(np.isfinite(a).all() and ((0.0 <= a) & (a <= 1.0)).all())


def _is_sorted_nondecreasing(a: np.ndarray) -> bool:
    """
    Non-decreasing order: a[i] <= a[i+1]
    """
    if not _is_ndarray_1d(a):
        return False

    if a.size <= 1:
        return True

    if not _is_real_numeric_dtype(a):
        return False

    return bool((a[:-1] <= a[1:]).all())


def _has_no_inf(a: np.ndarray) -> bool:
    """
    For densities:
    - empty array -> True
    - real numeric dtype
    - +inf / -inf are prohibited
    - NaN is allowed
    """
    if not _is_ndarray_1d(a):
        return False

    if a.size == 0:
        return True

    if not _is_real_numeric_dtype(a):
        return False

    return bool(not np.isinf(a).any())


def _finite_values_in_unit(a: np.ndarray) -> bool:
    """
    For densities:
    - check only finite values
    - if there are no finite values (e.g. all NaN), return True
    """
    if not _is_ndarray_1d(a):
        return False

    if a.size == 0:
        return True

    if not _is_real_numeric_dtype(a):
        return False

    finite_mask = np.isfinite(a)
    if not finite_mask.any():
        return True

    finite_vals = a[finite_mask]
    return bool(((finite_vals >= 0.0) & (finite_vals <= 1.0)).all())


# ---------- shape predicates for specific fields ----------

_is_positions_shape_1d = _make_is_shape_1d("positions")
_is_densities_shape_1d = _make_is_shape_1d("densities")


# ---------- typed aliases (contracts) ----------

Pos1D = Annotated[
    np.ndarray,
    Is[_is_positions_shape_1d] & Is[_is_unit_finite] & Is[_is_sorted_nondecreasing],
]

Density1D = Annotated[
    np.ndarray,
    Is[_is_densities_shape_1d] & Is[_has_no_inf] & Is[_finite_values_in_unit],
]


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
        require_equal_length(
            positions,
            densities,
            left_name="positions",
            right_name="densities",
            message="length mismatch between positions and densities",
        )

        pos = np.asarray(positions, dtype=np.float64)
        den = np.asarray(densities, dtype=np.float64)

        pos.setflags(write=False)
        den.setflags(write=False)

        self.positions.append(pos)
        self.densities.append(den)
        self.labels.append(label)