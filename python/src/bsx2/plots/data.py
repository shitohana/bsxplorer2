from __future__ import annotations
from dataclasses import dataclass, field
from typing import Annotated, Callable, List, Optional, TypeAlias

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
    if a.size == 0:
        return True

    if not _is_real_numeric_dtype(a):
        return False

    return bool(np.isfinite(a).all() and ((0.0 <= a) & (a <= 1.0)).all())


def _is_sorted_nondecreasing(a: np.ndarray) -> bool:
    """
    Non-decreasing order: a[i] <= a[i+1]
    """
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

# Internal container storage type: arrays are already validated at insert()
# boundaries and then frozen via setflags(write=False).
_StoredArray: TypeAlias = np.ndarray


@beartype
@dataclass
class DiscreteRegionData:
    positions: List[_StoredArray] = field(default_factory=list)   # each: (n_bins,)
    densities: List[_StoredArray] = field(default_factory=list)   # each: (n_bins,)
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

        positions.setflags(write=False)
        densities.setflags(write=False)

        pos_stored: _StoredArray = positions
        den_stored: _StoredArray = densities

        self.positions.append(pos_stored)
        self.densities.append(den_stored)
        self.labels.append(label)

    @beartype
    def clone(self, *, deep: bool = True) -> "DiscreteRegionData":
        """
        Clone container.

        deep=True  -> copies arrays
        deep=False -> reuses array objects (safe if arrays remain read-only)
        """
        out = DiscreteRegionData()

        if deep:
            out.positions = [np.array(p, dtype=np.float64, copy=True) for p in self.positions]
            out.densities = [np.array(d, dtype=np.float64, copy=True) for d in self.densities]
            for a in out.positions:
                a.setflags(write=False)
            for a in out.densities:
                a.setflags(write=False)
        else:
            out.positions = list(self.positions)
            out.densities = list(self.densities)

        out.labels = list(self.labels)
        return out

    @beartype
    def mean_density_per_region(self) -> np.ndarray:
        """
        Mean density for each region (ignores NaN).
        Returns NaN for empty regions or regions with no finite values.
        """
        out = np.full(len(self.densities), np.nan, dtype=np.float64)

        for i, d in enumerate(self.densities):
            if d.size == 0:
                continue
            finite = np.isfinite(d)
            if not np.any(finite):
                continue
            out[i] = float(np.mean(d[finite]))

        return out

    @beartype
    def nan_counts(self) -> np.ndarray:
        """
        Number of NaN values in densities for each region.
        """
        out = np.zeros(len(self.densities), dtype=np.int64)
        for i, d in enumerate(self.densities):
            if d.size == 0:
                out[i] = 0
            else:
                out[i] = int(np.isnan(d).sum())
        return out

    @beartype
    def fill_nan(
        self,
        value: float = 0.0,
        *,
        in_place: bool = False,
    ) -> "DiscreteRegionData":
        """
        Replace NaN in densities with a constant value.
        Positions and labels are preserved.
        """
        target = self if in_place else self.clone(deep=False)

        new_densities: list[_StoredArray] = []
        for d in target.densities:
            if d.size == 0 or not np.isnan(d).any():
                # keep original reference (already read-only)
                new_densities.append(d)
                continue
            d2 = np.asarray(np.where(np.isnan(d), value, d), dtype=np.float64)
            d2.setflags(write=False)
            new_densities.append(d2)

        target.densities = new_densities
        return target

    @beartype
    def filter(
        self,
        *,
        labels: Optional[List[str]] = None,
        predicate: Optional[Callable[[np.ndarray, np.ndarray, Optional[str], int], bool]] = None,
        in_place: bool = False,
    ) -> "DiscreteRegionData":
        """
        Filter regions by labels and/or predicate.

        predicate signature:
            (positions, densities, label, index) -> bool
        """
        label_set = set(labels) if labels is not None else None

        keep_idx: list[int] = []
        for i, (p, d, lbl) in enumerate(zip(self.positions, self.densities, self.labels)):
            if label_set is not None and lbl not in label_set:
                continue
            if predicate is not None and not predicate(p, d, lbl, i):
                continue
            keep_idx.append(i)

        target = self if in_place else DiscreteRegionData()

        new_positions = [self.positions[i] for i in keep_idx]
        new_densities = [self.densities[i] for i in keep_idx]
        new_labels = [self.labels[i] for i in keep_idx]

        if in_place:
            self.positions = new_positions
            self.densities = new_densities
            self.labels = new_labels
            return self

        target.positions = new_positions
        target.densities = new_densities
        target.labels = new_labels
        return target

    @beartype
    def sort(
        self,
        *,
        by: str = "label",  # "label" | "mean_density" | "length"
        reverse: bool = False,
        in_place: bool = False,
    ) -> "DiscreteRegionData":
        """
        Sort regions (stable sort) by label / mean_density / length.
        """
        if by not in {"label", "mean_density", "length"}:
            raise ValueError("sort.by must be one of: 'label', 'mean_density', 'length'")

        means = None
        if by == "mean_density":
            means = self.mean_density_per_region()

        def _key(i: int):
            if by == "label":
                lbl = self.labels[i]
                # None goes last by default in ascending
                return (lbl is None, "" if lbl is None else str(lbl))
            if by == "length":
                return int(self.positions[i].size)
            # by == "mean_density"
            v = float(means[i])  # type: ignore[index]
            # NaN goes last in ascending
            return (np.isnan(v), v if np.isfinite(v) else 0.0)

        order = sorted(range(len(self.positions)), key=_key, reverse=reverse)

        target = self if in_place else DiscreteRegionData()
        new_positions = [self.positions[i] for i in order]
        new_densities = [self.densities[i] for i in order]
        new_labels = [self.labels[i] for i in order]

        if in_place:
            self.positions = new_positions
            self.densities = new_densities
            self.labels = new_labels
            return self

        target.positions = new_positions
        target.densities = new_densities
        target.labels = new_labels
        return target
