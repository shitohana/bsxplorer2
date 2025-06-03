import numpy as np


class DimMismatchErr(ValueError):
    def __init__(self, expected: int, **kwargs):
        assert kwargs.__len__() > 0
        name, value = list(kwargs.items())[0]
        assert isinstance(value, np.ndarray)
        super().__init__(f"{name} should be {expected}D (got {np.ndim(value)})")


class LengthMismatchErr(ValueError):
    def __init__(self, **kwargs):
        err_string = ", ".join(f"{name} ({len(value)})" for name, value in kwargs.items())
        super().__init__(f"Lengths {err_string} do not match")


class ValueRangeErr(ValueError):
    def __init__(self, name: str, bounds: tuple[float, float] = (0, 1)):
        super().__init__(f"Values of {name} are out of bounds [{bounds[0]}; {bounds[1]}]")


class UnsortedErr(ValueError):
    def __init__(self, name: str, ascending: bool = True):
        if ascending:
            super().__init__(f"{name} should be sorted in ascending order")
        else:
            super().__init__(f"{name} should be sorted in descending order")