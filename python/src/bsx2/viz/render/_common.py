from __future__ import annotations

import holoviews as hv

from bsx2.validation import NanPolicy


def _hv_init() -> None:
    hv.extension("plotly")


__all__ = ["NanPolicy", "_hv_init"]
