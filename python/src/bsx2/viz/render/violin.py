from __future__ import annotations

import holoviews as hv
from beartype.typing import Optional

from bsx2.validation import validate_positive_int

from ..compute.data import DiscreteRegionData
from ..compute.distribution import DistributionMode, build_violin_distribution_data
from ..compute.metagene import MetageneProfileSegment
from ._common import NanPolicy, _hv_init


def violin_plot(
    drd: DiscreteRegionData,
    *,
    segments: list[MetageneProfileSegment] | None = None,
    n_windows: Optional[int] = None,
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    nan_policy: NanPolicy = NanPolicy.DROP,
    per_region: bool = False,
    distribution_mode: DistributionMode = "windows",
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
):
    """
    Build a HoloViews violin plot for normalized-profile methylation distributions.

    Parameters
    ----------
    drd
        Discrete normalized profiles to render.
    segments
        Optional normalized-profile segment layout used for binning.
    n_windows
        Number of output windows. Defaults to the total segment bin count.
    as_percent
        If `True`, scale methylation values to percentages.
    nan_fill
        Optional value used to replace NaNs before aggregation.
    nan_policy
        Policy controlling how NaN values are handled during distribution prep.
    per_region
        If `True`, keep individual region values instead of summarizing by bin.
    distribution_mode
        Grouping mode for metagene distributions. ``"windows"`` builds one
        violin per metagene window; ``"segments"`` summarizes one value per
        region in each named segment.
    title
        Optional plot title.
    width, height
        Optional plot size in pixels.

    Returns
    -------
    object
        HoloViews violin plot object.

    Notes
    -----
    For metagene interpretation, the primary grouped mode is
    ``per_region=False``. In that mode the violin summarizes distributions
    across the selected gene or region set by window or named segment.
    """
    _hv_init()
    data = build_violin_distribution_data(
        drd,
        segments=segments,
        n_windows=n_windows,
        as_percent=as_percent,
        nan_fill=nan_fill,
        nan_policy=nan_policy,
        per_region=per_region,
        distribution_mode=distribution_mode,
    )

    opts_kwargs = dict(
        xlabel=data.x_label,
        ylabel=data.y_label,
        show_legend=False,
        title=title or "Scaled region profile - Violin",
        box=True,
    )
    if as_percent:
        opts_kwargs["ylim"] = (0, 100)
    if width is not None:
        opts_kwargs["width"] = validate_positive_int(width, name="width")
    if height is not None:
        opts_kwargs["height"] = validate_positive_int(height, name="height")

    return hv.Violin(data.rows, kdims=["group"], vdims=["density"]).opts(**opts_kwargs)

