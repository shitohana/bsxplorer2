from __future__ import annotations

from collections.abc import Mapping, Sequence

from beartype.typing import Optional

from bsx2 import RegionReader

from .data import DiscreteRegionData
from .metagene import (
    AnnotProfileLayout,
    MetageneProfileSegment,
    collect_layout_parts_from_hcannot,
    compose_layout_drd,
    compute_discrete_regions,
    validate_segments,
)

PartInputMap = Mapping[str, tuple[Sequence[object], Sequence[str]]]


def _resolve_manual_layout(
    *,
    layout: AnnotProfileLayout | None,
    segments: Optional[Sequence[MetageneProfileSegment]],
    parts_order: Optional[Sequence[str]],
) -> tuple[list[str], list[MetageneProfileSegment]]:
    if layout is not None:
        resolved_parts = [part.name for part in layout.parts]
        resolved_segments = list(layout.segments if segments is None else segments)
    else:
        if segments is None or parts_order is None:
            raise ValueError(
                "manual metagene composition requires either layout=... or both "
                "segments=... and parts_order=..."
            )
        resolved_parts = list(parts_order)
        resolved_segments = list(segments)

    validate_segments(resolved_segments)
    if len(resolved_parts) != len(resolved_segments):
        raise ValueError("parts_order and segments must contain the same number of items")
    return resolved_parts, resolved_segments


def build_manual_metagene(
    reader: RegionReader,
    *,
    part_map: PartInputMap,
    layout: AnnotProfileLayout | None = None,
    segments: Optional[Sequence[MetageneProfileSegment]] = None,
    parts_order: Optional[Sequence[str]] = None,
    reverse_negative: bool = True,
) -> DiscreteRegionData:
    """
    Build one composed metagene from explicitly prepared per-part contigs.

    Parameters
    ----------
    reader
        `RegionReader`-like object used to query methylation data.
    part_map
        Mapping from part name to ``(contigs, labels)``. Labels should align
        across parts for the same biological feature.
    layout
        Optional ordered annotation layout. When provided, the final composed
        segments and part order are derived from it unless overridden via
        ``segments``.
    segments
        Optional final metagene segments. Required together with ``parts_order``
        when ``layout`` is omitted.
    parts_order
        Ordered list of part names to compose. Required when ``layout`` is
        omitted.
    reverse_negative
        Whether negative-strand contigs should be reversed before insertion.

    Returns
    -------
    DiscreteRegionData
        One composed metagene-compatible profile collection that can be passed
        directly to the standard ``bsx2.viz`` render functions.

    Notes
    -----
    This is the low-level layout entrypoint. Use it when the biological parts
    are already known or prepared externally and you want full control over the
    composition order.
    """

    resolved_parts, resolved_segments = _resolve_manual_layout(
        layout=layout,
        segments=segments,
        parts_order=parts_order,
    )

    drd_map: dict[str, DiscreteRegionData] = {}
    for part_name, segment in zip(resolved_parts, resolved_segments):
        contigs, labels = part_map.get(part_name, ([], []))
        if not contigs:
            continue
        drd_map[part_name] = compute_discrete_regions(
            reader,
            contigs,
            segments=[segment],
            reverse_negative=reverse_negative,
            labels=list(labels),
        )

    if not drd_map:
        return DiscreteRegionData()

    return compose_layout_drd(
        drd_map,
        segments=resolved_segments,
        parts_order=resolved_parts,
    )


def build_annotation_metagene(
    reader: RegionReader,
    annot: object,
    *,
    layout: AnnotProfileLayout,
    segments: Optional[Sequence[MetageneProfileSegment]] = None,
    reverse_negative: bool = True,
    limit: Optional[int] = None,
) -> DiscreteRegionData:
    """
    Build one composed metagene from annotation and a declared layout.

    This is the high-level metagene entrypoint for AW25-style
    ``RegionReader + HcAnnotStore`` workflows.

    Parameters
    ----------
    reader
        Reader used to query methylation values from ``report.bsx``.
    annot
        Annotation store compatible with ``collect_layout_parts_from_hcannot``.
    layout
        Ordered layout describing the biological parts to resolve per gene.
    segments
        Optional explicit rendering segments. Defaults to the layout-derived
        segments.
    reverse_negative
        Whether negative-strand features should be aligned to the shared
        biological orientation.
    limit
        Optional cap on the number of genes collected from the annotation.

    Returns
    -------
    DiscreteRegionData
        Composed metagene-compatible profiles ready for line, heatmap, box,
        violin, or downstream metagene renderers.
    """

    part_map = collect_layout_parts_from_hcannot(
        annot,
        layout=layout,
        limit=limit,
    )
    return build_manual_metagene(
        reader,
        part_map=part_map,
        layout=layout,
        segments=segments,
        reverse_negative=reverse_negative,
    )
