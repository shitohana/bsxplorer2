"""Single source of truth for report colors and legend placement.

Every figure in the validation reports (matplotlib static bundle, validation
summary, interactive HoloViews report) should pull colors from here so that the
same entity - a methylation context, an effect direction, a model-agreement
status, a confidence class, a DMR tier - is always drawn with the same color.

The base palette is Okabe-Ito, which stays distinguishable in grayscale print
and for the common forms of color-vision deficiency.
"""

from __future__ import annotations

import colorsys

# Okabe-Ito qualitative palette (color-vision-deficiency safe).
OKABE_ITO = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # green
    "#D55E00",  # vermillion
    "#CC79A7",  # purple
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
]

PALETTE = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "skyblue": "#56B4E9",
    "yellow": "#F0E442",
    "gray": "#6B7280",
    "lightgray": "#E5E7EB",
    "black": "#111111",
}

# Methylation contexts.
CONTEXT_COLORS = {
    "CG": PALETTE["blue"],
    "CHG": PALETTE["orange"],
    "CHH": PALETTE["green"],
    "NA": PALETTE["gray"],
}

# Effect direction.
DIRECTION_COLORS = {
    "hyper": PALETTE["red"],
    "hypo": PALETTE["blue"],
    "near_zero": PALETTE["gray"],
    "unknown": PALETTE["gray"],
    "mixed": PALETTE["purple"],
}

# Annotation seqname compatibility.
MATCH_COLORS = {
    "matched": PALETTE["green"],
    "unmatched": PALETTE["red"],
}

# GLM/GLMM model-agreement status.
STATUS_COLORS = {
    "GLMM_CONFIRMED": PALETTE["green"],
    "GLM_ONLY": PALETTE["orange"],
    "GLMM_ONLY": PALETTE["skyblue"],
    "DIRECTION_DISCORDANT": PALETTE["red"],
    "NOT_SIGNIFICANT": PALETTE["gray"],
    "INSUFFICIENT_DATA": PALETTE["lightgray"],
    "MODEL_FAILED": PALETTE["black"],
}

# Final confidence class.
CONFIDENCE_COLORS = {
    "HIGH_CONFIDENCE": PALETTE["green"],
    "MODERATE_CONFIDENCE": PALETTE["skyblue"],
    "WEAK_CONFIDENCE": PALETTE["orange"],
    "LOW_CONFIDENCE": PALETTE["red"],
    "NOT_EVALUATED": PALETTE["gray"],
}

# Final DMR tier.
TIER_COLORS = {
    "TIER_1_CONSENSUS_VALIDATED": PALETTE["green"],
    "TIER_2_SINGLE_CALLER_VALIDATED": PALETTE["skyblue"],
    "TIER_3_MULTI_CALLER_QC_LIMITED": PALETTE["orange"],
    "TIER_4_CALLER_DISCORDANT": PALETTE["purple"],
    "TIER_5_REJECTED_BY_AUDIT": PALETTE["red"],
    "CANDIDATE_ONLY": PALETTE["gray"],
}


def caller_color(name: str) -> str:
    """Deterministic, stable color for an arbitrary caller name.

    Known callers map to fixed Okabe-Ito entries; unknown callers get a
    reproducible hue derived from the name so the same caller is always the
    same color across figures and runs.
    """
    fixed = {
        "DSS": OKABE_ITO[0],
        "methylKit": OKABE_ITO[1],
        "BSmooth": OKABE_ITO[2],
        "comb-p": OKABE_ITO[3],
        "metilene": OKABE_ITO[4],
        "dmrseq": OKABE_ITO[5],
        "internal": PALETTE["gray"],
        "Regional Evidence": PALETTE["gray"],
    }
    if name in fixed:
        return fixed[name]
    digest = sum(ord(ch) for ch in str(name))
    hue = (digest % 360) / 360.0
    r, g, b = colorsys.hls_to_rgb(hue, 0.45, 0.55)
    return f"#{int(r * 255):02X}{int(g * 255):02X}{int(b * 255):02X}"


def caller_cmap(names) -> dict[str, str]:
    """Stable {caller: color} map for a set of caller names (dedup, ordered)."""
    return {str(name): caller_color(str(name)) for name in dict.fromkeys(names)}


def color_for(value: object, mapping: dict[str, str], default: str | None = None) -> str:
    """Look up a color, falling back to a neutral gray for unmapped values."""
    return mapping.get(str(value), default or PALETTE["gray"])


def legend_outside(ax, *, title: str | None = None, where: str = "right", ncol: int = 1):
    """Place the legend outside the axes so it never overlaps the data.

    ``where`` is ``"right"`` (default) or ``"bottom"``. Callers must save the
    figure with ``bbox_inches="tight"`` so the externalized legend is not
    clipped (the report helpers already do this).
    """
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None
    if where == "bottom":
        return ax.legend(
            handles,
            labels,
            title=title,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.18),
            ncol=ncol or max(1, len(labels)),
            frameon=False,
        )
    return ax.legend(
        handles,
        labels,
        title=title,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        ncol=ncol,
        frameon=False,
    )
