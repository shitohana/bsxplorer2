"""Candidate-selection funnel shared by the static and interactive reports.

The funnel narrates how DMR candidates are progressively filtered:

    all DMR calls -> consensus regions -> multi-caller (>=2) ->
    direction-consistent -> not rejected by audit -> GLMM-confirmed

Stages 2-5 are a strict nested subset chain over the consensus region set, so
the bars genuinely shrink. The first ("all calls") and last ("GLMM-confirmed")
stages are evaluated on different region sets and are flagged ``nested=False``
so the renderer can mark them and the caption can stay honest.

``compute_funnel_stages`` is pure pandas (no plotting deps) and is reused by
both the matplotlib bundle and the HoloViews interactive report.
"""

from __future__ import annotations

import pandas as pd


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "t", "y"}


def _num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def compute_funnel_stages(
    consensus: pd.DataFrame | None = None,
    glm_glmm: pd.DataFrame | None = None,
    caller_rows: pd.DataFrame | None = None,
) -> list[dict]:
    """Ordered funnel stages: ``label``, ``n``, ``nested``, ``note``."""
    stages: list[dict] = []

    if caller_rows is not None and not caller_rows.empty:
        stages.append(
            {
                "label": "All DMR calls",
                "n": int(len(caller_rows)),
                "nested": False,
                "note": "raw calls across callers, before clustering",
            }
        )

    if consensus is not None and not consensus.empty:
        stages.append(
            {
                "label": "Consensus regions",
                "n": int(len(consensus)),
                "nested": True,
                "note": "overlapping calls clustered into loci",
            }
        )
        # Thread a single running mask through stages 3-5 so each stage is a
        # strict subset of the previous one. Computing any stage against the
        # full consensus set (instead of the running subset) breaks the funnel:
        # e.g. caller-discordant regions dropped at "Direction-consistent" would
        # reappear under "Not rejected by audit", making the bars grow again.
        running = pd.Series(True, index=consensus.index)
        if "n_callers_supporting" in consensus.columns:
            running = running & (_num(consensus["n_callers_supporting"]) >= 2)
            stages.append(
                {
                    "label": "Multi-caller (>=2)",
                    "n": int(running.sum()),
                    "nested": True,
                    "note": "supported by at least two callers",
                }
            )
        # Direction-consistency and the audit filter are merged into a single
        # stage: on this tier scheme the only multi-caller regions the audit
        # rejects are exactly the direction-discordant ones, so a separate
        # "audit" bar would just duplicate the direction-consistent count. We
        # apply both filters and emit one stage, labelled by whatever applies.
        direction_present = "direction_consensus" in consensus.columns
        tier_present = "final_dmr_tier" in consensus.columns
        if direction_present:
            running = running & (consensus["direction_consensus"].astype(str) == "same")
        if tier_present:
            rejected = {"TIER_5_REJECTED_BY_AUDIT", "CANDIDATE_ONLY"}
            running = running & (~consensus["final_dmr_tier"].astype(str).isin(rejected))
        if direction_present and tier_present:
            label, note = (
                "Direction-consistent & audit-passed",
                "consistent hyper/hypo direction across callers and not rejected by audit (TIER_5/candidate-only excluded)",
            )
        elif direction_present:
            label, note = "Direction-consistent", "consistent hyper/hypo direction across callers"
        elif tier_present:
            label, note = "Not rejected by audit", "excludes TIER_5 and candidate-only regions"
        else:
            label = None
        if label is not None:
            stages.append({"label": label, "n": int(running.sum()), "nested": True, "note": note})

    if glm_glmm is not None and not glm_glmm.empty and "confirmed_by_glmm" in glm_glmm.columns:
        confirmed = int(glm_glmm["confirmed_by_glmm"].map(_truthy).sum())
        stages.append(
            {
                "label": "GLMM-confirmed",
                "n": confirmed,
                "nested": False,
                "note": "confirmatory CpG-level GLMM on selected top-N candidates (separate set)",
            }
        )

    return stages


def build_funnel_figure(stages: list[dict], *, title: str = "Candidate selection funnel"):
    """Build (do not save) the funnel as centered trapezoid steps (matplotlib).

    Returns the figure (caller is responsible for saving/closing) or ``None``.
    """
    if not stages:
        return None
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Rectangle

    from dmr_validation_framework.reports.palette import PALETTE

    counts = [max(int(stage["n"]), 0) for stage in stages]
    labels = [f"{i + 1}. {stage['label']}" for i, stage in enumerate(stages)]
    colors = [
        PALETTE["blue"]
        if stage["nested"]
        else (PALETTE["green"] if stage["label"] == "GLMM-confirmed" else PALETTE["gray"])
        for stage in stages
    ]
    n = len(stages)
    max_c = max(counts) or 1
    band, gap = 1.0, 0.45
    # Stage 1 at the top: assign the largest y to the first stage.
    centers = [(n - 1 - i) * (band + gap) for i in range(n)]
    half = [c / 2 for c in counts]

    fig, ax = plt.subplots(figsize=(8.5, max(3.5, n * 0.95)))
    # Connectors (trapezoids) between consecutive stages, drawn first.
    for i in range(n - 1):
        y_bottom_i = centers[i] - band / 2
        y_top_next = centers[i + 1] + band / 2
        ax.add_patch(
            Polygon(
                [
                    (-half[i], y_bottom_i),
                    (half[i], y_bottom_i),
                    (half[i + 1], y_top_next),
                    (-half[i + 1], y_top_next),
                ],
                closed=True,
                facecolor=PALETTE["lightgray"],
                edgecolor="none",
                alpha=0.6,
            )
        )
    # Stage blocks. Narrow bars cannot hold their label inside (white text would
    # spill onto the white background and look truncated), so when the label is
    # wider than the bar we keep the count inside and place the label to the
    # right in dark text. ``char_w`` is an approximate per-character width in data
    # units derived from the axis span.
    char_w = 0.010 * max_c
    for center, count, label, color in zip(centers, counts, labels, colors):
        ax.add_patch(
            Rectangle((-count / 2, center - band / 2), count, band, facecolor=color, edgecolor="white", linewidth=1.0)
        )
        if len(label) * char_w <= count:
            ax.text(0, center, f"{label}\n{count:,}", ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        else:
            ax.text(0, center, f"{count:,}", ha="center", va="center", fontsize=9, color="white", fontweight="bold")
            ax.text(
                count / 2 + max_c * 0.012,
                center,
                label,
                ha="left",
                va="center",
                fontsize=9,
                color=PALETTE["gray"],
                fontweight="bold",
            )
    ax.set_xlim(-max_c * 0.62, max_c * 0.62)
    ax.set_ylim(min(centers) - band, max(centers) + band)
    ax.axis("off")
    ax.set_title(title)
    fig.tight_layout()
    return fig


def render_funnel_matplotlib(stages: list[dict], out_path, *, title: str = "Candidate selection funnel") -> bool:
    """Render the funnel to a file. Returns success."""
    fig = build_funnel_figure(stages, title=title)
    if fig is None:
        return False
    import matplotlib.pyplot as plt
    from pathlib import Path

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return True
