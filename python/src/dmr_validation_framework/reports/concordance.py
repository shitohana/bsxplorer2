"""Pairwise caller concordance from the UpSet membership matrix.

For every pair of callers we compute the Jaccard index over the consensus
regions they support: ``|A and B| / |A or B|``. This answers "which callers
tend to call the same loci" - a compact, functional companion to the UpSet
plot. Shared by the matplotlib bundle and the interactive report.
"""

from __future__ import annotations

import pandas as pd

META_COLUMNS = {"region_id", "context"}


def _bool_columns(membership: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    set_cols = [c for c in membership.columns if c not in META_COLUMNS]
    if not set_cols:
        return pd.DataFrame(), []
    matrix = membership[set_cols].apply(lambda s: pd.to_numeric(s, errors="coerce").fillna(0) > 0)
    return matrix, set_cols


def compute_caller_concordance(membership: pd.DataFrame) -> pd.DataFrame:
    """Long-form pairwise Jaccard table: caller_a, caller_b, jaccard, n_intersect."""
    matrix, callers = _bool_columns(membership)
    if not callers:
        return pd.DataFrame()
    rows: list[dict] = []
    for a in callers:
        for b in callers:
            inter = int((matrix[a] & matrix[b]).sum())
            union = int((matrix[a] | matrix[b]).sum())
            rows.append(
                {
                    "caller_a": a,
                    "caller_b": b,
                    "jaccard": (inter / union) if union else float("nan"),
                    "n_intersect": inter,
                }
            )
    return pd.DataFrame(rows)


def build_concordance_figure(long_df: pd.DataFrame, *, title: str = "Caller concordance (Jaccard)"):
    """Build (do not save) a matplotlib heatmap of the concordance matrix."""
    if long_df.empty:
        return None
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    callers = list(dict.fromkeys(long_df["caller_a"].tolist()))
    grid = long_df.pivot(index="caller_a", columns="caller_b", values="jaccard").reindex(
        index=callers, columns=callers
    )
    counts = long_df.pivot(index="caller_a", columns="caller_b", values="n_intersect").reindex(
        index=callers, columns=callers
    )
    values = grid.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(max(5, len(callers) * 0.9), max(4, len(callers) * 0.8)))
    im = ax.imshow(values, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(callers)))
    ax.set_yticks(range(len(callers)))
    ax.set_xticklabels(callers, rotation=45, ha="right")
    ax.set_yticklabels(callers)
    for i in range(len(callers)):
        for j in range(len(callers)):
            jac = values[i, j]
            if np.isfinite(jac):
                n = counts.to_numpy()[i, j]
                ax.text(
                    j,
                    i,
                    f"{jac:.2f}\n(n={int(n)})",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="white" if jac < 0.6 else "black",
                )
    fig.colorbar(im, ax=ax, label="Jaccard overlap")
    ax.set_title(title)
    fig.tight_layout()
    return fig
