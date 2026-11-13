"""The closing-pattern figure: bonds coloured by their coupling change at the
closing pattern (red raised, blue lowered), circles sized by the weight of the
closing orbital on each atom."""
from __future__ import annotations

import numpy as np

RED, BLUE, INK = "#e34948", "#2a78d6", "#1f2937"


def _cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("gapclosure", [BLUE, "#f0efec", RED])


def _horizontal(P):
    P = P - P.mean(0)
    if len(P) < 2:
        return P
    _, _, Vt = np.linalg.svd(P, full_matrices=False)
    Q = P @ Vt.T
    return Q[:, :2] if Q.shape[1] >= 2 else np.c_[Q, np.zeros(len(Q))]


def plot_closing(result, ax=None, *, title=True, colorbar=True, orbital=True,
                 labels=False, scale=None, horizontal=True, vmax=None, lw=None):
    """Draw the closing pattern of a ClosureResult. Returns the matplotlib Axes."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    sk = result.skeleton
    P = np.asarray(sk.layout(), float)
    if horizontal:
        P = _horizontal(P)
    span = np.ptp(P, axis=0) + 1.4
    if ax is None:
        w = min(9.0, max(3.0, 0.55 * span[0]))
        h = max(2.2, w * span[1] / span[0]) + (0.7 if colorbar else 0)
        fig, ax = plt.subplots(figsize=(w, h))
    fig = ax.figure
    rho = result.rho
    vmax = vmax or max(0.05, round(min(rho, 1.0) + 0.005, 2))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    cmap = _cmap()
    # bond width in points: about 1/12 of a bond length on screen
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    pts_per_unit = 72 * min(bbox.width / span[0], bbox.height / span[1])
    lw = lw or max(1.0, min(6.0, pts_per_unit / 9))
    for (u, v), s in result.pattern.items():
        c = cmap(norm(s * rho)) if rho < 1 else "#9aa3ad"
        ax.plot(*P[[u, v]].T, color=c, lw=lw, solid_capstyle="round", zorder=2)
    ax.scatter(P[:, 0], P[:, 1], s=(lw * 0.9) ** 2, color=INK, zorder=3, lw=0)
    if orbital and result.orbital is not None:
        area = (0.9 * pts_per_unit) ** 2 * 2.2          # weight 0.15 ~ circle of 0.6 bond lengths
        ax.scatter(P[:, 0], P[:, 1], s=area * result.orbital, facecolor="none",
                   edgecolor=INK, lw=0.7, zorder=4)
    if labels:
        for a, (x, y) in enumerate(P):
            ax.text(x, y, str(a), fontsize=6, ha="center", va="center", zorder=5,
                    bbox=dict(boxstyle="circle,pad=0.1", fc="white", lw=0))
    ax.set_aspect("equal")
    ax.set_xlim(P[:, 0].min() - 0.7, P[:, 0].max() + 0.7)
    ax.set_ylim(P[:, 1].min() - 0.7, P[:, 1].max() + 0.7)
    ax.axis("off")
    if title:
        nm = sk.name or "skeleton"
        txt = f"{nm},  ρ = {rho:.4f}" if result.exact else \
            f"{nm},  {result.lower:.4f} ≤ ρ ≤ {result.upper:.4f}"
        ax.set_title(txt, fontsize=10, loc="left")
    if colorbar and rho < 1:
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.06, pad=0.02, shrink=0.6)
        cb.set_label("relative change of the coupling (red raised, blue lowered)", fontsize=8)
        cb.ax.tick_params(labelsize=7)
        cb.outline.set_linewidth(0.4)
    return ax
