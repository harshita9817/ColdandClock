#!/usr/bin/env python3
"""Panels C and D: DEG state transitions across all tested genes."""
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Polygon, Rectangle

from figure1_style import DATA, FS, panel_figure, save

STATES = ["Down", "Not DEG", "Up"]
COLORS = {"C": {"Up": "#E78AC3", "Down": "#AD1457", "Not DEG": "#D9D9D9"},
          "D": {"Up": "#66C2A5", "Down": "#006D5B", "Not DEG": "#D9D9D9"}}
TITLE = {"C": "Dawn vs Control", "D": "Dusk vs Control"}
BAR = 0.12


def bezier(p0, p1, p2, p3, n=40):
    t = np.linspace(0, 1, n)[:, None]
    p0, p1, p2, p3 = map(np.asarray, (p0, p1, p2, p3))
    return (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3


def ribbon(ax, x0, x1, y0a, y0b, y1a, y1b, color, alpha):
    dx = x1 - x0
    top = bezier((x0, y0a), (x0 + dx * 0.45, y0a), (x1 - dx * 0.45, y1a), (x1, y1a))
    bottom = bezier((x1, y1b), (x1 - dx * 0.45, y1b), (x0 + dx * 0.45, y0b), (x0, y0b))
    ax.add_patch(Polygon(np.vstack([top, bottom]), closed=True, facecolor=color, alpha=alpha, lw=0))


def draw(ax, panel, ylabel=True):
    folder = DATA / f"panel_{panel}"
    strata = pd.read_csv(folder / "strata_counts.csv")
    flows = pd.read_csv(folder / "flow_counts.csv")
    colors = COLORS[panel]
    times = sorted(strata["time_h"].unique())
    n_genes = int(strata[strata["time_h"] == times[0]]["genes"].sum())
    gap = max(n_genes * 0.015, 25)
    xpos = {t: i + 1 for i, t in enumerate(times)}

    spans = {}
    for t in times:
        y = 0
        for state in STATES:
            g = int(strata[(strata["time_h"] == t) & (strata["state"] == state)]["genes"].iloc[0])
            spans[t, state] = (y, y + g)
            y += g + gap

    for left, right in zip(times, times[1:]):
        lo = {s: spans[left, s][0] for s in STATES}
        ro = {s: spans[right, s][0] for s in STATES}
        for _, f in flows[(flows["from_time_h"] == left) & (flows["genes"] > 0)].iterrows():
            s, d, g = f["from_state"], f["to_state"], f["genes"]
            ribbon(ax, xpos[left] + BAR / 2, xpos[right] - BAR / 2, lo[s], lo[s] + g, ro[d], ro[d] + g,
                   colors[s], 0.16 if s == "Not DEG" else 0.31)
            lo[s] += g
            ro[d] += g

    halo = [pe.withStroke(linewidth=1.2, foreground="white")]
    for t in times:
        for state in STATES:
            y0, y1 = spans[t, state]
            ax.add_patch(Rectangle((xpos[t] - BAR / 2, y0), BAR, y1 - y0, facecolor=colors[state],
                                   edgecolor="white", lw=0.3, zorder=3))
            middle = y1 - y0 >= n_genes * 0.09
            ax.text(xpos[t], (y0 + y1) / 2 if middle else y1 + gap * 0.18, f"{int(round(y1 - y0))}",
                    ha="center", va="center" if middle else "bottom", fontsize=FS - 1,
                    fontweight="bold", zorder=4, path_effects=halo)

    ax.set_xlim(0.65, len(times) + 0.35)
    ax.set_ylim(0, n_genes + gap * 2 + n_genes * 0.085)  # room for counts above the top bars
    # DUSK samples were collected 12 h later; the axis shows hours after cold onset for both.
    ax.set_xticks(list(xpos.values()), [t - times[0] + 3 for t in times])
    ax.set_xlabel("Hours after cold onset", labelpad=2)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.text(0.05, 1.015, TITLE[panel], transform=ax.transAxes, fontweight="bold", va="bottom", ha="left")
    handles = [Patch(facecolor=colors[s], label=s) for s in ("Up", "Not DEG", "Down")]
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 0.985), ncol=3, frameon=False,
              handlelength=0.9, handleheight=0.9, handletextpad=0.35, columnspacing=0.9, borderaxespad=0.05, borderpad=0)


def main():
    for panel in "CD":
        fig, ax = panel_figure(panel)
        draw(ax, panel)
        save(fig, f"Figure1{panel}")


if __name__ == "__main__":
    main()
