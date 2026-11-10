#!/usr/bin/env python3
"""Earlier panels I and J (replaced by panel_I_J_dotplot.py): enrichment of every tested GO BP term."""
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from figure1_style import DATA, FS, INK, INK2, MUTED, ROOT, SURFACE, panel_figure, save

THRESH = -np.log10(0.05)
HOURS = {"I": 11, "J": 15}
# One marker shape; colour alone separates the groups, drawn in this order.
GROUPS = [("Ribosomal", "#2a78d6", 13), ("Photosynthesis", "#138a5b", 13), ("Light harvesting", "#eda100", 22)]


def draw(ax, panel, ylabel=True):
    folder = DATA / f"panel_{panel}"
    d = pd.read_csv(folder / "go_term_scores.csv")
    counts = pd.read_csv(folder / "term_counts_and_axis_limit.csv").iloc[0]
    lim = counts["axis_limit"]
    x, y = d["DAWN_signed_log10_padj"], d["DUSK_signed_log10_padj"]
    for sx, sy in ((1, -1), (-1, 1)):
        ax.fill_between([0, sx * lim], 0, sy * lim, color="#f4f3ee", zorder=0, lw=0)
    ax.set_xscale("symlog", linthresh=2)
    ax.set_yscale("symlog", linthresh=2)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    for t in (-THRESH, THRESH):
        ax.axvline(t, color=MUTED, lw=0.5, ls=":")
        ax.axhline(t, color=MUTED, lw=0.5, ls=":")
    ax.axvline(0, color=INK2, lw=0.5)
    ax.axhline(0, color=INK2, lw=0.5)
    m = d["plot_group"] == "other BP term"
    ax.scatter(x[m], y[m], s=2.5, color=MUTED, lw=0, zorder=1, rasterized=True)
    for label, color, size in GROUPS:
        m = d["plot_group"] == label
        edge = INK if label == "Light harvesting" else SURFACE
        ax.scatter(x[m], y[m], s=size, marker="o", color=color, edgecolor=edge, lw=0.5, zorder=3)
    ticks = [-50, -10, -2, 0, 2, 10, 50]
    labels = ["50", "10", "2", "0", "2", "10", "50"]
    ax.set_xticks(ticks, labels)
    ax.set_yticks(ticks, labels)
    ax.minorticks_off()
    ax.set_xlabel("DAWN $-log_{10}$(adj. P)\n← Down | Up DEGs →", labelpad=1)
    if ylabel:
        ax.set_ylabel("DUSK $-log_{10}$(adj. P)\n← Down | Up DEGs →", labelpad=1)
    ax.set_title(f"{HOURS[panel]} h after onset", fontweight="bold", pad=3)
    ax.text(0.97, 0.97, f"opposite {int(counts['opposite_direction_terms'])}\nsame {int(counts['same_direction_terms'])}",
            transform=ax.transAxes, ha="right", va="top", fontsize=FS, color=INK2, zorder=6)


def add_key(fig, x_in, y_in, ncol=3, anchor="upper left"):
    """Small colour key, placed by its corner in inches from the figure's lower left."""
    W, H = fig.get_size_inches()
    handles = [Line2D([], [], marker="o", ls="", color=color, ms=4 if label != "Light harvesting" else 4.6,
                      markeredgecolor=INK if label == "Light harvesting" else SURFACE, markeredgewidth=0.4,
                      label=label) for label, color, _ in GROUPS]
    fig.legend(handles=handles, loc=anchor, ncol=ncol, frameon=False, bbox_to_anchor=(x_in / W, y_in / H),
               borderaxespad=0, borderpad=0, handletextpad=0.1, columnspacing=0.8, labelspacing=0.15,
               handlelength=1.0)


def main():
    for panel in "IJ":
        fig, ax = panel_figure(panel)
        draw(ax, panel)
        width, height = fig.get_size_inches()
        add_key(fig, 0.55, height - 0.02, ncol=1)
        save(fig, f"Figure1{panel}", ROOT / "earlier_versions")


if __name__ == "__main__":
    main()
