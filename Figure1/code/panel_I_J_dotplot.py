#!/usr/bin/env python3
"""Unused option for panels I and J (the figure uses panel_I_J_bars.py): a dot plot of the three GO terms.

Dot area = fold enrichment (observed / expected DEGs). Filled = significantly enriched.
Writes data/panel_I_J/Figure1I_J_dotplot.
"""
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from figure1_style import DATA, FS, INK, INK2, MUTED, SURFACE, save

OUT = DATA / "panel_I_J"
TERMS = [("ribosome biogenesis", "Ribosome\nbiogenesis"),
         ("photosynthesis, light harvesting", "Photosynthesis,\nlight harvesting"),
         ("cellular response to cold", "Cellular response\nto cold")]
COLORS = {"Up": "#D55E00", "Down": "#0072B2"}
SIZE = 38  # points^2 per unit of fold enrichment
W, H = 3.25, 2.18  # size in the assembled figure, inches
LEFT, GAP, WIDTH = 1.20, 0.03, 1.01
BOTTOM, HEIGHT = 0.38, 1.50


def draw(ax, d, hours, labels):
    d = d[d["hours_after_onset"] == hours]
    y = 0
    for i, (term, label) in enumerate(TERMS):
        if i:
            ax.axhline(y - 0.5, color=MUTED, lw=0.5, zorder=0, clip_on=False)
        for direction in ("Up", "Down"):
            if labels:
                ax.text(-0.5, y, direction, ha="right", va="center", fontsize=FS - 1, color=INK2)
            for x, treatment in enumerate(("DAWN", "DUSK")):
                r = d[(d["go_name"] == term) & (d["treatment"] == treatment) & (d["direction"] == direction)].iloc[0]
                area = SIZE * r["fold_enrichment"]
                if r["stars"]:
                    ax.scatter(x, y, s=area, color=COLORS[direction], lw=0, zorder=3, clip_on=False)
                    ax.text(x + 0.27, y, f"{r['fold_enrichment']:.1f}×{r['stars']}", ha="left", va="center",
                            fontsize=FS - 2, color=INK, zorder=4)
                else:
                    ax.scatter(x, y, s=area, facecolor=SURFACE, edgecolor=MUTED, lw=0.6, zorder=2, clip_on=False)
            y += 1
        if labels:
            n = int(d[d["go_name"] == term]["go_population_count"].iloc[0])
            ax.text(-1.13, y - 1.5, f"{label}\n({n} genes)", ha="right", va="center", fontsize=FS, color=INK,
                    linespacing=1.05)
    ax.set_xlim(-0.45, 1.85)
    ax.set_ylim(y - 0.5, -0.5)
    ax.set_yticks([])
    for side in ("left", "bottom", "top", "right"):
        ax.spines[side].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=1)
    ax.set_xticks([0, 1], ["DAWN", "DUSK"], fontsize=FS - 1, color=INK2)
    ax.set_title(f"{hours} h after onset", fontweight="bold", fontsize=FS, pad=2, x=0.42)


def add_key(fig, W_):
    handles = [Line2D([], [], marker="o", ls="", color=COLORS["Up"], ms=5, mew=0, label="Up DEGs"),
               Line2D([], [], marker="o", ls="", color=COLORS["Down"], ms=5, mew=0, label="Down DEGs"),
               Line2D([], [], marker="o", ls="", markerfacecolor=SURFACE, markeredgecolor=MUTED, ms=5, mew=0.6,
                      label="Not enriched")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.17 / H),
               borderaxespad=0, borderpad=0, handletextpad=0.1, columnspacing=0.9, fontsize=FS - 1)
    fig.text(0.5, 0.02 / H, "Number, dot area: fold enrichment.  * adj. P < 0.05; *** adj. P < 0.001",
             ha="center", va="bottom", fontsize=FS - 1.5, color=INK2)


def letter(fig, x_in, text):
    fig.text(x_in / W, (H - 0.01) / H, text, fontsize=FS + 2, fontweight="bold", va="top", ha="left", color=INK)


def main():
    d = pd.read_csv(OUT / "selected_go_terms.csv", keep_default_na=False)
    fig = plt.figure(figsize=(W, H), facecolor=SURFACE)
    for j, (panel, hours) in enumerate((("I", 11), ("J", 15))):
        x0 = LEFT + j * (WIDTH + GAP)
        draw(fig.add_axes((x0 / W, BOTTOM / H, WIDTH / W, HEIGHT / H)), d, hours, labels=j == 0)
        letter(fig, 0.0 if j == 0 else x0 - 0.12, panel)
    add_key(fig, W)
    save(fig, "Figure1I_J_dotplot", folder=OUT)


if __name__ == "__main__":
    main()
