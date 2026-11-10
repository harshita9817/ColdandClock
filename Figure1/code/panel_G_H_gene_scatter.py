#!/usr/bin/env python3
"""Panels G and H: log2(DAWN/CONT) vs log2(DUSK/CONT) of every tested gene, against the time-matched control."""
import numpy as np
import pandas as pd

from figure1_style import DATA, FS, INK2, SURFACE, panel_figure, save

ONE_C, NOT_C = "#a9a8a2", "#e6e5e0"
# Colours of the four DEG groups, with a darker shade of each for its count.
PALETTE = {"DAWN up / DUSK down": ("#e9a23b", "#ba822f"), "DAWN down / DUSK up": ("#2a9d8f", "#227e72"),
           "both up": ("#b8324f", "#93283f"), "both down": ("#3f6fc4", "#32599d")}
LIM = 8
HOURS = {"G": 11, "H": 15}


def corner(ax, text, x, y, ha, va, color):
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=FS, fontweight="bold", color=color,
            zorder=6, bbox=dict(boxstyle="round,pad=0.12", fc=SURFACE, ec="none", alpha=0.85))


def draw(ax, panel, ylabel=True):
    d = pd.read_csv(DATA / f"panel_{panel}" / "gene_log2fc_and_category.csv")
    x, y = np.clip(d["DAWN_log2FC"], -LIM, LIM), np.clip(d["DUSK_log2FC"], -LIM, LIM)
    style = [("not DEG", NOT_C, 1.0), ("DEG in one treatment", ONE_C, 1.5)]
    style += [(name, PALETTE[name][0], 2.4) for name in ("both up", "both down", "DAWN up / DUSK down", "DAWN down / DUSK up")]
    for name, color, size in style:
        m = d["category"] == name
        ax.scatter(x[m], y[m], s=size, color=color, lw=0, zorder=2, rasterized=True)
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-LIM, LIM)
    ax.axhline(0, color=INK2, lw=0.5, zorder=1)
    ax.axvline(0, color=INK2, lw=0.5, zorder=1)
    ax.set_xticks([-8, -4, 0, 4, 8])
    ax.set_yticks([-8, -4, 0, 4, 8])
    ax.set_xlabel("$log_2$(Dawn/Control)", labelpad=1)
    if ylabel:
        ax.set_ylabel("$log_2$(Dusk/Control)", labelpad=3)
    ax.set_title(f"{HOURS[panel]} h after onset", fontweight="bold", pad=3)
    n = d["category"].value_counts()
    corner(ax, f"{n['DAWN down / DUSK up']:,}", 0.03, 0.97, "left", "top", PALETTE["DAWN down / DUSK up"][1])
    corner(ax, f"{n['DAWN up / DUSK down']:,}", 0.97, 0.03, "right", "bottom", PALETTE["DAWN up / DUSK down"][1])
    corner(ax, f"{n['both up']:,}", 0.97, 0.97, "right", "top", PALETTE["both up"][1])
    corner(ax, f"{n['both down']:,}", 0.03, 0.03, "left", "bottom", PALETTE["both down"][1])


def main():
    for panel in "GH":
        fig, ax = panel_figure(panel)
        draw(ax, panel)
        save(fig, f"Figure1{panel}")


if __name__ == "__main__":
    main()
