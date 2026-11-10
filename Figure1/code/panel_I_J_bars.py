#!/usr/bin/env python3
"""Panels I and J: percent of each chosen GO term's genes that are up or down DEGs, as mirrored bars.

The black tick on each bar is the percent expected by chance, that treatment's overall DEG rate.
Writes data/panel_I_J/Figure1I_J_bars (inserted into the figure) and panels/Figure1I, Figure1J.
"""
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import pandas as pd

from figure1_style import DATA, FS, INK, INK2, MUTED, SURFACE, save

OUT = DATA / "panel_I_J"
TERMS = [("ribosome biogenesis", "Ribosome\nbiogenesis"),
         ("photosynthesis, light harvesting", "Photosynthesis,\nlight harvesting"),
         ("cellular response to cold", "Cellular response\nto cold")]
COLORS = {"DAWN": {"Up": "#E78AC3", "Down": "#AD1457"}, "DUSK": {"Up": "#66C2A5", "Down": "#006D5B"}}
LIM = 70
BAR = 0.34
W, H = 3.25, 2.18  # size in the assembled figure, inches
LEFT, GAP, WIDTH = 1.22, 0.03, 1.00
BOTTOM, TOP = 0.56, 0.20  # margins below and above the axes


def draw(ax, d, hours, labels):
    d = d[d["hours_after_onset"] == hours]
    ax.axvline(0, color=INK2, lw=0.6, zorder=2)
    for i, (term, label) in enumerate(TERMS):
        if i:
            ax.axhline(i - 0.5, color=MUTED, lw=0.5, zorder=0)
        for offset, treatment in ((-0.2, "DAWN"), (0.2, "DUSK")):
            y = i + offset
            for direction, sign in (("Up", 1), ("Down", -1)):
                r = d[(d["go_name"] == term) & (d["treatment"] == treatment) & (d["direction"] == direction)].iloc[0]
                value = sign * r["percent_of_term_genes"]
                ax.barh(y, value, height=BAR, color=COLORS[treatment][direction], lw=0, zorder=1)
                expected = sign * r["percent_expected"]
                ax.plot([expected, expected], [y - BAR / 2 - 0.02, y + BAR / 2 + 0.02], color=INK, lw=0.7, zorder=3,
                        solid_capstyle="butt", path_effects=[pe.withStroke(linewidth=1.5, foreground=SURFACE)])
                if r["stars"]:
                    end = max(abs(value), abs(expected))
                    ax.text(sign * (end + 1.5), y + 0.07, r["stars"], ha="left" if sign > 0 else "right",
                            va="center", fontsize=FS, color=INK, zorder=4)
            if labels:
                ax.text(-LIM - 1.5, y, treatment.capitalize(), ha="right", va="center", fontsize=FS - 1, color=INK2)
        if labels:
            n = int(d[d["go_name"] == term]["go_population_count"].iloc[0])
            ax.text(-0.34, i, f"{label}\n({n} genes)", transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=FS, color=INK, linespacing=1.05)
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(len(TERMS) - 0.5, -0.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ticks = [-50, -25, 0, 25, 50]
    ax.set_xticks(ticks, [str(abs(t)) for t in ticks])
    ax.set_xlabel("← Down | Up →", labelpad=1)
    ax.set_title(f"{hours} h after onset", fontweight="bold", fontsize=FS, pad=3)


def add_notes(fig, width, height, center_x, key_y=0.02, title_y=0.15):
    """Shared x-axis title and a one-line key under the axes, positions in inches."""
    fig.text(center_x / width, title_y / height, "% of the term's genes that are DEGs", ha="center", va="bottom",
             fontsize=FS, color=INK)
    fig.text(1.0, key_y / height, "Black tick: expected by chance. Enriched: * adj. P < 0.05; *** adj. P < 0.001",
             ha="right", va="bottom", fontsize=FS - 1.5, color=INK2)


def letter(fig, x_in, height, text, width):
    fig.text(x_in / width, (height - 0.01) / height, text, fontsize=FS + 2, fontweight="bold", va="top", ha="left",
             color=INK)


def main():
    d = pd.read_csv(OUT / "selected_go_terms.csv", keep_default_na=False)
    axes_h = H - BOTTOM - TOP
    # Assembled pair: I holds the term and treatment labels, J shares them.
    fig = plt.figure(figsize=(W, H), facecolor=SURFACE)
    for j, (panel, hours) in enumerate((("I", 11), ("J", 15))):
        x0 = LEFT + j * (WIDTH + GAP)
        draw(fig.add_axes((x0 / W, BOTTOM / H, WIDTH / W, axes_h / H)), d, hours, labels=j == 0)
        letter(fig, 0.0 if j == 0 else x0 - 0.02, H, panel, W)
    add_notes(fig, W, H, LEFT + WIDTH + GAP / 2)
    save(fig, "Figure1I_J_bars", folder=OUT)
    # Stand-alone panels, each with its own labels.
    for panel, hours in (("I", 11), ("J", 15)):
        w = W - WIDTH - GAP
        fig = plt.figure(figsize=(w, H), facecolor=SURFACE)
        draw(fig.add_axes((LEFT / w, BOTTOM / H, WIDTH / w, axes_h / H)), d, hours, labels=True)
        letter(fig, 0.0, H, panel, w)
        fig.text((LEFT + WIDTH / 2) / w, 0.15 / H, "% of term's genes", ha="center", va="bottom", fontsize=FS, color=INK)
        fig.text(0.0, 0.02 / H, "Tick: expected. * adj. P < 0.05; *** adj. P < 0.001", ha="left", va="bottom",
                 fontsize=FS - 1.5, color=INK2)
        save(fig, f"Figure1{panel}")


if __name__ == "__main__":
    main()
