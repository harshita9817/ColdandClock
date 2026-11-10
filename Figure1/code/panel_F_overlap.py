#!/usr/bin/env python3
"""Panel F: expression-adjusted same-direction overlap of DAWN and DUSK DEGs."""
import pandas as pd

from figure1_style import DATA, panel_figure, save

COLORS = {"Up": "#D55E00", "Down": "#0072B2"}


def draw(ax, panel="F", ylabel=True):
    d = pd.read_csv(DATA / "panel_F" / "overlap_observed_expected.csv")
    hours = sorted(d["hours_after_onset"].unique())
    ax.axhline(0, color="#8C8C8C", lw=0.6, ls="--")
    y = "dodecile_log2_observed_expected"
    last = d[d["hours_after_onset"] == hours[-1]].set_index("direction")[y]
    nudge = {"Up": 0.0, "Down": 0.0}
    if abs(last["Up"] - last["Down"]) < 0.4:
        high = last.idxmax()
        nudge = {high: 0.18, ("Down" if high == "Up" else "Up"): -0.18}
    for direction, color in COLORS.items():
        s = d[d["direction"] == direction].sort_values("hours_after_onset")
        sig = s["dodecile_bh_adjusted_p"] < 0.05
        ax.plot(s["hours_after_onset"], s[y], color=color, lw=1.0, zorder=2)
        ax.plot(s["hours_after_onset"][sig], s[y][sig], ls="", marker="o", ms=3, color=color, mew=0.8, zorder=3)
        ax.plot(s["hours_after_onset"][~sig], s[y][~sig], ls="", marker="o", ms=3, mfc="white", mec=color,
                mew=0.8, zorder=3)
        ax.text(hours[-1] + 0.6, s[y].iloc[-1] + nudge[direction], direction, color=color, va="center",
                clip_on=False)
    ax.set_xlim(hours[0] - 1, hours[-1] + 0.3)
    ax.set_ylim(-2.15, 2.15)
    ax.set_xticks(hours)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_xlabel("Hours after cold onset", labelpad=2)
    ax.set_ylabel("Overlap $log_2$(O/E)", labelpad=2)


def main():
    fig, ax = panel_figure("F")
    draw(ax)
    save(fig, "Figure1F")


if __name__ == "__main__":
    main()
