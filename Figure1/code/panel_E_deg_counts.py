#!/usr/bin/env python3
"""Panel E: DEG counts at matched durations after cold onset."""
import pandas as pd

from figure1_style import DATA, panel_figure, save

COLORS = {"DAWN up": "#E78AC3", "DAWN down": "#AD1457", "DUSK up": "#66C2A5", "DUSK down": "#006D5B"}


def draw(ax, panel="E", ylabel=True, xlabel=True):
    d = pd.read_csv(DATA / "panel_E" / "deg_counts_matched_duration.csv")
    hours = sorted(d["hours_after_onset"].unique())
    handles = {}
    for direction in ("Up", "Down"):
        s = d[d["direction"] == direction].sort_values("hours_after_onset")
        name = direction.lower()
        handles[f"DAWN {name}"] = ax.plot(s["hours_after_onset"], s["dawn_deg"], color=COLORS[f"DAWN {name}"],
                                          marker="o", ms=2.8, lw=1.0, ls="-", label=f"Dawn {name}")[0]
        handles[f"DUSK {name}"] = ax.plot(s["hours_after_onset"], s["dusk_deg"], color=COLORS[f"DUSK {name}"],
                                          marker="^", ms=3, lw=1.0, ls="--", label=f"Dusk {name}")[0]
    ax.set_xlim(hours[0] - 1, hours[-1] + 0.3)  # same range as panel F
    ax.set_ylim(0, 6000)
    ax.set_xticks(hours)
    ax.set_yticks([0, 2000, 4000, 6000], ["0", "2,000", "4,000", "6,000"])
    if xlabel:
        ax.set_xlabel("Hours after cold onset", labelpad=2)
    ax.set_ylabel("DEGs", labelpad=2)
    ax.legend(handles=[handles[k] for k in COLORS], loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2,
              frameon=False, handlelength=1.6, handletextpad=0.25, columnspacing=1.0, labelspacing=0.15,
              borderaxespad=0.05, borderpad=0)


def main():
    fig, ax = panel_figure("E")
    draw(ax)
    save(fig, "Figure1E")


if __name__ == "__main__":
    main()
