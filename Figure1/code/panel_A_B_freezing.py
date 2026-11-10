#!/usr/bin/env python3
"""Panels A and B: freezing injury curves with per-set LT50 insets."""
import numpy as np
import pandas as pd

from figure1_style import DATA, add_axes, panel_figure, save

STYLES = {
    "Sorghum dawn control": dict(color="#AD1457", marker="o", mfc="white", label="Control Dawn"),
    "Sorghum dawn cold day 1": dict(color="#E78AC3", marker="^", mfc="#E78AC3", label="Cold Dawn"),
    "Sorghum dusk control": dict(color="#1D2557", marker="s", mfc="white", label="Control Dusk"),
    "Sorghum dusk cold day 1": dict(color="#2A9DA6", marker="o", mfc="white", label="Cold Dusk"),
}


def logistic(x, top, bottom, midpoint, k):
    return bottom + (top - bottom) / (1 + np.exp(k * (x - midpoint)))


def draw(ax, panel, ylabel=True):
    folder = DATA / f"panel_{panel}"
    points = pd.read_csv(folder / "points_mean_sem.csv")
    fits = pd.read_csv(folder / "curve_fit_parameters.csv")
    lt50 = pd.read_csv(folder / "inset_lt50_per_set.csv")
    stars = pd.read_csv(folder / "inset_statistics.csv")["Significance"].iloc[0]
    groups = list(fits["group"])

    ax.set_xlim(1, -4.5)
    ax.set_ylim(-10, 120)
    ax.set_xticks([0, -1, -2, -3, -4])
    ax.set_yticks(range(0, 101, 20))
    ax.set_xlabel("Temperature (°C)", labelpad=2)
    if ylabel:
        ax.set_ylabel("Relative injury (%)", labelpad=2)
    x = np.linspace(0.3, -4.3, 300)
    for _, fit in fits.iterrows():
        s = STYLES[fit["group"]]
        ax.plot(x, logistic(x, fit["top"], fit["bottom"], fit["midpoint"], fit["k"]),
                color=s["color"], lw=1.1, zorder=2)
    for group in groups:
        s = STYLES[group]
        d = points[points["group"] == group]
        ax.errorbar(d["temperature"], d["mean"], yerr=d["sem"], ls="", color=s["color"],
                    elinewidth=0.6, capsize=1.2, capthick=0.6, zorder=3)
        ax.plot(d["temperature"], d["mean"], ls="", marker=s["marker"], mfc=s["mfc"],
                mec=s["color"], mew=0.8, ms=3.4, zorder=4)
    handles = [ax.plot([], [], color=STYLES[g]["color"], lw=1.1, marker=STYLES[g]["marker"],
                       mfc=STYLES[g]["mfc"], mec=STYLES[g]["color"], mew=0.8, ms=3.4,
                       label=STYLES[g]["label"])[0] for g in groups]
    ax.legend(handles=handles, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.98), ncol=2,
              handlelength=1.8, handletextpad=0.4, columnspacing=1.5, borderaxespad=0.05, borderpad=0)

    inset = ax.inset_axes([0.19, 0.46, 0.17, 0.50])
    values = lt50["lt50_C"].to_numpy()
    span = values.max() - values.min()
    limits = (values.min() - 0.1 * span, values.max() + 0.4 * span)
    step = 0.1 if span <= 0.3 else 0.2 if span <= 0.6 else 0.5
    ticks = np.arange(np.ceil(limits[0] / step), np.floor(limits[1] / step) + 1) * step
    for i, group in enumerate(groups, start=1):
        color = STYLES[group]["color"]
        low, q1, med, q3, high = np.percentile(lt50.loc[lt50["group"] == group, "lt50_C"], [0, 25, 50, 75, 100])
        inset.plot([i, i], [low, q1], color=color, lw=0.6)
        inset.plot([i, i], [q3, high], color=color, lw=0.6)
        inset.plot([i - 0.2, i + 0.2], [low, low], color=color, lw=0.6)
        inset.plot([i - 0.2, i + 0.2], [high, high], color=color, lw=0.6)
        inset.fill_between([i - 0.35, i + 0.35], q1, q3, facecolor=color, alpha=0.35, lw=0)
        inset.plot([i - 0.35, i + 0.35, i + 0.35, i - 0.35, i - 0.35], [q1, q1, q3, q3, q1], color=color, lw=0.6)
        inset.plot([i - 0.35, i + 0.35], [med, med], color=color, lw=1.0)
    inset.set_xlim(0.4, 2.6)
    inset.set_ylim(*limits)
    inset.set_xticks([1, 2], ["", ""])
    inset.set_yticks(ticks, [f"{t:.1f}".replace("-", "−") for t in ticks])
    inset.set_ylabel("$LT_{50}$ (°C)", labelpad=1)
    inset.patch.set_alpha(0)
    inset.text(1.5, limits[1], stars, ha="center", va="top")


def main():
    for panel in "AB":
        fig, ax = panel_figure(panel)
        draw(ax, panel)
        save(fig, f"Figure1{panel}")


if __name__ == "__main__":
    main()
