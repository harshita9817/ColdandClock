"""Shared paths, fonts, colours and panel layout for Figure 1."""
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "source"
DATA = ROOT / "data"
PANELS = ROOT / "panels"

FS = 8
INK, INK2, MUTED, SURFACE = "#0b0b0b", "#52514e", "#c9c8c2", "#ffffff"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "Nimbus Sans"],
    "font.size": FS, "axes.labelsize": FS, "axes.titlesize": FS,
    "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": FS,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
    "xtick.major.pad": 2, "ytick.major.pad": 2,
    "axes.spines.top": False, "axes.spines.right": False,
    "mathtext.default": "regular", "svg.fonttype": "none",
})

FIG_W = 6.5  # figures/README.md: maximum width 6.5 in, 300 dpi PNG, PNG and SVG only

# Stand-alone panel size in inches. box = whole panel; axes = left, bottom, width, height.
LAYOUT = {
    "A": {"box": (3.32, 1.70), "axes": (0.42, 0.36, 2.85, 1.12)},
    "B": {"box": (3.32, 1.70), "axes": (0.42, 0.36, 2.85, 1.12)},
    "C": {"box": (3.25, 2.04), "axes": (0.06, 0.36, 3.12, 1.46)},
    "D": {"box": (3.25, 2.04), "axes": (0.06, 0.36, 3.12, 1.46)},
    "E": {"box": (3.00, 2.23), "axes": (0.52, 0.36, 2.10, 1.55)},
    "F": {"box": (3.00, 1.99), "axes": (0.52, 0.36, 2.10, 1.43)},
    "G": {"box": (1.76, 1.86), "axes": (0.42, 0.36, 1.30, 1.30)},
    "H": {"box": (1.76, 1.86), "axes": (0.42, 0.36, 1.30, 1.30)},
    "I": {"box": (1.92, 2.45), "axes": (0.55, 0.49, 1.30, 1.30)},
    "J": {"box": (1.92, 2.45), "axes": (0.55, 0.49, 1.30, 1.30)},
}

# Assembled figure, in inches. FIG_H is the total height. Each panel is
# (letter x, letter top measured down from the figure top, axes left,
#  axes top measured down from the figure top, axes width, axes height).
# E sits over F with one shared x-axis title; I sits under G (11 h); J under H (15 h).
FIG_H = 7.75
ROW_C = 1.70   # top of the C, D row
LOWER = 3.74   # top of the E-J block
PLACEMENT = {
    "A": (0.00, 0.00, 0.42, 0.22, 2.85, 1.12),
    "B": (3.30, 0.00, 3.57, 0.22, 2.85, 1.12),
    "C": (0.00, ROW_C, 0.06, ROW_C + 0.22, 3.12, 1.46),
    "D": (3.27, ROW_C, 3.33, ROW_C + 0.22, 3.12, 1.46),
    "E": (0.00, LOWER, 0.52, LOWER + 0.32, 2.10, 1.55),
    "F": (0.00, LOWER + 2.07, 0.52, LOWER + 2.22, 2.10, 1.43),
    "G": (3.05, LOWER, 3.60, LOWER + 0.20, 1.30, 1.30),
    "H": (4.92, LOWER, 5.14, LOWER + 0.20, 1.30, 1.30),
    "I": (3.05, LOWER + 2.03, 3.60, LOWER + 2.22, 1.30, 1.30),
    "J": (4.92, LOWER + 2.03, 5.14, LOWER + 2.22, 1.30, 1.30),
}
NO_Y_TITLE = {"B", "H", "J"}   # share the title of the panel on their left
NO_X_TITLE = {"E"}             # shares the title of F below it
KEY_I_J = (6.47, LOWER + 1.86)  # top-right corner of the I, J key, down from the figure top


def add_axes(fig, panel, x0=0.0, y0=0.0):
    """Add the panel's axes with the panel box's lower-left corner at (x0, y0) in."""
    W, H = fig.get_size_inches()
    left, bottom, width, height = LAYOUT[panel]["axes"]
    return fig.add_axes(((x0 + left) / W, (y0 + bottom) / H, width / W, height / H))


def add_letter(fig, panel, x0=0.0, y0=0.0):
    W, H = fig.get_size_inches()
    _, box_h = LAYOUT[panel]["box"]
    fig.text((x0 + 0.02) / W, (y0 + box_h - 0.02) / H, panel, fontsize=FS + 2,
             fontweight="bold", va="top", ha="left", color=INK)


def write_svg_with_arial(fig, path):
    """Name Arial in the SVG; Liberation Sans has the same letter widths."""
    fig.savefig(path, dpi=600, facecolor=SURFACE)
    text = path.read_text()
    path.write_text(re.sub(r"Liberation Sans", "Arial", text))


def save(fig, stem, folder=PANELS):
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f"{stem}.png", dpi=300, facecolor=SURFACE)
    write_svg_with_arial(fig, folder / f"{stem}.svg")
    plt.close(fig)
    print("wrote", folder / f"{stem}.png and .svg")


def panel_figure(panel):
    """Stand-alone figure for one panel."""
    w, h = LAYOUT[panel]["box"]
    fig = plt.figure(figsize=(w, h), facecolor=SURFACE)
    ax = add_axes(fig, panel)
    add_letter(fig, panel)
    return fig, ax
