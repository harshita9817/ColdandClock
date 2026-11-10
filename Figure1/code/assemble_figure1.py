#!/usr/bin/env python3
"""Code-only layout (earlier scatter I and J), saved to earlier_versions/.

The current Figure1.svg holds hand edits, so this does not write it. insert_I_J_panels.py takes the
updated G and H from this layout, because their points are embedded images.
"""
import matplotlib.pyplot as plt

import panel_A_B_freezing as ab
import panel_C_D_alluvial as cd
import panel_E_deg_counts as e
import panel_F_overlap as f
import panel_G_H_gene_scatter as gh
import panel_I_J_go_scatter as ij
from figure1_style import (FIG_H, FIG_W, FS, INK, KEY_I_J, NO_X_TITLE, NO_Y_TITLE, PLACEMENT, ROOT, SURFACE,
                           save)

DRAW = {"A": ab.draw, "B": ab.draw, "C": cd.draw, "D": cd.draw, "E": e.draw, "F": f.draw,
        "G": gh.draw, "H": gh.draw, "I": ij.draw, "J": ij.draw}


def main(stem="Figure1_code_layout_scatter_IJ", folder=ROOT / "earlier_versions"):
    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor=SURFACE)
    for panel, (letter_x, letter_top, left, top, width, height) in PLACEMENT.items():
        ax = fig.add_axes((left / FIG_W, (FIG_H - top - height) / FIG_H, width / FIG_W, height / FIG_H))
        if panel in NO_X_TITLE:
            DRAW[panel](ax, panel, panel not in NO_Y_TITLE, xlabel=False)
        else:
            DRAW[panel](ax, panel, panel not in NO_Y_TITLE)
        fig.text((letter_x + 0.02) / FIG_W, (FIG_H - letter_top - 0.01) / FIG_H, panel, fontsize=FS + 2,
                 fontweight="bold", va="top", ha="left", color=INK)
    ij.add_key(fig, KEY_I_J[0], FIG_H - KEY_I_J[1], anchor="upper right")  # one key shared by I and J
    print(f"Figure 1: {FIG_W} x {FIG_H} in ({FIG_W * 25.4:.0f} x {FIG_H * 25.4:.0f} mm)")
    save(fig, stem, folder)


if __name__ == "__main__":
    main()
