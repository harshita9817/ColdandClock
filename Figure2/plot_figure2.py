#!/usr/bin/env python3
"""Figure 2: frozen diurnal progression creates apparent cold responses.

A: conceptual model. Controls continue to cycle while chilled plants hold the
   expression level present at chilling onset.
B: expression of the sorghum SFR2 ortholog (Sobic.005G217000) in controls and
   after dawn- and dusk-onset chilling. Hours are measured from dawn onset.

The per-sample TPM values plotted in B are written to
Figure2_SFR2_expression.csv.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm

FIGDIR = Path(__file__).resolve().parent
ROOT = FIGDIR.parent.parent
GENE = "Sobic.005G217000"

# Nimbus Sans (URW Helvetica clone) if installed; otherwise matplotlib falls back to its default sans.
for f in ["NimbusSans-Regular.otf", "NimbusSans-Bold.otf"]:
    p = Path("/usr/share/fonts/opentype/urw-base35") / f
    if p.exists():
        fm.fontManager.addfont(str(p))
installed = {f.name for f in fm.fontManager.ttflist}
FONT = next((f for f in ["Nimbus Sans", "Helvetica", "Arial"] if f in installed), "DejaVu Sans")
FS = 9
plt.rcParams.update({
    "font.family": FONT,
    "font.size": FS, "axes.labelsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS,
    "legend.fontsize": FS, "text.color": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
})
COL = {"CONT": "#4D4D4D", "DAWN": "#C77700", "DUSK": "#7B3294"}
LAB = {"CONT": "Control", "DAWN": "Dawn", "DUSK": "Dusk"}
LW = 2.0

# ---- Extract the plotted gene from the full TPM matrix -------------------------
meta = pd.read_csv(ROOT / "data/provided/sorghum_cold_timecourse_sample_metadata.csv")
tpm = pd.read_csv(ROOT / "data/externalsourcerequired/kallisto_tpms.csv.gz", index_col=0)
row = tpm.loc[[i for i in tpm.index if i.split(".")[0] + "." + i.split(".")[1] == GENE]]
if len(row) != 1:
    raise SystemExit(f"Expected one row for {GENE}, found {len(row)}")
row = row.iloc[0]
# Time = hours from dawn onset. Controls span a full cycle; chilled samples start 3 h after their onset.
times = {"CONT": [-1] + list(range(3, 40, 4)), "DAWN": list(range(3, 28, 4)), "DUSK": list(range(15, 40, 4))}
parts = []
for k, ts in times.items():
    d = meta[(meta.Treatment == k) & meta.Time.isin(ts)][["Sample name", "Treatment", "Time"]].copy()
    d["tpm"] = row[d["Sample name"]].values
    parts.append(d)
data = pd.concat(parts, ignore_index=True)
data.insert(0, "gene_id", GENE)
data.to_csv(FIGDIR / "Figure2_SFR2_expression.csv", index=False)

# ---- Plot ----------------------------------------------------------------------
fig, (a, b) = plt.subplots(1, 2, figsize=(6.5, 2.2), gridspec_kw={"width_ratios": [1, 1], "wspace": 0.32})

t = np.linspace(-4, 22, 400)
a.plot(t, np.cos(2 * np.pi * t / 24), color=COL["CONT"], lw=LW, solid_capstyle="round", label="Control")
a.plot([0, 12], [1, 1], color=COL["DAWN"], lw=LW, solid_capstyle="round", label="Chilled")
a.annotate("", xy=(12, -0.93), xytext=(12, 0.93),
           arrowprops=dict(arrowstyle="<->", color="black", lw=1.0, shrinkA=0, shrinkB=0, mutation_scale=9))
a.set_xticks([]); a.set_yticks([]); a.set_xlim(-4.5, 22.5); a.set_ylim(-1.12, 1.12)
a.set_xlabel("Time"); a.set_ylabel("Expression")
a.legend(frameon=False, loc="lower left", handlelength=1.6, borderaxespad=0.1, labelspacing=0.3)

for s, e in [(13, 23), (37, 40.6)]:  # dark periods (lights off ZT14, on ZT0; dawn onset at ZT1)
    b.axvspan(s, e, color="#EDEDED", zorder=0, lw=0)
for k in ["CONT", "DAWN", "DUSK"]:
    d = data[data.Treatment == k]
    g = d.groupby("Time").tpm.mean()
    b.scatter(d.Time, d.tpm, s=10, color=COL[k], alpha=0.3, lw=0, zorder=2)
    b.plot(g.index, g.values, color=COL[k], lw=LW, solid_capstyle="round", solid_joinstyle="round",
           label=LAB[k], zorder=3)
b.set_ylim(0, 500); b.set_xlim(-2, 40.6)
b.set_xticks([-1, 7, 15, 23, 31, 39]); b.set_yticks([0, 100, 200, 300, 400, 500])
b.set_xlabel("Hours"); b.set_ylabel(f"{GENE} Exp. (TPM)")
b.legend(frameon=False, loc="upper right", handlelength=1.6, borderaxespad=0.1, labelspacing=0.3)

for ax, letter in [(a, "A"), (b, "B")]:
    ax.text(0.03, 0.97, letter, transform=ax.transAxes, fontsize=FS, fontweight="bold", va="top", ha="left")

fig.savefig(FIGDIR / "Figure2.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(FIGDIR / "Figure2.svg", bbox_inches="tight", facecolor="white")
print("Figure 2 written to", FIGDIR)
