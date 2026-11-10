# Figure 2: frozen diurnal progression creates apparent cold responses

`plot_figure2.py` renders the two-panel figure:

- **A:** conceptual model. Controls (gray) continue to cycle while chilled
  plants (orange) hold the expression level present at chilling onset; the
  arrow marks the resulting apparent difference 12 hours later.
- **B:** expression of the sorghum *SFR2* ortholog (Sobic.005G217000) in
  controls and after dawn- and dusk-onset chilling. Hours are measured from
  dawn onset; lines connect means, points are individual samples, and gray
  shading marks the dark period.

The script reads the sample metadata and the cold-experiment TPM matrix, writes
the plotted per-sample values to `Figure2_SFR2_expression.csv`, and saves
`Figure2.png` (300 dpi) and `Figure2.svg`. Text uses Nimbus Sans if it is
installed (URW base35 fonts) and otherwise falls back to a default sans font.

```sh
python3 figures/Figure2/plot_figure2.py
```
