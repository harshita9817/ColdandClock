# Figure 4: raffinose pathway induction and raffinose accumulation

- **A:** diurnal-adjusted expression change 27 hours after dawn- and
  dusk-onset chilling for myo-inositol-3-phosphate synthase
  (Sobic.001G472800), inositol monophosphatase (Sobic.009G048900),
  galactinol synthase (*SbGolS2*, Sobic.002G423600), and raffinose synthase
  (Sobic.003G052300), read from
  `data/generatable/diurnal_adjusted_response_categories.csv.gz`. Arrow color
  shows the response category (dusk-specific, shared, or not classified).
  Plotted values are written to `Figure4_pathway_expression.csv`.
- **B:** raffinose change 27 hours after each onset relative to time-matched
  controls, estimated from five growth trials.

Run the model first, then the plot:

```sh
Rscript figures/Figure4/fit_figure4_raffinose.R
python3 figures/Figure4/plot_figure4.py
```

## Files

- `Figure4_raffinose_samples.csv`: one row per GC-MS injection (270 rows).
  `log2_raffinose` is the log2 raffinose peak height normalized to sample
  weight; it is missing for three failed injections. `T` is hours from dawn
  onset (27 hours after dusk onset is `T39`), `date` is the growth trial,
  `run` the GC-MS run, `order` the scaled injection order, `rib` the scaled
  ribitol internal standard, and `cell` the treatment-by-time group. Peak
  integration and weight normalization were done upstream of this repository.
- `fit_figure4_raffinose.R`: robust (bisquare M-estimator) regression of
  `log2_raffinose` on treatment-by-time cell, growth trial, GC-MS run,
  injection order, and ribitol. Writes `Figure4_raffinose_model.csv` with the
  chilled-minus-control contrasts for each onset and their dusk-minus-dawn
  difference (log2 change, 95% CI, two- and one-sided p, percent change).
- `plot_figure4.py`: renders `Figure4.png` (300 dpi) and `Figure4.svg`.
  Panel B bar heights and significance stars are read from the model output
  (chilled vs control two-sided; dusk vs dawn one-sided).
