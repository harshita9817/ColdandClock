# Fig. S4: analyses repeated with a long-day diurnal reference

The Figure 3 diurnal-progression analyses and the residual-response
classification repeated with an independent BTx623 leaf diel series grown under
14 h light/10 h dark (JGI Project 1051031), in place of the 12 h/12 h series of
Lai et al. (2020).

- **A--D:** as in Figure 3A--D.
- **E:** genes assigned to each residual-response category: shared, dawn-specific,
  and dusk-specific induction at 27 h, and light-cold- and dark-cold-specific
  induction from the 3- and 15-h samples.

## Run

```sh
bash figures/SFigure4_longday_reference/run_SFigure4_longday_reference.sh
```

The driver copies the repository to a temporary directory, replaces the
12L:12D reference input with `longday_diurnal_tpms.csv.gz`, and re-runs the
unmodified pipeline there: `build_external_diurnal_reference.R`,
`generate_diurnal_adjusted_response_table.R`, the residual-response GO and
regulon enrichment scripts, and `plot_SFigure4_longday_reference.R`. The
repository's own `data/generatable/` files are not modified. It needs the same
environment as the main analysis (R with MetaCycle 1.2.1 and limma 3.68.4; the
Python packages in `requirements.txt`), plus the locally generated
`data/generatable/deseq2_timepoint_contrasts.csv.gz` and the GO ontology file
in `data/externalsourcerequired/`.

The committed outputs were generated with R 4.5.3 and limma 3.68.5 (the pinned
3.68.4 was not available for this platform). Before use, this environment was
checked by regenerating the committed 12L:12D reference and response table:
every value agreed to within 1e-11 and every category call was identical.

## Files

- `longday_diurnal_tpms.csv.gz`: TPM for the 24 libraries from the first
  complete 24-h cycle (ZT0--ZT21, three replicates per phase), in the layout
  read by `build_external_diurnal_reference.R` (column `SorghumN` is phase
  ZT = 3N mod 24 in replicate series ceiling(3N/24)).
- `longday_diurnal_samples.csv`: maps each column to its JGI library ID,
  replicate, hour, and ZT.
- `plot_SFigure4_longday_reference.R`: `figures/Figure3/plot_figure3.R` with
  Dawn/Dusk legend labels, panel E, and S4 output names.
- `trim_longday_outputs.py`: writes the small tables below from the full
  long-day outputs.
- `SFigure4_*.csv`: plotted values (as for Figure 3) and panel E counts.
- `longday_response_calls.csv`: genes called in any residual-response category.
- `longday_go_enrichment_significant.csv`,
  `longday_regulon_enrichment_significant.csv`: significant enrichment results
  for the long-day categories (15,835-gene universe).
