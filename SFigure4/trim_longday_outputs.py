#!/usr/bin/env python3
"""Reduce the long-day pipeline outputs to small tables kept beside Fig. S4.

Usage: trim_longday_outputs.py WORK_GENERATABLE_DIR FIGURE_DIR
Writes longday_response_calls.csv (genes called in any residual-response
category, with their key statistics) and the significant rows of the GO and
regulon enrichment tests for those categories.
"""
import sys
from pathlib import Path
import pandas as pd

gen, out = Path(sys.argv[1]), Path(sys.argv[2])
c = pd.read_csv(gen / "diurnal_adjusted_response_categories.csv.gz", low_memory=False)
called = (c.called_category != "not_called") | c.collection_phase_specific_category.isin(
    ["light_cold_specific_induction", "dark_cold_specific_induction"])
cols = ["gene_id", "hours_after_onset", "called_category", "category_q_value",
        "collection_phase_specific_category", "collection_phase_category_q_value",
        "mean_adjusted_dawn_response_log2", "mean_adjusted_dusk_response_log2", "mean_external_tpm"]
c.loc[called, cols].to_csv(out / "longday_response_calls.csv", index=False)
for name, f in [("go", "go_enrichment_diurnal_adjusted_responses.tsv.gz"),
                ("regulon", "regulon_enrichment_diurnal_adjusted_responses.tsv.gz")]:
    e = pd.read_csv(gen / f, sep="\t")
    e[e.significant_enrichment].to_csv(out / f"longday_{name}_enrichment_significant.csv", index=False)
print("Trimmed long-day outputs written to", out)
