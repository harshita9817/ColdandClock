#!/usr/bin/env bash
# Rebuild Figure 2.
#
#   bash code/run_all.sh          redraw the figure from the committed data (about 15 s)
#   bash code/run_all.sh --all    first rerun every upstream step, then redraw
#
# --all overwrites the files in data/generatable/. The upstream steps check their
# pinned package versions (MetaCycle 1.2.1, limma 3.68.4) and stop if they differ.
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ "${1:-}" == "--all" ]]; then
  Rscript "$PWD/code/deseq2_timepoint_contrasts.R"              # -> deseq2_timepoint_contrasts.csv.gz
  Rscript "$PWD/code/build_external_diurnal_reference.R"        # -> external_diurnal_reference.tsv.gz
  Rscript "$PWD/code/generate_diurnal_adjusted_response_table.R" # -> diurnal_adjusted_response_categories.csv.gz
  python3 code/prepare_go_tested_sets.py                        # -> go_tested_sets.tsv.gz
fi
Rscript "$PWD/code/plot_figure2.R"
