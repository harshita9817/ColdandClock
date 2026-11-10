#!/usr/bin/env bash
# Fig. S4: repeat the Figure 3 and residual-response analyses with the long-day
# (14L:10D) JGI BTx623 diel series as the external diurnal reference.
#
# The repository is copied to a temporary directory, the 12L:12D reference input
# (data/externalsourcerequired/sorghum_2020_diurnal_tpms.csv.gz) is replaced by
# longday_diurnal_tpms.csv.gz from this directory, and the unmodified pipeline
# scripts are re-run there. Outputs are copied back into this directory; the
# repository's own data/generatable files are never touched.
set -euo pipefail
FIGDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$FIGDIR/../.." && pwd)"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

rsync -a --exclude .git --exclude figures "$ROOT/" "$WORK/"
mkdir -p "$WORK/figures/SFigure4_longday_reference"
cp "$FIGDIR/longday_diurnal_tpms.csv.gz" "$WORK/data/externalsourcerequired/sorghum_2020_diurnal_tpms.csv.gz"
cp "$FIGDIR/plot_SFigure4_longday_reference.R" "$WORK/figures/SFigure4_longday_reference/"

Rscript "$WORK/scripts/build_external_diurnal_reference.R"
Rscript "$WORK/scripts/generate_diurnal_adjusted_response_table.R"
python3 "$WORK/scripts/prepare_go_tested_sets.py" --mode residual-response
python3 "$WORK/scripts/go_enrichment_timepoint_contrasts.py" --mode residual-response
python3 "$WORK/scripts/prepare_regulon_tested_sets.py" --mode residual-response
python3 "$WORK/scripts/regulon_enrichment_timepoint_contrasts.py" --mode residual-response
Rscript "$WORK/figures/SFigure4_longday_reference/plot_SFigure4_longday_reference.R"

cp "$WORK"/figures/SFigure4_longday_reference/SFigure4_* "$FIGDIR/"
python3 "$FIGDIR/trim_longday_outputs.py" "$WORK/data/generatable" "$FIGDIR"
echo "Fig. S4 outputs copied to $FIGDIR"
