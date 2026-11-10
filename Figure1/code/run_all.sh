#!/usr/bin/env bash
# Rebuild all panel data and every panel, then put I and J into the hand-edited figure.
# Panels A-H in Figure1.svg are the hand-edited copies: they do not update from this script.
set -euo pipefail
cd "$(dirname "$0")"
python3 prepare_panel_data.py
python3 panel_A_B_freezing.py
python3 panel_C_D_alluvial.py
python3 panel_E_deg_counts.py
python3 panel_F_overlap.py
python3 panel_G_H_gene_scatter.py
python3 panel_I_J_bars.py
python3 assemble_figure1.py
python3 insert_I_J_panels.py
