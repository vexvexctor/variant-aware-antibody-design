#!/usr/bin/env bash
# Controls and ablations (A-F). Analyses run from released CSVs; generation needs a cluster.
set -euo pipefail
cd "$(dirname "$0")/.."
echo "A temperature-matched   -> src/analysis/temperature_control.py   results/controls/TEMPERATURE_CONTROL.md"
echo "B edit-distance null    -> src/analysis/edit_matched_null.py     results/controls/lift_over_random_null.md"
echo "C objective ablation    -> configs/steering.yaml (steer-agg)     results/controls/EXP4_CVAR_MATRIX.md"
echo "D panel-size M=1/3/6/12 -> src/analysis/panel_size.py            results/controls/panel_size_summary.csv"
echo "E budget-matched search -> src/analysis/budget_matched_search.py results/controls/EXP1_BUDGET.md"
echo "F RFantibody (de-novo)  -> src/generators/rfantibody_*.py        graded on its OWN backbone"
python3 src/analysis/paired_test.py 2>/dev/null || true
