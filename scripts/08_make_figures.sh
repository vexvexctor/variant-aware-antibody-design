#!/usr/bin/env bash
# Regenerate every figure from released result CSVs (no hard-coded numbers).
set -euo pipefail
cd "$(dirname "$0")/.."
python3 src/analysis/plots/plot_evaluator_agreement.py
python3 src/analysis/plots/plot_generator_by_family.py
python3 src/analysis/plots/plot_panel_size.py
ls -1 figures/
