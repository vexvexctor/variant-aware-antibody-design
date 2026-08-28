#!/usr/bin/env bash
# Secondary analyses: evaluator agreement, affinity/retention, trajectories, biology, compute.
set -euo pipefail
cd "$(dirname "$0")/.."
cd src && python3 -m analysis.evaluator_agreement \
    --pool ../data/processed/shared_1785_design_pool.csv \
    --outdir ../results/evaluator_agreement && cd ..
echo "affinity vs retention -> src/analysis/affinity_retention.py (+ _ancova.py)"
echo "trajectories          -> results/trajectories/ (DESCRIPTIVE: partial capture, 121 runs)"
echo "residue features      -> src/analysis/residue_features.py, epitope_charge.py"
echo "compute cost          -> src/analysis/compute_cost.py -> results/compute/overhead_by_arm.csv"
