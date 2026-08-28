#!/usr/bin/env bash
# Generalization: out-of-family evaluator, post-cutoff antigens, natural variants.
set -euo pipefail
cd "$(dirname "$0")/.."
echo "Pythia-PPI (post-selection only, paired deltas):"
echo "  src/evaluators/pythia_ppi.py -> results/generalization/PYTHIA_OUT_OF_FAMILY.txt"
echo "Post-cutoff / low-similarity antigens: data/manifests/postcutoff_targets.csv"
echo "Natural SARS-CoV-2 RBD variants:       data/manifests/natural_variants.csv"
python3 src/analysis/pythia_beatnative.py --help >/dev/null 2>&1 || true
