#!/usr/bin/env bash
# Unguided baselines: each frozen generator alone (GPU, cluster-scale).
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VAAD_ROOT:?}"
TARGETS=${1:-data/manifests/robustbench_targets.csv}
echo "Runs the frozen generator with NO reward. See configs/generators.yaml for per-generator"
echo "settings and scripts/03_run_variant_aware_design.sh for the guided counterpart."
echo "Baseline temperature is 0.2 and steered sampling is 1.0 -- the temperature-matched"
echo "control in 05_run_controls.sh runs both at 1.0. Quote both framings."
python3 src/steering/variant_aware.py --help >/dev/null
echo "(entry point: src/steering/variant_aware.py with --mechanism svdd --svdd-k 1)"
