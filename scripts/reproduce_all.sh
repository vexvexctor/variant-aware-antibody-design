#!/usr/bin/env bash
# Reproduce the paper. Default rebuilds statistics+figures from released results.
set -euo pipefail
cd "$(dirname "$0")/.."
MODE=${1:---from-precomputed}
case "$MODE" in
  --from-precomputed)
    echo "== rebuilding paper statistics and figures from released result files =="
    python3 -m unittest discover -s tests -q
    bash scripts/07_run_secondary_analyses.sh
    bash scripts/08_make_figures.sh
    echo "== done: results/ and figures/ regenerated =="
    ;;
  --full)
    echo "== FULL pipeline: multi-GPU, days of compute, staged dependencies required =="
    : "${VAAD_ROOT:?set VAAD_ROOT}"
    bash scripts/00_download_dependencies.sh
    bash scripts/01_build_robustbench.sh
    bash scripts/02_generate_baselines.sh
    echo "03: run src/steering/variant_aware.py per target (array job) -- see REPRODUCIBILITY.md"
    echo "04: then scripts/04_evaluate_heldout.sh"
    bash scripts/05_run_controls.sh
    bash scripts/06_run_generalization.sh
    bash scripts/07_run_secondary_analyses.sh
    bash scripts/08_make_figures.sh
    ;;
  *) echo "usage: reproduce_all.sh [--from-precomputed|--full]"; exit 2 ;;
esac
