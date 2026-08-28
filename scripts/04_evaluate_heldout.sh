#!/usr/bin/env bash
# Held-out evaluation on the 24 variants never seen during design.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VAAD_ROOT:?}"
MANIFEST=${1:?usage: 04_evaluate_heldout.sh <manifest.tsv> <valdir> <outdir>}
VALDIR=${2:?}; OUTDIR=${3:?}
mkdir -p "$OUTDIR"
# H3-DDG is the primary held-out evaluator (BA-DDG is the steering objective; never grade on it)
python3 src/evaluators/h3_ddg.py --manifest "$MANIFEST" --valdir "$VALDIR" --outdir "$OUTDIR"
