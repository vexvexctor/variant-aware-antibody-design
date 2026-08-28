#!/usr/bin/env bash
# Variant-aware guided design (GPU, cluster-scale).
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VAAD_ROOT:?}"
TARGET=${1:?usage: 03_run_variant_aware_design.sh <TARGET_ID> [OUT.json]}
OUT=${2:-data/processed/${TARGET}_design.json}
python3 src/steering/variant_aware.py \
  --target "$TARGET" \
  --clustered-csv "$VAAD_ROOT/results/cluster_subsets/$TARGET.csv" \
  --aacdb-dir "$VAAD_ROOT/datasets/AACDB/complex_structure" \
  --features-dir "$VAAD_ROOT/work/complex_features" \
  --proteinmpnn-dir "$VAAD_TOOLS/ProteinMPNN" \
  --proteinmpnn-weights "$VAAD_TOOLS/ProteinMPNN/vanilla_model_weights/v_48_002.pt" \
  --auto-cdr-h3 \
  --mechanism svdd --base antifold --oracle baddg_ddg \
  --steer-agg cvar --cvar-alpha 0.2 \
  --max-cluster 6 --n-designs 6 --svdd-k 6 \
  --holdout-by omega --eval-cluster 24 --holdout-seed 0 \
  --out "$OUT"
