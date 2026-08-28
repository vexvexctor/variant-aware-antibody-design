#!/usr/bin/env bash
# Rebuild RobustBench manifests from AACDB (expensive; released manifests are in data/).
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VAAD_ROOT:?}"
echo "[1/4] extract antigen sequences with the historical rule"
python3 src/benchmark/extract_antigens.py all data/processed/AACDB_all_antigens.fasta
echo "[2/4] cluster full AACDB + benchmark subset (mmseqs2)"
bash src/benchmark/run_aacdb_mmseqs.sh data/processed/AACDB_all_antigens.fasta data/sequences/antigens.fasta
echo "[3/4] scan structurally eligible candidates"
python3 src/benchmark/scan_candidates.py
echo "[4/4] select the 200-target panel and build 6/24 variant panels"
python3 src/benchmark/select_targets.py
python3 src/benchmark/build_panels.py --help >/dev/null
echo "done -> data/manifests/, data/splits/"
