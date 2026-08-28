#!/usr/bin/env bash
# Stage third-party models/data. Nothing here is redistributed by this repository.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VAAD_ROOT:?set VAAD_ROOT to a writable data root}"
: "${VAAD_TOOLS:=$VAAD_ROOT/tools}"
mkdir -p "$VAAD_ROOT"/{datasets,results,work} "$VAAD_TOOLS"
cat <<'TXT'
Fetch these yourself; licences differ and none permit blanket redistribution:

  AACDB structures   -> $VAAD_ROOT/datasets/AACDB/complex_structure/*.pdb
  ProteinMPNN        -> $VAAD_TOOLS/ProteinMPNN            (MIT)          github.com/dauparas/ProteinMPNN
  AntiFold           -> $VAAD_TOOLS/antifold                              github.com/oxpig/AntiFold
  AbMPNN checkpoint  -> $VAAD_TOOLS/abmpnn/abmpnn.pt        (check terms)
  ESM-IF1 / ESM-1v   -> via fair-esm (torch hub cache)      (MIT)
  Pythia-PPI         -> $VAAD_TOOLS/pythia
  RFantibody         -> $VAAD_TOOLS/RFantibody
  FoldX 5            -> $VAAD_TOOLS/foldx  PROPRIETARY: obtain your own academic licence.
                        Put molecules/ beside the binary. Not redistributable.
  MM-GBSA            -> OpenMM + pdbfixer in a venv; export OPENMM_PLUGIN_DIR=<env>/lib/plugins

mmseqs2 must be on PATH for benchmark clustering.
TXT
