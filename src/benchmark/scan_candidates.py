#!/usr/bin/env python3
"""Phase 1 step A: enumerate eligible design targets for the expanded panel.

For every target that has BOTH a cluster_subsets escape library and precomputed
complex_features, check the cheap CPU-only eligibility conditions and emit one TSV row:

  target  n_variants  n_positions  omega_span  n_ab_chains  cdr_len  cdr_seq  vh_seq  ag_seq

Eligibility (hard):
  - complex_features .pt loads and has S / chain_encoding_all
  - AACDB pdb exists
  - cdr_h3_positions() resolves a Cys->WGxG anchored CDR-H3 (this is the filter that
    killed the naive 95-106 window on many targets)
  - >= MIN_VARIANTS escape variants with >= MIN_POSITIONS distinct mutated antigen positions
    (needed for a position-disjoint steer/eval split)

Sequences are emitted so the next step can mmseqs-cluster antigens (70%) and VH (90%).
"""
from __future__ import annotations
import csv, os, sys, traceback
from multiprocessing import Pool

import torch

sys.path.insert(0, f"{VAAD_ROOT}/projects/mutation_sampling/lyra/scripts/nos_diffusion")
sys.path.insert(0, f"{VAAD_ROOT}/projects/mutation_sampling/lyra/src")
from cluster_steer import cdr_h3_positions, parse_chains, to_seq  # noqa: E402
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

SUBSETS = f"{VAAD_ROOT}/results/cluster_subsets"
FEATS = f"{VAAD_ROOT}/scratch/complex_features"
PDBS = f"{VAAD_ROOT}/datasets/AACDB/data_zip/complex_structure"
MIN_VARIANTS = 24
MIN_POSITIONS = 4          # matches the loosest target in the existing 51-panel (3OGO_GB, 7ZRA_FC)

csv.field_size_limit(10_000_000)


def scan(target):
    try:
        pt = os.path.join(FEATS, f"{target}.pt")
        pdb = os.path.join(PDBS, f"{target}.pdb")
        if not os.path.exists(pt):
            return ("skip", target, "no_features")
        if not os.path.exists(pdb):
            return ("skip", target, "no_pdb")

        # --- escape library stats (streamed; some subsets are 100k+ rows) ---
        n_var, positions, omegas = 0, set(), []
        with open(os.path.join(SUBSETS, f"{target}.csv"), newline="") as fh:
            r = csv.reader(fh)
            hdr = next(r)
            pi = hdr.index("mutated_positions")
            oi = hdr.index("omega_normalized")
            ai = hdr.index("accepted_or_rejected") if "accepted_or_rejected" in hdr else None
            for row in r:
                if ai is not None and len(row) > ai and row[ai] != "accepted":
                    continue
                n_var += 1
                if len(row) > pi:
                    positions.update(row[pi].split("|"))
                if len(row) > oi:
                    try:
                        omegas.append(float(row[oi]))
                    except ValueError:
                        pass
        if n_var < MIN_VARIANTS:
            return ("skip", target, f"few_variants:{n_var}")
        if len(positions) < MIN_POSITIONS:
            return ("skip", target, f"few_positions:{len(positions)}")

        # --- structure / CDR-H3 registration ---
        base = torch.load(pt, map_location="cpu", weights_only=False)
        if isinstance(base, dict) and "S" not in base and len(base) == 1:
            base = list(base.values())[0]
        n_ab = len(parse_chains(target))
        cdr = cdr_h3_positions(base, n_ab)
        S = base["S"][0].tolist()
        ce = base["chain_encoding_all"][0].tolist()
        seq = to_seq(S)
        cdr_seq = "".join(seq[i] for i in cdr)

        ab_vals = sorted(set(ce))[:n_ab]
        cdr_chain = ce[cdr[0]]
        vh = "".join(seq[i] for i in range(len(S)) if ce[i] == cdr_chain)
        ag = "".join(seq[i] for i in range(len(S)) if ce[i] not in ab_vals)
        if len(ag) < 30:
            return ("skip", target, f"tiny_antigen:{len(ag)}")

        span = (max(omegas) - min(omegas)) if omegas else 0.0
        return ("ok", target, (n_var, len(positions), f"{span:.6g}", n_ab,
                               len(cdr), cdr_seq, vh[:400], ag[:1200]))
    except Exception as e:
        return ("skip", target, f"err:{type(e).__name__}:{e}"[:120])


def main():
    targets = sorted(f[:-4] for f in os.listdir(SUBSETS) if f.endswith(".csv"))
    print(f"[scan] {len(targets)} targets with escape libraries", flush=True)
    out = f"{VAAD_ROOT}/scratch/genzoo/cand/candidates.tsv"
    rej = f"{VAAD_ROOT}/scratch/genzoo/cand/rejected.tsv"
    nok = 0
    with Pool(48) as p, open(out, "w") as fo, open(rej, "w") as fr:
        fo.write("target\tn_variants\tn_positions\tomega_span\tn_ab_chains\tcdr_len\tcdr_seq\tvh_seq\tag_seq\n")
        for i, (status, target, payload) in enumerate(p.imap_unordered(scan, targets, chunksize=4), 1):
            if status == "ok":
                nok += 1
                fo.write(target + "\t" + "\t".join(str(x) for x in payload) + "\n")
            else:
                fr.write(f"{target}\t{payload}\n")
            if i % 100 == 0:
                print(f"  {i}/{len(targets)} scanned, {nok} eligible", flush=True)
    print(f"[scan] DONE eligible={nok}/{len(targets)} -> {out}", flush=True)


if __name__ == "__main__":
    main()
