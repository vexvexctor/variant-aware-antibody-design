#!/usr/bin/env python3
"""Build the antigen-escape panel JSON for one or more targets WITHOUT steering.

The native-fragility gate (native_fragility_gate.py) only needs antibody_chains,
antigen_chains, and the antigen_variants panel — not any CDR design (native fragility
is design-independent). This writes exactly that, reusing the tested build_cluster
mutation-stamping. Pure CPU, no diffusion model, no GPU — so it's cheap to run across
many candidate targets for the F4 fragility screen (PAPER_PIPELINE_REPORT.md §12.7).

One pass over the (large) clustered CSV buckets rows by complex_id for all requested
targets, then builds each panel.

Usage
-----
    PYTHONPATH=src python3 scripts/nos_diffusion/build_panel.py \
        --targets 7XYZ_HLA 8ABC_HLB ... \
        --clustered-csv $RESULTS/clustered.csv \
        --features-dir $SCRATCH/complex_features \
        --aacdb-dir .../complex_structure \
        --max-cluster 24 --out-dir $RESULTS/panels
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import torch

from cluster_steer import build_cluster, parse_chains, to_seq, N_AA  # noqa: F401


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", required=True)
    ap.add_argument("--clustered-csv", required=True)
    ap.add_argument("--features-dir", required=True)
    ap.add_argument("--aacdb-dir", required=True)
    ap.add_argument("--max-cluster", type=int, default=24)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    want = set(args.targets)

    # one fast pass: bucket minimal rows (complex_id, mutation_string, omega) by target
    buckets = {t: [] for t in want}
    with open(args.clustered_csv, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        ci = header.index("complex_id")
        mi = header.index("mutation_string")
        oi = header.index("omega_normalized")
        for n, row in enumerate(reader, 1):
            if len(row) > ci and row[ci] in want:
                buckets[row[ci]].append(
                    {"complex_id": row[ci], "mutation_string": row[mi], "omega_normalized": row[oi]})
            if n % 5_000_000 == 0:
                print(f"  scanned {n:,} rows...", file=sys.stderr, flush=True)

    ok, bad = [], []
    for t in args.targets:
        rows = buckets[t]
        pdb_path = f"{args.aacdb_dir}/{t}.pdb"
        feats_path = f"{args.features_dir}/{t}.pt"
        if not rows or not os.path.exists(pdb_path) or not os.path.exists(feats_path):
            print(f"[{t}] SKIP (rows={len(rows)}, pdb={os.path.exists(pdb_path)}, "
                  f"feats={os.path.exists(feats_path)})", flush=True)
            bad.append(t)
            continue
        base = torch.load(feats_path, map_location="cpu", weights_only=False)
        antibody = parse_chains(t)
        chains_str = t.split("_")[-1]
        antigen = [c for c in chains_str if c not in antibody]
        cluster, ag_len, aligned, (n_parse, n_map, n_wt), labels = build_cluster(
            base, t, rows, args.max_cluster, pdb_path, "cpu")
        if not cluster:
            print(f"[{t}] SKIP (empty cluster; aligned={aligned}, "
                  f"skipped parse {n_parse}/resnum {n_map}/wt {n_wt})", flush=True)
            bad.append(t)
            continue
        antigen_variants = [{"name": "WT", "mutations": []}]
        for (name, toks), (_, omega) in zip(labels, cluster):
            antigen_variants.append({
                "name": name, "omega": round(omega, 5),
                "mutations": [{"chain": c, "wt": wt, "resnum": rn, "icode": ic, "mut": mut}
                              for (c, wt, rn, ic, mut) in toks]})
        out = {"target": t, "pdb": os.path.basename(pdb_path),
               "antibody_chains": antibody, "antigen_chains": antigen,
               "antigen_variants": antigen_variants, "alignment_ok": bool(aligned)}
        out_path = f"{args.out_dir}/{t}_panel.json"
        with open(out_path, "w") as fhj:
            json.dump(out, fhj, indent=2)
        print(f"[{t}] panel: {len(antigen_variants)-1} mutants, aligned={aligned} -> {out_path}", flush=True)
        ok.append(t)

    print(f"\nbuilt {len(ok)} panels, skipped {len(bad)}"
          + (f" ({', '.join(bad)})" if bad else ""))


if __name__ == "__main__":
    main()
