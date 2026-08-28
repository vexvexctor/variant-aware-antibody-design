#!/usr/bin/env python3
"""Collect RFantibody designs into the shared design-JSON schema.

Reads the ProteinMPNN-sequenced RFdiffusion outputs, pulls out the CDR-H3 residues (RFantibody
carries them as `REMARK PDBinfo-LABEL: <resnum> H3`), and writes the same design JSON every
other arm produces so the graders and aggregator run unmodified.

The H3 length was pinned to the native length at diffusion time, so the designed loop threads
onto the native CDR positions 1:1 -- without that the native-backbone grading protocol would
not be defined at all.

Two grading protocols are supported, per the campaign decision:
  * native-backbone (this JSON) -- comparable to every other arm, but does not credit the
    conformation RFantibody designed;
  * own-structure -- the designed complexes are left in place and pointed at by `designed_pdbs`
    in the JSON, for FoldX/co-fold scoring on their own backbones.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import time

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def h3_sequence(pdb_path):
    """CDR-H3 one-letter sequence from an RFantibody HLT-format output."""
    labels = set()
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("REMARK PDBinfo-LABEL:") and line.rstrip().endswith("H3"):
                labels.add(int(line.split(":")[1].split()[0]))
    if not labels:
        return None
    seen, seq = set(), []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM") or line[21] != "H":
                continue
            rn = int(line[22:26])
            key = (rn, line[26].strip())
            if rn in labels and key not in seen:
                seen.add(key)
                seq.append(AA3.get(line[17:20].strip(), "X"))
    return "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--seq-dir", required=True, help="ProteinMPNN output dir (HLT pdbs)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-designs", type=int, default=6)
    ap.add_argument("--started", type=float, default=None, help="chain start time for cost accounting")
    args = ap.parse_args()

    tmpl = json.load(open(args.template))
    target, native_cdr = tmpl["target"], tmpl["native_cdr"]

    seqs, kept = [], []
    for f in sorted(glob.glob(os.path.join(args.seq_dir, "*.pdb"))):
        s = h3_sequence(f)
        if not s:
            continue
        if len(s) != len(native_cdr):
            print(f"[{target}] skip {os.path.basename(f)}: H3 len {len(s)} != native "
                  f"{len(native_cdr)} (cannot thread onto the native backbone)", flush=True)
            continue
        if any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in s):
            continue
        seqs.append(s)
        kept.append(os.path.abspath(f))
        if len(seqs) >= args.n_designs:
            break
    if not seqs:
        raise SystemExit(f"[{target}] RFantibody produced no length-matched CDR-H3")
    while len(seqs) < args.n_designs:
        seqs.append(seqs[-1])
        kept.append(kept[-1])

    names = [f"cluster_top{i+1}" for i in range(args.n_designs - 1)] + ["single"]
    designs = [{"name": nm, "cdr_seq": sq, "mpnn_mean": None, "mpnn_worst": None,
                "wt_ddg": None, "baddg_grade": None} for nm, sq in zip(names, seqs)]
    out = dict(tmpl)
    out.update(antibody_designs=designs, oracle="baseline_rfantibody",
               generator="rfantibody", mechanism="baseline_denovo", base="rfantibody",
               template_json=os.path.basename(args.template), alignment_ok=True,
               designed_pdbs=kept[:len(designs)],
               note=("RFdiffusion-Ab designs a NEW CDR-H3 backbone and re-docks the Fv; this JSON "
                     "grades its SEQUENCE on the native complex. Own-structure grading uses "
                     "designed_pdbs."))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    wall = (time.time() - args.started) if args.started else None
    with open(args.out[:-5] + ".cost.json", "w") as fh:
        json.dump({"tag": "baseline/rfantibody", "target": target,
                   "wall_s": round(wall, 3) if wall else None, "cuda_s": 0.0,
                   "counters": {"designs_out": len(designs), "reward_evals": 0},
                   "sections": {}, "host": {"node": os.uname().nodename,
                                            "slurm_job": os.environ.get("SLURM_JOB_ID", "")},
                   "meta": {"n_designs": args.n_designs}}, fh, indent=2)
    print(f"[{target}] rfantibody: native={native_cdr} -> {seqs} -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
