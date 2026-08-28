#!/usr/bin/env python3
"""DiffAb CDR-H3 BASELINE generator (antibody sequence-structure diffusion).

DiffAb co-designs the CDR-H3 sequence AND its backbone conformation, conditioned on the rest of
the Fv and on the antigen — a genuinely different generator family from the inverse-folding
models, which only ever emit a sequence for the fixed native backbone.

Two consequences worth stating plainly:
  * The CDR-H3 SEQUENCE it produces is graded like every other arm (threaded on the native
    complex), which keeps it in the same table but does not credit its designed conformation.
  * Its designed structures are kept in `--keep-pdb-dir` so the same designs can additionally be
    graded on their own backbones, which is the fair-to-DiffAb protocol.

Runs under the rebuilt env (DiffAb's shipped torch stops at sm_86 and hangs on H100 PTX JIT):
  <genzoo>/diffab_env/bin/python diffab_baseline.py --template T.json --out OUT.json
with ${VAAD_ROOT}/tools/diffab/env/bin on PATH (hmmscan for ANARCI renumbering).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
DIFFAB = f"{VAAD_ROOT}/tools/diffab"


def chain_seq(pdb_path, chain):
    seen, order = {}, []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[21] == chain:
                key = (int(line[22:26]), line[26].strip())
                if key not in seen:
                    seen[key] = AA3.get(line[17:20].strip(), "X")
                    order.append(key)
    return "".join(seen[k] for k in order)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--config", required=True, help="DiffAb yml restricted to H_CDR3")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--keep-pdb-dir", default=None,
                    help="copy the designed complexes here for own-structure grading")
    ap.add_argument("--n-designs", type=int, default=6)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    tmpl = json.load(open(args.template))
    target, native_cdr = tmpl["target"], tmpl["native_cdr"]
    ab = list(tmpl["antibody_chains"])
    heavy = tmpl["cdr_residues"][0]["chain"]
    light = ab[1] if len(ab) > 1 else None
    t0 = time.perf_counter()

    os.makedirs(args.workdir, exist_ok=True)
    cmd = [sys.executable, "design_pdb.py", os.path.join(args.pdb_dir, f"{target}.pdb"),
           "--heavy", heavy, "-c", args.config, "-o", args.workdir, "-t", target,
           "-b", str(args.n_designs), "-s", str(args.seed)]
    if light:
        cmd += ["--light", light]
    r = subprocess.run(cmd, cwd=DIFFAB, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-2000:] + "\n" + r.stderr[-2000:] + "\n")
        raise SystemExit(f"[{target}] DiffAb failed (rc={r.returncode})")

    runs = sorted(glob.glob(os.path.join(args.workdir, f"*{target}*", "*")) +
                  glob.glob(os.path.join(args.workdir, "*", f"{target}*")))
    rundir = None
    for c in sorted(glob.glob(os.path.join(args.workdir, "**", "reference.pdb"), recursive=True)):
        rundir = os.path.dirname(c)
    if rundir is None:
        raise SystemExit(f"[{target}] no DiffAb run directory under {args.workdir}")

    ref = chain_seq(os.path.join(rundir, "reference.pdb"), heavy)
    if native_cdr not in ref:
        raise SystemExit(f"[{target}] native CDR {native_cdr!r} not found in DiffAb reference "
                         f"heavy chain (renumbering mismatch)")
    if ref.count(native_cdr) != 1:
        raise SystemExit(f"[{target}] native CDR is ambiguous in the reference heavy chain")
    lo = ref.index(native_cdr)
    hi = lo + len(native_cdr)

    seqs, kept = [], []
    for f in sorted(glob.glob(os.path.join(rundir, "H_CDR3", "0*.pdb"))):
        s = chain_seq(f, heavy)
        if len(s) != len(ref):                      # length-changing design: skip, cannot thread
            continue
        cdr = s[lo:hi]
        if cdr and all(c in "ACDEFGHIKLMNPQRSTVWY" for c in cdr):
            seqs.append(cdr)
            kept.append(f)
        if len(seqs) >= args.n_designs:
            break
    if not seqs:
        raise SystemExit(f"[{target}] DiffAb produced no threadable CDR-H3")
    while len(seqs) < args.n_designs:
        seqs.append(seqs[-1])

    if args.keep_pdb_dir:                            # for the own-structure grading protocol
        os.makedirs(args.keep_pdb_dir, exist_ok=True)
        for i, f in enumerate(kept):
            shutil.copy(f, os.path.join(args.keep_pdb_dir, f"{target}_diffab_{i}.pdb"))

    names = [f"cluster_top{i+1}" for i in range(args.n_designs - 1)] + ["single"]
    designs = [{"name": nm, "cdr_seq": sq, "mpnn_mean": None, "mpnn_worst": None,
                "wt_ddg": None, "baddg_grade": None} for nm, sq in zip(names, seqs)]
    out = dict(tmpl)
    out.update(antibody_designs=designs, oracle="baseline_diffab", generator="diffab",
               mechanism="baseline_diffusion", base="diffab", seed=args.seed,
               template_json=os.path.basename(args.template), alignment_ok=True,
               note="DiffAb co-designs sequence+backbone; graded here on the native backbone")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    wall = time.perf_counter() - t0
    with open(args.out[:-5] + ".cost.json", "w") as fh:
        json.dump({"tag": "baseline/diffab", "target": target, "wall_s": round(wall, 3),
                   "cuda_s": 0.0, "counters": {"designs_out": len(designs), "reward_evals": 0},
                   "sections": {}, "host": {"node": os.uname().nodename,
                                            "slurm_job": os.environ.get("SLURM_JOB_ID", "")},
                   "meta": {"n_designs": args.n_designs}}, fh, indent=2)
    print(f"[{target}] diffab: native={native_cdr} -> {seqs} ({wall:.1f}s) -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
