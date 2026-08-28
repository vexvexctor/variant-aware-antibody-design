#!/usr/bin/env python3
"""ESM-IF1 constrained AUTOREGRESSIVE baseline generator (the model's native inference mode).

ESM-IF1 decodes sequences autoregressively, so sampling CDR-H3 positions independently from a
marginal table discards the coupling the model actually provides and understates it as a
baseline. Here the heavy chain is decoded left-to-right with every non-CDR residue FORCED to
its native identity and only the CDR-H3 positions sampled -- the standard constrained-design
use of an autoregressive inverse-folding model, with the light chain and antigen as complex
context.

Only n_cdr forward passes are needed, not one per residue: every position before the CDR is
pinned to native, so the prefix is fixed and each pass yields the logits for the next CDR slot.

Runs under if_env:
  <if_env>/bin/python esmif_ar_baseline.py --template T.json --pdb-dir <pdbs> --out OUT.json
"""
from __future__ import annotations
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

import argparse
import json
import os
import time

import numpy as np
import torch

AA = "ACDEFGHIKLMNPQRSTVWY"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-designs", type=int, default=6)
    ap.add_argument("--sampling-temp", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--torch-hub", default=f"{VAAD_ROOT}/scratch/hf_home/torchhub")
    args = ap.parse_args()

    os.environ.setdefault("TORCH_HOME", args.torch_hub)
    import esm
    import esm.inverse_folding as inv
    import biotite.structure as bs

    tmpl = json.load(open(args.template))
    target, native_cdr = tmpl["target"], tmpl["native_cdr"]
    cdr_res = tmpl["cdr_residues"]
    n_cdr = len(cdr_res)
    chains = list(tmpl["antibody_chains"]) + list(tmpl["antigen_chains"])
    h_chain = cdr_res[0]["chain"]
    t0 = time.perf_counter()

    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.eval().to(device)

    structure = inv.util.load_structure(os.path.join(args.pdb_dir, f"{target}.pdb"), chains)
    coords, native_seqs = inv.multichain_util.extract_coords_from_complex(structure)
    all_coords = inv.multichain_util._concatenate_coords(coords, h_chain)
    h_seq = native_seqs[h_chain]

    sub = structure[structure.chain_id == h_chain]
    keys = [(int(sub.res_id[i]), str(sub.ins_code[i]).strip()) for i in bs.get_residue_starts(sub)]
    k2i = {k: i for i, k in enumerate(keys)}
    slots = []
    for i, r in enumerate(cdr_res):
        rn = r.get("resnum")
        j = None if rn is None else k2i.get((int(rn), str(r.get("icode") or "").strip()))
        if j is not None and h_seq[j] != native_cdr[i]:
            j = None
        slots.append(j)

    cols = [alphabet.get_idx(a) for a in AA]
    bc = inv.util.CoordBatchConverter(alphabet)
    torch.manual_seed(args.seed)

    def sample_one():
        seq = list(h_seq)
        out = list(native_cdr)
        for i, j in enumerate(slots):
            if j is None:                                   # unregisterable slot: keep native
                continue
            c, conf, _s, tok, pad = bc([(all_coords, None, "".join(seq))], device=device)
            with torch.no_grad():
                lg, _ = model.forward(c, pad, conf, tok[:, :-1])
            logits = lg[0, :, j].float()
            p = torch.softmax(logits[cols] / max(args.sampling_temp, 1e-3), dim=0)
            aa = AA[int(torch.multinomial(p, 1))]
            out[i] = aa
            seq[j] = aa                                     # commit -> conditions the next slot
        return "".join(out)

    seqs, seen = [], set()
    for _ in range(max(3 * args.n_designs, 18)):
        if len(seqs) >= args.n_designs:
            break
        s = sample_one()
        if s not in seen:
            seen.add(s); seqs.append(s)
    while len(seqs) < args.n_designs:
        seqs.append(seqs[-1] if seqs else native_cdr)

    names = [f"cluster_top{i+1}" for i in range(args.n_designs - 1)] + ["single"]
    designs = [{"name": nm, "cdr_seq": sq, "mpnn_mean": None, "mpnn_worst": None,
                "wt_ddg": None, "baddg_grade": None} for nm, sq in zip(names, seqs)]
    o = dict(tmpl)
    o.update(antibody_designs=designs, oracle="baseline_esmif", generator="esmif",
             mechanism="baseline_autoregressive", base="esmif",
             sampling_temp=args.sampling_temp, seed=args.seed,
             template_json=os.path.basename(args.template),
             alignment_ok=all(j is not None for j in slots),
             note="constrained autoregressive decoding: non-CDR residues forced to native")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(o, fh, indent=2)

    wall = time.perf_counter() - t0
    with open(args.out[:-5] + ".cost.json", "w") as fh:
        json.dump({"tag": "baseline_ar/esmif", "target": target, "wall_s": round(wall, 3),
                   "cuda_s": 0.0, "counters": {"designs_out": len(designs), "reward_evals": 0,
                                               "base_forwards": len(seqs) * n_cdr},
                   "host": {"node": os.uname().nodename,
                            "slurm_job": os.environ.get("SLURM_JOB_ID", "")},
                   "sections": {}, "meta": {"n_designs": args.n_designs}}, fh, indent=2)
    print(f"[{target}] esmif AR: native={native_cdr} -> {seqs} ({wall:.1f}s) -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
