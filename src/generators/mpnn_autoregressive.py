#!/usr/bin/env python3
"""ProteinMPNN-family AUTOREGRESSIVE baseline generator (the model's native inference mode).

Why this exists: sampling each CDR position independently from ProteinMPNN's marginal table is
NOT how ProteinMPNN is used -- it is an autoregressive decoder, and independent-marginal
sampling throws away all of its cross-position coupling. Grading our steered arm against a
crippled baseline would inflate our margin, so the "generator alone" arm has to use
`model.sample()`: sequential decoding over the CDR-H3 positions with the rest of the complex
(and the antigen) held fixed as context. AntiFold is deliberately NOT run this way -- its
published inference IS independent sampling from a per-position logit table, so its
independent-marginal baseline is the faithful one.

Pool convention matches every other arm exactly (same names, same panel, same pool size), so
graders and aggregators are untouched.

    python3 mpnn_ar_baseline.py --template T_template.json --features-dir <feats> \
        --proteinmpnn-dir vendor/ProteinMPNN --weights <ckpt> --model-name proteinmpnn_ar \
        --out T_designs_base_proteinmpnn_r1.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_steer import AA, N_AA, to_seq  # noqa: E402
from profile_util import Profiler  # noqa: E402

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--features-dir", required=True)
    ap.add_argument("--proteinmpnn-dir", required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--model-name", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-designs", type=int, default=6)
    ap.add_argument("--sampling-temp", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    tmpl = json.load(open(args.template))
    target, native_cdr = tmpl["target"], tmpl["native_cdr"]
    cdr_flat = tmpl["cdr_positions_flat"]
    prof = Profiler(tag=f"baseline_ar/{args.model_name}", target=target,
                    meta={"n_designs": args.n_designs, "sampling_temp": args.sampling_temp})

    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
    from lyra_mutants.Model.mpnn_adapter import ProteinMPNNAdapter

    device = "cuda" if torch.cuda.is_available() else "cpu"
    with prof.section("load"):
        adapter = ProteinMPNNAdapter(args.proteinmpnn_dir, args.weights,
                                     args.features_dir, cdr_flat).to(device)
        for p in adapter.model.parameters():
            p.requires_grad_(False)
        f = adapter._features(target)

    seq = to_seq(f["S"][0].tolist())
    struct_wt = "".join(seq[i] for i in cdr_flat)
    if struct_wt != native_cdr:
        print(f"[{target}] WARNING struct CDR {struct_wt!r} != template {native_cdr!r}", flush=True)

    L = f["X"].shape[1]
    # design ONLY the CDR-H3 positions; everything else (incl. the antigen) is fixed context
    chain_M_pos = torch.zeros(1, L, device=device)
    chain_M_pos[0, cdr_flat] = 1.0
    chain_M = torch.ones(1, L, device=device)
    omit_AAs_np = np.array([aa in "X" for aa in ALPHABET], dtype=np.float32)
    bias_AAs_np = np.zeros(len(ALPHABET), dtype=np.float32)
    zeros_res = torch.zeros(1, L, 21, device=device, dtype=torch.float32)
    pssm_coef = torch.zeros(1, L, device=device)
    pssm_bias = torch.zeros(1, L, 21, device=device)
    pssm_odds = torch.zeros(1, L, 21, device=device)

    torch.manual_seed(args.seed)
    seqs, seen = [], set()
    with prof.section("sample"):
        for _ in range(max(6 * args.n_designs, 40)):
            if len(seqs) >= args.n_designs:
                break
            randn = torch.randn(chain_M.shape, device=device)
            with torch.no_grad():
                out = adapter.model.sample(
                    f["X"], randn, f["S"], chain_M, f["chain_encoding_all"], f["residue_idx"],
                    mask=f["mask"], temperature=max(args.sampling_temp, 1e-3),
                    omit_AAs_np=omit_AAs_np, bias_AAs_np=bias_AAs_np, chain_M_pos=chain_M_pos,
                    omit_AA_mask=None, pssm_coef=pssm_coef, pssm_bias=pssm_bias, pssm_multi=0.0,
                    pssm_log_odds_flag=False, pssm_log_odds_mask=pssm_odds,
                    pssm_bias_flag=False, bias_by_res=zeros_res)
            prof.bump("base_forwards")
            S = out["S"][0]
            cdr = "".join(AA[int(S[i])] if int(S[i]) < N_AA else native_cdr[k]
                          for k, i in enumerate(cdr_flat))
            if cdr not in seen:
                seen.add(cdr)
                seqs.append(cdr)
    while len(seqs) < args.n_designs:                       # diversity exhausted: allow dups
        seqs.append(seqs[-1] if seqs else native_cdr)

    names = [f"cluster_top{i+1}" for i in range(args.n_designs - 1)] + ["single"]
    designs = [{"name": nm, "cdr_seq": sq, "mpnn_mean": None, "mpnn_worst": None,
                "wt_ddg": None, "baddg_grade": None} for nm, sq in zip(names, seqs)]
    out = dict(tmpl)
    out.update(antibody_designs=designs, oracle=f"baseline_{args.model_name}",
               generator=args.model_name, mechanism="baseline_autoregressive",
               base=args.model_name, sampling_temp=args.sampling_temp, seed=args.seed,
               weights=args.weights, template_json=os.path.basename(args.template),
               alignment_ok=(struct_wt == native_cdr),
               note="native autoregressive decoding (model.sample), CDR-H3 only")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"[{target}] {args.model_name} AR: native={native_cdr} -> {seqs} -> {args.out}", flush=True)
    prof.bump("designs_out", len(designs))
    prof.write(args.out)


if __name__ == "__main__":
    main()
