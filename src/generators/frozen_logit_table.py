#!/usr/bin/env python3
"""Model-agnostic one-shot inverse-folding BASELINE generator from a base-logits .npz.

Every inverse-folding model in the zoo (AntiFold, ProteinMPNN, AbMPNN, LigandMPNN, ESM-IF1)
exports the same per-CDR-position logit table, so its unsteered baseline pool can be produced
by ONE script instead of one bespoke generator per model. That matters for fairness as much as
for effort: identical sampling convention (argmax first, then temperature sampling without
replacement), identical pool size, identical panel, identical output schema — so the only
thing that differs between two baseline arms is the model that produced the logits.

This is the "base alone" column of the base x {alone, +our objective} matrix; the "+ours"
column is design_for_validation.py --mechanism svdd --base <model> --base-logits <same npz>,
which samples from this very same table and only reweights it by the escape reward. The two
arms are therefore exactly matched on prior, panel, and pool size, and differ only in the
objective — which is the claim under test.

    python3 npz_baseline.py --template T_template.json --base-logits T_esmif_logits.npz \
        --model-name esmif --out T_designs_base_esmif_r1.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from profile_util import Profiler  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--base-logits", required=True)
    ap.add_argument("--model-name", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-designs", type=int, default=6)
    ap.add_argument("--sampling-temp", type=float, default=0.20,
                    help="AntiFold's own convention, kept identical across models")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    tmpl = json.load(open(args.template))
    target = tmpl["target"]
    native_cdr = tmpl["native_cdr"]
    n_cdr = len(native_cdr)
    prof = Profiler(tag=f"baseline/{args.model_name}", target=target,
                    meta={"n_designs": args.n_designs, "sampling_temp": args.sampling_temp})

    with prof.section("load_logits"):
        d = np.load(args.base_logits, allow_pickle=True)
        logits = d["logits"].astype(np.float64)
        usable = d["usable"].astype(bool)
        npz_native = str(d["native_cdr"])
    if logits.shape[0] != n_cdr:
        raise SystemExit(f"[{target}] npz n_cdr={logits.shape[0]} != template {n_cdr}")
    if npz_native != native_cdr:
        raise SystemExit(f"[{target}] npz native_cdr={npz_native!r} != template {native_cdr!r}")

    sample_idx = [i for i in range(n_cdr) if usable[i] and np.isfinite(logits[i]).all()]
    align_ok = bool(usable.all())
    print(f"[{target}] {args.model_name}: redesign {len(sample_idx)}/{n_cdr} positions "
          f"align_ok={align_ok}", flush=True)

    names = [f"cluster_top{i+1}" for i in range(args.n_designs - 1)] + ["single"]

    if not sample_idx:
        # Registration-broken target: emit an all-native pool (a tie vs native) and exit 0, so
        # the matched panel stays complete and downstream afterok chains are not broken.
        # Mirrors antifold_baseline.py's graceful degradation exactly.
        print(f"[{target}] WARNING no registerable CDR position; emitting all-native pool", flush=True)
        seqs = [native_cdr] * args.n_designs
    else:
        with prof.section("sample"):
            lg = torch.tensor(logits[sample_idx])
            probs = torch.softmax(lg / max(args.sampling_temp, 0.001), dim=1)
            torch.manual_seed(args.seed)

            def assemble(tokens):
                chars = list(native_cdr)
                for j, i in enumerate(sample_idx):
                    chars[i] = AA[tokens[j]]
                return "".join(chars)

            seqs, seen, tries = [], set(), 0
            s0 = assemble(lg.argmax(dim=1).tolist())
            seqs.append(s0); seen.add(s0)
            while len(seqs) < args.n_designs and tries < 2000:
                tries += 1
                s = assemble(torch.multinomial(probs, 1).squeeze(-1).tolist())
                if s not in seen:
                    seen.add(s); seqs.append(s)
            while len(seqs) < args.n_designs:                 # diversity exhausted: allow dups
                seqs.append(assemble(torch.multinomial(probs, 1).squeeze(-1).tolist()))
            prof.bump("base_forwards", 0)                     # table is precomputed; see export cost

    designs = [{"name": nm, "cdr_seq": sq, "mpnn_mean": None, "mpnn_worst": None,
                "wt_ddg": None, "baddg_grade": None} for nm, sq in zip(names, seqs)]
    out = dict(tmpl)
    out.update(antibody_designs=designs, oracle=f"baseline_{args.model_name}",
               generator=args.model_name, mechanism="baseline", base=args.model_name,
               sampling_temp=args.sampling_temp, seed=args.seed,
               base_logits=os.path.basename(args.base_logits),
               template_json=os.path.basename(args.template), alignment_ok=align_ok)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"[{target}] wrote {len(designs)} {args.model_name} designs -> {args.out}", flush=True)
    prof.bump("designs_out", len(designs))
    prof.write(args.out)


if __name__ == "__main__":
    main()
