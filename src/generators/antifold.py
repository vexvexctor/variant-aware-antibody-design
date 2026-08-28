#!/usr/bin/env python3
"""AntiFold CDR-H3 mutant-optimization BASELINE generator (BUCKET A).

Stands up AntiFold (OPIG antibody inverse-folding, ESM-IF1/GVP backbone fine-tuned on
SAbDab; antigen-conditioned via custom_chain_mode) as a SOTA baseline to compare against
Native and against our steered-diffusion (STEER) designs, graded by the SAME oracles
(FoldX / H3-DDG / Boltz-2 cofold).

Design choice — matched panel:
  For each target we MIRROR the existing STEER design JSON
  (results/validation/<target>_designs_f5_worst_wt1_r1.json): we reuse its cdr_residues,
  cdr_positions_flat, native_cdr, antibody/antigen chains and — crucially — its EXACT
  antigen_variants escape panel (WT + steer + held-out eval). We ONLY replace
  antibody_designs with AntiFold-sampled CDR-H3 loops. This guarantees the FoldX/H3/cofold
  graders run UNMODIFIED and the comparison is matched-budget (same pool size, same panel).

AntiFold usage:
  * build a one-row pdbs_csv (pdb, Hchain=antibody_chains[0], Lchain=antibody_chains[1],
    Agchain*=antigen chains) and run get_pdbs_logits(..., custom_chain_mode=True) on the
    AACDB complex PDB -> per-residue amino-acid logits over the WHOLE complex, so the CDR
    logits are conditioned on the antigen structure.
  * match the STEER JSON's cdr_residues (chain, resnum, icode) to the AntiFold df_logits
    rows (pdb_chain, pdb_pos, insertion code parsed from pdb_posins); verify wildtype ==
    native_cdr as an alignment guard (same guard the graders use).
  * temperature-sample --n-designs CDR loops with AntiFold's own convention
    softmax(logits / max(t,0.001)) (default t=0.20). First design = argmax (deterministic
    AntiFold-optimal), remaining sampled. Named cluster_top1..cluster_top5 + single to match
    the STEER pool schema exactly (graders + best-of-pool aggregators are name-driven).

Runs under the AntiFold venv:
  ${VAAD_ROOT}/tools/antifold/venv/bin/python antifold_baseline.py \
      --template  ${VAAD_ROOT}/results/validation/3SE9_HLG_designs_f5_worst_wt1_r1.json \
      --pdb-dir   ${VAAD_ROOT}/datasets/AACDB/data_zip/complex_structure \
      --out       ${VAAD_ROOT}/results/validation/3SE9_HLG_designs_antifold_r1.json
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd
import torch
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

AA = list("ACDEFGHIKLMNPQRSTVWY")


def parse_posins(posins):
    """'100A' -> (100, 'A'); '100' -> (100, '')."""
    s = str(posins).strip()
    m = re.match(r"^(-?\d+)([A-Za-z]?)$", s)
    if not m:
        return None, ""
    return int(m.group(1)), m.group(2).upper()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True,
                    help="STEER design JSON to mirror (schema + antigen panel + CDR residues)")
    ap.add_argument("--pdb-dir", required=True,
                    help="dir with <target>.pdb complex structures (AACDB)")
    ap.add_argument("--out", required=True, help="output design JSON path")
    ap.add_argument("--n-designs", type=int, default=6,
                    help="pool size to emit (matched to STEER: 5 cluster_top + 1 single)")
    ap.add_argument("--sampling-temp", type=float, default=0.20,
                    help="AntiFold CDR sampling temperature (their default 0.20)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--model-path", default="", help="AntiFold weights (default: models/model.pt)")
    args = ap.parse_args()

    from antifold.antiscripts import get_pdbs_logits, load_model

    tmpl = json.load(open(args.template))
    target = tmpl["target"]
    ab_chains = list(tmpl["antibody_chains"])
    ag_chains = list(tmpl["antigen_chains"])
    cdr_res = tmpl["cdr_residues"]
    native_cdr = tmpl["native_cdr"]
    pdb_name = tmpl.get("pdb", target + ".pdb")
    pdb_stem = os.path.splitext(pdb_name)[0]

    # --- build the one-row pdbs_csv: every non-pdb column must contain "chain" ---
    row = {"pdb": pdb_stem, "Hchain": ab_chains[0]}
    extra = []
    if len(ab_chains) > 1:
        row["Lchain"] = ab_chains[1]
        extra = ab_chains[2:] + ag_chains
    else:
        extra = ag_chains  # nanobody: only H; antigen chains follow
    for i, c in enumerate(extra):
        row[f"Agchain_{i}"] = c
    pdbs_csv = pd.DataFrame([row])
    print(f"[{target}] pdbs_csv columns={list(pdbs_csv.columns)} row={row}", flush=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[{target}] loading AntiFold model on {device} ...", flush=True)
    model = load_model(args.model_path)

    df_logits_list = get_pdbs_logits(
        model, pdbs_csv, args.pdb_dir, custom_chain_mode=True,
        batch_size=1, save_flag=False, seed=args.seed)
    dfl = df_logits_list[0]

    # --- match STEER cdr_residues -> AntiFold logit rows by (chain, resnum, icode) ---
    key2row = {}
    for _, r in dfl.iterrows():
        rn, ic = parse_posins(r["pdb_posins"])
        if rn is None:
            continue
        key2row[(str(r["pdb_chain"]), rn, ic)] = r

    # Only redesign CDR positions that (a) have a PDB resnum in the STEER JSON, (b) map to an
    # AntiFold logit row, AND (c) whose structural wildtype matches native_cdr. Positions that
    # are unmapped in the STEER JSON (resnum=None; a numbering gap) or missing from the AntiFold
    # logits are KEPT NATIVE — this is the matched behavior: the shared FoldX grader skips
    # mutations at unmapped residues and the H3 grader errors the whole target for both arms, so
    # redesigning only the structurally-present CDR positions keeps the comparison apples-to-apples.
    sampleable, struct_wt = [], []          # per position: logit-vector-or-None, structural WT char
    n_unmapped, n_wtmismatch = 0, 0
    for i, cr in enumerate(cdr_res):
        rn = cr.get("resnum")
        nat = native_cdr[i]
        if rn is None:
            sampleable.append(None); struct_wt.append("-"); n_unmapped += 1; continue
        key = (str(cr["chain"]), int(rn), (cr.get("icode") or "").strip().upper())
        r = key2row.get(key)
        if r is None:
            sampleable.append(None); struct_wt.append("?"); n_unmapped += 1; continue
        wt = str(r["pdb_res"])
        struct_wt.append(wt)
        if wt != nat:
            sampleable.append(None); n_wtmismatch += 1; continue
        sampleable.append(np.asarray([float(r[a]) for a in AA], dtype=np.float64))

    sample_idx = [i for i, v in enumerate(sampleable) if v is not None]
    align_ok = (n_unmapped == 0 and n_wtmismatch == 0)
    print(f"[{target}] native_cdr(json)   = {native_cdr}", flush=True)
    print(f"[{target}] native_cdr(struct) = {''.join(struct_wt)}  align_ok={align_ok} "
          f"redesign {len(sample_idx)}/{len(cdr_res)} positions "
          f"(unmapped={n_unmapped}, wt_mismatch={n_wtmismatch})", flush=True)
    if not sample_idx:
        # Registration-broken target (CDR map disagrees with native_cdr at every position — the
        # "numbering differs / needs ANARCI" gotcha). The shared FoldX/H3 graders can't grade this
        # for the STEER arm either (native_cdr is used as the FoldX mutation wt). Degrade gracefully:
        # emit an all-native pool (a tie vs native) and exit 0 so the matched panel stays complete
        # and the grading chain's afterok dependency is not broken. Flagged align_ok=False.
        print(f"[{target}] WARNING: no CDR position is registerable "
              f"(unmapped={n_unmapped}, wt_mismatch={n_wtmismatch}); emitting all-native pool.",
              flush=True)
        names = [f"cluster_top{i+1}" for i in range(args.n_designs - 1)] + ["single"]
        designs = [{"name": nm, "cdr_seq": native_cdr, "mpnn_mean": None, "mpnn_worst": None,
                    "wt_ddg": None, "baddg_grade": None} for nm in names]
        out = dict(tmpl)
        out.update(antibody_designs=designs, oracle="antifold", generator="antifold-0.3.1",
                   antifold_sampling_temp=args.sampling_temp, antifold_seed=args.seed,
                   template_json=os.path.basename(args.template), alignment_ok=False,
                   antifold_note="registration_failed_all_native")
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"wrote all-native fallback pool -> {args.out}", flush=True)
        return

    # --- temperature-sample the pool: design 0 = argmax (AntiFold-optimal), rest sampled ---
    logits = torch.tensor(np.stack([sampleable[i] for i in sample_idx]))   # (n_sample, 20)
    temp_probs = torch.softmax(logits / max(args.sampling_temp, 0.001), dim=1)
    torch.manual_seed(args.seed)

    def _assemble(tokens):
        """tokens: index per sampleable position -> full-length CDR (native at kept positions)."""
        chars = list(native_cdr)
        for j, i in enumerate(sample_idx):
            chars[i] = AA[tokens[j]]
        return "".join(chars)

    seqs, seen = [], set()
    argmax_seq = _assemble(logits.argmax(dim=1).tolist())
    seqs.append(argmax_seq)
    seen.add(argmax_seq)
    tries = 0
    while len(seqs) < args.n_designs and tries < 2000:
        tries += 1
        s = _assemble(torch.multinomial(temp_probs, 1).squeeze(-1).tolist())
        if s not in seen:
            seen.add(s)
            seqs.append(s)
    # if diversity exhausted (short/near-deterministic CDR), pad with resampled (dups allowed)
    while len(seqs) < args.n_designs:
        seqs.append(_assemble(torch.multinomial(temp_probs, 1).squeeze(-1).tolist()))

    names = [f"cluster_top{i+1}" for i in range(args.n_designs - 1)] + ["single"]
    designs = [{"name": nm, "cdr_seq": sq,
                "mpnn_mean": None, "mpnn_worst": None, "wt_ddg": None, "baddg_grade": None}
               for nm, sq in zip(names, seqs)]

    out = dict(tmpl)  # mirror the STEER schema (antigen panel, cdr map, chains, etc.)
    out["antibody_designs"] = designs
    out["oracle"] = "antifold"
    out["generator"] = "antifold-0.3.1"
    out["antifold_sampling_temp"] = args.sampling_temp
    out["antifold_seed"] = args.seed
    out["template_json"] = os.path.basename(args.template)
    out["alignment_ok"] = bool(tmpl.get("alignment_ok", True) and align_ok)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)

    print(f"\n[{target}] native  {native_cdr}")
    for d in designs:
        nmis = sum(a != b for a, b in zip(d["cdr_seq"], native_cdr))
        print(f"  {d['name']:<13} {d['cdr_seq']}  ({nmis} muts vs native)")
    print(f"\nwrote {len(designs)} AntiFold designs x {len(out['antigen_variants'])} "
          f"antigen variants -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
