#!/usr/bin/env python3
"""Export AntiFold per-CDR-position logits to an .npz for use as a steering PRIOR.

This is the OFFLINE half of AntiFold-as-base-distribution (product-of-experts steering).
It runs under the AntiFold venv, extracts AntiFold's antigen-conditioned amino-acid logits
for exactly the CDR positions in a STEER design JSON, row-aligns them to that JSON's
`cdr_residues` order, and dumps them to disk. The steering process (a different venv:
torch + ProteinMPNN + BA-DDG) then loads the .npz as plain numpy and adds an AntiFold
log-prob term to the diffusion energy WITHOUT importing AntiFold. This sidesteps the
venv/dependency split that blocks an in-process import.

The extraction logic MIRRORS antifold_baseline.py (same get_pdbs_logits call, same
(chain,resnum,icode) matching, same wt==native_cdr alignment guard) so the exported table
is row-consistent with the baseline generator and with the FoldX/H3 graders.

Output .npz (keyed by the STEER JSON basename, written next to it by convention):
  logits   : float64 [n_cdr, 20]  raw AntiFold logits, AA order "ACDEFGHIKLMNPQRSTVWY",
                                   row i == cdr_residues[i]; UNUSABLE rows are NaN.
  usable   : bool    [n_cdr]       True where the position mapped to an AntiFold row AND
                                   its structural WT matches native_cdr (else the prior
                                   must contribute 0 at that position).
  native_cdr : str                 from the template (json)
  struct_wt  : str                 structural WT per position (AntiFold pdb_res; '-' unmapped,
                                   '?' missing-from-logits)
  aa_order   : str                 "ACDEFGHIKLMNPQRSTVWY"
  align_ok   : bool                all positions usable
  target     : str

Runs under the AntiFold venv:
  ${VAAD_ROOT}/tools/antifold/venv/bin/python antifold_logits_export.py \
      --template ${VAAD_ROOT}/results/validation/3SE9_HLG_designs_f5_worst_wt1_r1.json \
      --pdb-dir  ${VAAD_ROOT}/datasets/AACDB/data_zip/complex_structure \
      --out      ${VAAD_ROOT}/results/validation/3SE9_HLG_antifold_logits.npz
"""
from __future__ import annotations

import argparse
import json
import os
import re

import numpy as np
import pandas as pd
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

AA = list("ACDEFGHIKLMNPQRSTVWY")
AA_ORDER = "".join(AA)


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
                    help="STEER design JSON to mirror (cdr_residues + native_cdr + chains)")
    ap.add_argument("--pdb-dir", required=True,
                    help="dir with <target>.pdb complex structures (AACDB)")
    ap.add_argument("--out", required=True, help="output .npz path")
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
    n_cdr = len(cdr_res)

    # --- build the one-row pdbs_csv: every non-pdb column must contain "chain" ---
    row = {"pdb": pdb_stem, "Hchain": ab_chains[0]}
    if len(ab_chains) > 1:
        row["Lchain"] = ab_chains[1]
        extra = ab_chains[2:] + ag_chains
    else:
        extra = ag_chains  # nanobody: only H; antigen chains follow
    for i, c in enumerate(extra):
        row[f"Agchain_{i}"] = c
    pdbs_csv = pd.DataFrame([row])

    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[{target}] loading AntiFold on {device}; pdbs_csv row={row}", flush=True)
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

    # Row i aligns to cdr_res[i]. A position is USABLE only if it maps to an AntiFold row
    # AND its structural WT matches native_cdr[i] (same guard the FoldX/H3 graders apply).
    # Unusable rows -> NaN logits + usable=False so the steering prior contributes 0 there.
    logits = np.full((n_cdr, len(AA)), np.nan, dtype=np.float64)
    usable = np.zeros(n_cdr, dtype=bool)
    struct_wt = []
    n_unmapped, n_wtmismatch = 0, 0
    for i, cr in enumerate(cdr_res):
        rn = cr.get("resnum")
        nat = native_cdr[i]
        if rn is None:
            struct_wt.append("-"); n_unmapped += 1; continue
        key = (str(cr["chain"]), int(rn), (cr.get("icode") or "").strip().upper())
        r = key2row.get(key)
        if r is None:
            struct_wt.append("?"); n_unmapped += 1; continue
        wt = str(r["pdb_res"])
        struct_wt.append(wt)
        if wt != nat:
            n_wtmismatch += 1; continue
        logits[i] = np.asarray([float(r[a]) for a in AA], dtype=np.float64)
        usable[i] = True

    struct_wt = "".join(struct_wt)
    align_ok = bool(n_unmapped == 0 and n_wtmismatch == 0)
    print(f"[{target}] native_cdr(json)   = {native_cdr}", flush=True)
    print(f"[{target}] native_cdr(struct) = {struct_wt}  align_ok={align_ok}  "
          f"usable {int(usable.sum())}/{n_cdr} "
          f"(unmapped={n_unmapped}, wt_mismatch={n_wtmismatch})", flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    np.savez(
        args.out,
        logits=logits,
        usable=usable,
        native_cdr=native_cdr,
        struct_wt=struct_wt,
        aa_order=AA_ORDER,
        align_ok=align_ok,
        target=target,
    )
    print(f"[{target}] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
