#!/usr/bin/env python3
"""Export ESM-IF1 per-CDR-position logits to the shared base-logits .npz.

ESM-IF1 (esm_if1_gvp4_t16_142M_UR50) is an autoregressive structure->sequence model, so its
"per-position distribution" is the next-token distribution under teacher forcing on the native
sequence — the closest analogue to the one-shot marginals AntiFold/ProteinMPNN expose, and the
same quantity its own scoring API integrates. The antibody heavy chain is the target chain and
the light chain + antigen enter as complex context (multichain_util._concatenate_coords), so
ESM-IF1 sees the same information as the other bases.

Runs under its own venv (torch-geometric/torch-sparse/torch-scatter conflict with lyra_env):
  ${VAAD_ROOT}/scratch/genzoo/if_env/bin/python export_esmif_logits.py \
      --template T_template.json --pdb-dir <AACDB> --out T_esmif_logits.npz
"""
from __future__ import annotations
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

import argparse
import json
import os

import numpy as np
import torch

AA = "ACDEFGHIKLMNPQRSTVWY"


def residue_keys(structure, chain_id):
    """(res_id, ins_code) per residue of `chain_id`, in the exact order ESM-IF1's coord/seq
    extraction produces — i.e. over the backbone-filtered atom array."""
    import biotite.structure as bs
    sub = structure[structure.chain_id == chain_id]
    starts = bs.get_residue_starts(sub)
    return [(int(sub.res_id[i]), str(sub.ins_code[i]).strip()) for i in starts]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--torch-hub", default=f"{VAAD_ROOT}/scratch/hf_home/torchhub")
    args = ap.parse_args()

    os.environ.setdefault("TORCH_HOME", args.torch_hub)
    import esm
    import esm.inverse_folding as inv

    tmpl = json.load(open(args.template))
    target = tmpl["target"]
    cdr_res = tmpl["cdr_residues"]
    native_cdr = tmpl["native_cdr"]
    n_cdr = len(cdr_res)
    chains = list(tmpl["antibody_chains"]) + list(tmpl["antigen_chains"])
    h_chain = cdr_res[0]["chain"]
    pdb = os.path.join(args.pdb_dir, f"{target}.pdb")

    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)

    structure = inv.util.load_structure(pdb, chains)
    coords, native_seqs = inv.multichain_util.extract_coords_from_complex(structure)
    if h_chain not in coords:
        raise SystemExit(f"[{target}] heavy chain {h_chain} absent from {pdb} (chains={list(coords)})")

    keys = residue_keys(structure, h_chain)
    h_seq = native_seqs[h_chain]
    if len(keys) != len(h_seq):
        raise SystemExit(f"[{target}] residue-key/seq length mismatch {len(keys)} vs {len(h_seq)}")
    key_to_i = {k: i for i, k in enumerate(keys)}

    all_coords = inv.multichain_util._concatenate_coords(coords, h_chain)
    batch_converter = inv.util.CoordBatchConverter(alphabet)
    bc, confidence, strs, tokens, padding_mask = batch_converter(
        [(all_coords, None, h_seq)], device=device)
    with torch.no_grad():
        logits, _ = model.forward(bc, padding_mask, confidence, tokens[:, :-1])
    logits = logits[0].float().cpu().numpy()               # [vocab, L]

    cols = [alphabet.get_idx(a) for a in AA]
    out = np.full((n_cdr, 20), np.nan, dtype=np.float64)
    usable = np.zeros(n_cdr, dtype=bool)
    struct_wt = []
    for i, r in enumerate(cdr_res):
        rn = r.get("resnum")
        if rn is None:                                     # unmapped in the template's CDR map
            struct_wt.append("-")
            continue
        k = (int(rn), str(r.get("icode") or "").strip())
        j = key_to_i.get(k)
        if j is None or j >= logits.shape[1]:
            struct_wt.append("-")
            continue
        struct_wt.append(h_seq[j])
        if h_seq[j] != native_cdr[i]:                      # same wt-match guard as every other base
            continue
        out[i] = logits[cols, j]
        usable[i] = True
    struct_wt = "".join(struct_wt)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez(args.out, logits=out, usable=usable, native_cdr=native_cdr, struct_wt=struct_wt,
             aa_order=AA, align_ok=bool(usable.all()), target=target, model="esmif")
    print(f"[{target}] esmif: native={native_cdr} struct={struct_wt} "
          f"usable {int(usable.sum())}/{n_cdr} -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
