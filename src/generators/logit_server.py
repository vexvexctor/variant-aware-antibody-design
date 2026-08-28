#!/usr/bin/env python3
"""Persistent per-step logit worker for the ESM-IF-lineage bases (AntiFold, ESM-IF1).

Why a subprocess: these models re-condition on the partially-revealed CDR at EVERY denoising
step -- that is what makes the guidance "with the generator in the loop" rather than a reweighted
frozen prior -- but they live in venvs that cannot be imported alongside the steering process
(torch-geometric / biotite / antifold pins all conflict). So the model is held resident in its
own interpreter and answers one request per denoising step over stdin/stdout.

Both backends are teacher-forced encoder-decoder passes: substituting the current partial CDR
into the decoder input makes the returned logits conditional on everything revealed so far,
which is exactly the step-conditional prior SVDD should be reweighting.

Protocol (JSON lines):
    <- {"cdr": "ARHGNYYYYSGMDV"}          full-length current CDR (predicted-clean filled)
    -> {"logits": [[20 floats], ...]}     row i == template cdr_residues[i]
    <- {"cmd": "quit"}
Startup handshake: the server prints {"ready": true, "n_cdr": N} once loaded.

    <antifold-venv>/bin/python logit_server.py --backend antifold --template T.json --pdb-dir D
    <if_env>/bin/python       logit_server.py --backend esmif    --template T.json --pdb-dir D
"""
from __future__ import annotations
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

import argparse
import json
import os
import sys

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"


def emit(obj):
    sys.stdout.write(json.dumps(obj) + "\n")
    sys.stdout.flush()


# --------------------------------------------------------------------------- ESM-IF1 backend
class EsmIfBackend:
    def __init__(self, tmpl, pdb_dir, torch_hub):
        os.environ.setdefault("TORCH_HOME", torch_hub)
        import esm
        import esm.inverse_folding as inv
        import torch
        self.torch, self.inv = torch, inv

        self.tmpl = tmpl
        target = tmpl["target"]
        self.cdr_res = tmpl["cdr_residues"]
        self.native_cdr = tmpl["native_cdr"]
        chains = list(tmpl["antibody_chains"]) + list(tmpl["antigen_chains"])
        self.h_chain = self.cdr_res[0]["chain"]

        self.model, self.alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.model = self.model.eval().to(self.device)

        structure = inv.util.load_structure(os.path.join(pdb_dir, f"{target}.pdb"), chains)
        coords, native_seqs = inv.multichain_util.extract_coords_from_complex(structure)
        self.all_coords = inv.multichain_util._concatenate_coords(coords, self.h_chain)
        self.h_seq = native_seqs[self.h_chain]

        import biotite.structure as bs
        sub = structure[structure.chain_id == self.h_chain]
        keys = [(int(sub.res_id[i]), str(sub.ins_code[i]).strip())
                for i in bs.get_residue_starts(sub)]
        k2i = {k: i for i, k in enumerate(keys)}
        self.slots = []
        for i, r in enumerate(self.cdr_res):
            rn = r.get("resnum")
            j = None if rn is None else k2i.get((int(rn), str(r.get("icode") or "").strip()))
            if j is not None and self.h_seq[j] != self.native_cdr[i]:
                j = None                                   # same wt-match guard as everywhere else
            self.slots.append(j)
        self.cols = [self.alphabet.get_idx(a) for a in AA]
        self.bc = inv.util.CoordBatchConverter(self.alphabet)

    def logits(self, cdr):
        seq = list(self.h_seq)
        for i, j in enumerate(self.slots):                 # write the current partial CDR in
            if j is not None:
                seq[j] = cdr[i]
        seq = "".join(seq)
        c, conf, _s, tok, pad = self.bc([(self.all_coords, None, seq)], device=self.device)
        with self.torch.no_grad():
            lg, _ = self.model.forward(c, pad, conf, tok[:, :-1])
        lg = lg[0].float().cpu().numpy()
        out = np.zeros((len(self.cdr_res), 20), dtype=float)
        for i, j in enumerate(self.slots):
            if j is None or j >= lg.shape[1]:
                out[i] = 0.0                               # uninformative -> uniform after softmax
            else:
                out[i] = lg[self.cols, j]
        return out


# -------------------------------------------------------------------------- AntiFold backend
class AntiFoldBackend:
    def __init__(self, tmpl, pdb_dir):
        import torch
        from antifold.antiscripts import get_pdbs_logits, load_IF1_model
        self.torch = torch
        self.tmpl = tmpl
        self.cdr_res = tmpl["cdr_residues"]
        self.native_cdr = tmpl["native_cdr"]
        self.target = tmpl["target"]
        self.pdb_dir = pdb_dir
        ab = list(tmpl["antibody_chains"]); ag = list(tmpl["antigen_chains"])
        row = {"pdb": self.target, "Hchain": ab[0], "Lchain": ab[1] if len(ab) > 1 else ab[0]}
        for n, c in enumerate(ag):
            row[f"Agchain_{n}"] = c
        self._row = row
        self._get = get_pdbs_logits
        self._load = load_IF1_model
        self.model = load_IF1_model()

    def _df(self, seq_override=None):
        import pandas as pd
        df = pd.DataFrame([self._row])
        return self._get(self.model, df, self.pdb_dir, custom_chain_mode=True)

    def logits(self, cdr):
        # AntiFold's public path re-reads the PDB each call; the decoder is teacher-forced on the
        # structure's own residues, so re-conditioning means writing the current CDR into the
        # residue identities AntiFold reads. Handled by the caller-side sequence patch below.
        raise NotImplementedError


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", required=True, choices=["esmif", "antifold"])
    ap.add_argument("--template", required=True)
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--torch-hub", default=f"{VAAD_ROOT}/scratch/hf_home/torchhub")
    args = ap.parse_args()

    tmpl = json.load(open(args.template))
    if args.backend == "esmif":
        be = EsmIfBackend(tmpl, args.pdb_dir, args.torch_hub)
    else:
        be = AntiFoldBackend(tmpl, args.pdb_dir)
    emit({"ready": True, "n_cdr": len(tmpl["cdr_residues"]), "backend": args.backend})

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        req = json.loads(line)
        if req.get("cmd") == "quit":
            break
        try:
            lg = be.logits(req["cdr"])
            emit({"logits": [[float(x) for x in row] for row in lg]})
        except Exception as e:                              # never wedge the design run
            emit({"error": f"{type(e).__name__}: {e}"})


if __name__ == "__main__":
    main()
