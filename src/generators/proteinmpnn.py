#!/usr/bin/env python3
"""ProteinMPNN-family CDR-H3 BASELINE generator (vanilla-ProteinMPNN / AbMPNN).

Mirror of antifold_baseline.py's I/O + schema, but the CDR-H3 loops are produced by
fixed-backbone, antigen-conditioned ProteinMPNN sampling instead of AntiFold logits.

Design choice — matched panel (identical to antifold_baseline.py):
  For each target we MIRROR the template design JSON
  (results/validation/<target>_designs_refold_refold_r1.json): reuse its cdr_residues,
  cdr_positions_flat, native_cdr, antibody/antigen chains and its EXACT antigen_variants
  escape panel (WT + steer + held-out eval, split-tagged). We ONLY replace antibody_designs
  with ProteinMPNN-sampled CDR-H3 loops. The FoldX/H3 graders then run UNMODIFIED and the
  comparison is matched-budget (same panel, same pool size, same native reference).

ProteinMPNN usage (standard "redesign a loop" recipe, fixed backbone):
  * parse the AACDB complex PDB with ProteinMPNN's own parse_PDB.
  * mark the antibody heavy chain (antibody_chains[0]) as the single DESIGN (masked) chain and
    every other chain (light chain + antigen chain(s)) as VISIBLE (fixed to native sequence).
    Autoregressive decoding therefore conditions the CDR logits on the antigen structure AND
    the fixed native sequence of the rest of the complex.
  * fix EVERY heavy-chain position except the CDR-H3 residues via fixed_positions_dict, so only
    the CDR-H3 loop is resampled (everything else stays native, exactly like the AntiFold arm).
  * temperature-sample --n-designs loops: design 0 = near-argmax (temp 1e-5, ProteinMPNN-optimal),
    the rest sampled at --sampling-temp (default 0.20). Named cluster_top1..cluster_top{n-1} +
    single to match the pool schema (graders + best-of-pool aggregators are name-driven).

Registration guard (same intent as AntiFold arm): the heavy chain is re-numbered exactly the way
ProteinMPNN's parse_PDB_biounits numbers residues, so each template cdr_residue (chain,resnum,icode)
maps to a concrete flattened index; the structural wildtype at that index must equal native_cdr[i].
Positions that don't map or whose WT disagrees are KEPT NATIVE. If NO CDR position is registerable
we emit an all-native pool (a tie vs native) and exit 0 so the matched panel stays complete.

Runs under any torch env (system python3 or the antifold venv); ProteinMPNN needs only numpy+torch:
  python3 scripts/nos_diffusion/proteinmpnn_baseline.py \
      --template /data/.../<target>_designs_refold_refold_r1.json \
      --pdb-dir  /data/.../AACDB/data_zip/complex_structure \
      --out      /data/.../<target>_designs_mpnn_r1.json \
      --model-path /data/.../vendor/ProteinMPNN/vanilla_model_weights/v_48_002.pt \
      --n-designs 6 --sampling-temp 0.20 --seed 42
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import torch

# ProteinMPNN vendor dir (protein_mpnn_utils importable) -----------------------------------------
MPNN_DIR = os.environ.get(
    "PROTEINMPNN_DIR",
    f"{VAAD_ROOT}/projects/mutation_sampling/lyra/vendor/ProteinMPNN",
)
if MPNN_DIR not in sys.path:
    sys.path.insert(0, MPNN_DIR)

from protein_mpnn_utils import ProteinMPNN, parse_PDB, tied_featurize  # noqa: E402
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

# ProteinMPNN's amino-acid alphabet (21st = X). Index of a residue in S maps through this.
MPNN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"


def chain_numbering(pdb_path: str, chain: str):
    """Replicate ProteinMPNN parse_PDB_biounits numbering for one chain.

    Returns a list parallel to the parsed chain sequence: each entry is (pdb_resnum, icode) for a
    real residue, or (None, "") for a gap-filler that parse_PDB inserts for missing residue numbers.
    parse_PDB_biounits stores resn internally as int(pdb_resnum)-1, iterates range(min,max+1), and
    for each present resn appends residues in sorted-insertion-code order (gap => single filler).
    """
    xyz_resn = {}  # internal_resn -> {icode: resname3}
    min_resn, max_resn = 10**9, -(10**9)
    with open(pdb_path, "rb") as fh:
        for raw in fh:
            line = raw.decode("utf-8", "ignore").rstrip()
            if line[:6] == "HETATM" and line[17:20] == "MSE":
                line = line.replace("HETATM", "ATOM  ").replace("MSE", "MET")
            if line[:4] != "ATOM":
                continue
            if line[21:22] != chain:
                continue
            resn_field = line[22:27].strip()
            if not resn_field:
                continue
            if resn_field[-1].isalpha():
                resa = resn_field[-1]
                resn = int(resn_field[:-1]) - 1
            else:
                resa = ""
                resn = int(resn_field) - 1
            min_resn = min(min_resn, resn)
            max_resn = max(max_resn, resn)
            xyz_resn.setdefault(resn, {}).setdefault(resa, line[17:20])
    order = []
    if max_resn < min_resn:
        return order
    for resn in range(min_resn, max_resn + 1):
        if resn in xyz_resn:
            for icode in sorted(xyz_resn[resn]):
                order.append((resn + 1, icode))  # back to 1-based PDB numbering
        else:
            order.append((None, ""))  # gap filler (parse_PDB emits a '-'/X here)
    return order


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True,
                    help="template design JSON to mirror (schema + antigen panel + CDR residues)")
    ap.add_argument("--pdb-dir", required=True, help="dir with <target>.pdb complex structures (AACDB)")
    ap.add_argument("--out", required=True, help="output design JSON path")
    ap.add_argument("--model-path", required=True, help="ProteinMPNN checkpoint .pt (vanilla or AbMPNN)")
    ap.add_argument("--n-designs", type=int, default=6, help="pool size to emit (5 cluster_top + 1 single)")
    ap.add_argument("--sampling-temp", type=float, default=0.20, help="ProteinMPNN sampling temperature")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--oracle-label", default="proteinmpnn", help="metadata oracle tag (graders ignore)")
    ap.add_argument("--generator-label", default="proteinmpnn-vanilla", help="metadata generator tag")
    args = ap.parse_args()

    tmpl = json.load(open(args.template))
    target = tmpl["target"]
    ab_chains = list(tmpl["antibody_chains"])
    ag_chains = list(tmpl["antigen_chains"])
    cdr_res = tmpl["cdr_residues"]
    native_cdr = tmpl["native_cdr"]
    pdb_name = tmpl.get("pdb", target + ".pdb")
    pdb_path = os.path.join(args.pdb_dir, pdb_name)

    names = [f"cluster_top{i+1}" for i in range(args.n_designs - 1)] + ["single"]

    def write_out(designs, align_ok, note=None):
        out = dict(tmpl)
        out["antibody_designs"] = designs
        out["oracle"] = args.oracle_label
        out["generator"] = args.generator_label
        out["mpnn_model"] = os.path.basename(args.model_path)
        out["mpnn_sampling_temp"] = args.sampling_temp
        out["mpnn_seed"] = args.seed
        out["template_json"] = os.path.basename(args.template)
        out["alignment_ok"] = bool(tmpl.get("alignment_ok", True) and align_ok)
        if note:
            out["mpnn_note"] = note
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=2)

    def all_native_pool(note):
        designs = [{"name": nm, "cdr_seq": native_cdr, "mpnn_mean": None, "mpnn_worst": None,
                    "wt_ddg": None, "baddg_grade": None} for nm in names]
        write_out(designs, align_ok=False, note=note)

    # --- parse the complex with ProteinMPNN; heavy chain designed, rest visible/native ----------
    h_chain = ab_chains[0]
    other_chains = ab_chains[1:] + ag_chains
    want_chains = [h_chain] + other_chains
    pdb_dicts = parse_PDB(pdb_path, input_chain_list=want_chains)
    pdb = pdb_dicts[0]
    name = pdb["name"]

    parsed_h = pdb.get(f"seq_chain_{h_chain}")
    if parsed_h is None:
        print(f"[{target}] WARNING: heavy chain {h_chain} not parsed from {pdb_path}; all-native.",
              flush=True)
        all_native_pool("heavy_chain_not_parsed")
        return

    # numbering parallel to parsed_h -------------------------------------------------------------
    numbering = chain_numbering(pdb_path, h_chain)
    if len(numbering) != len(parsed_h):
        print(f"[{target}] WARNING: numbering len {len(numbering)} != parsed seq len "
              f"{len(parsed_h)}; all-native.", flush=True)
        all_native_pool("numbering_length_mismatch")
        return
    key2idx = {}
    for idx, (rn, ic) in enumerate(numbering):
        if rn is not None:
            key2idx[(rn, (ic or "").strip().upper())] = idx

    # map template CDR residues -> heavy-chain flattened index, guard structural WT == native_cdr
    cdr_idx_in_h, struct_wt = [], []
    n_unmapped, n_wtmismatch = 0, 0
    for i, cr in enumerate(cdr_res):
        rn = cr.get("resnum")
        nat = native_cdr[i]
        if rn is None:
            cdr_idx_in_h.append(None); struct_wt.append("-"); n_unmapped += 1; continue
        key = (int(rn), (cr.get("icode") or "").strip().upper())
        idx = key2idx.get(key)
        if idx is None:
            cdr_idx_in_h.append(None); struct_wt.append("?"); n_unmapped += 1; continue
        wt = parsed_h[idx]
        struct_wt.append(wt)
        if wt != nat:
            cdr_idx_in_h.append(None); n_wtmismatch += 1; continue
        cdr_idx_in_h.append(idx)

    sample_pos = [(i, idx) for i, idx in enumerate(cdr_idx_in_h) if idx is not None]  # (cdr_i, h_idx)
    align_ok = (n_unmapped == 0 and n_wtmismatch == 0)
    print(f"[{target}] native_cdr(json)   = {native_cdr}", flush=True)
    print(f"[{target}] native_cdr(struct) = {''.join(struct_wt)}  align_ok={align_ok} "
          f"redesign {len(sample_pos)}/{len(cdr_res)} positions "
          f"(unmapped={n_unmapped}, wt_mismatch={n_wtmismatch})", flush=True)
    if not sample_pos:
        print(f"[{target}] WARNING: no CDR position registerable; emitting all-native pool.",
              flush=True)
        all_native_pool("registration_failed_all_native")
        return

    # --- ProteinMPNN featurize: masked=[H], visible=others; fix all H except CDR-H3 -------------
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[{target}] loading ProteinMPNN {os.path.basename(args.model_path)} on {device} ...",
          flush=True)
    ckpt = torch.load(args.model_path, map_location=device)
    model = ProteinMPNN(num_letters=21, node_features=128, edge_features=128, hidden_dim=128,
                        num_encoder_layers=3, num_decoder_layers=3, augment_eps=0.0,
                        k_neighbors=ckpt["num_edges"])
    model.to(device)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()

    chain_id_dict = {name: ([h_chain], other_chains)}   # (designable, fixed/visible)
    h_len = len(parsed_h)
    design_h_idx = set(idx for _, idx in sample_pos)
    fixed_pos_list = [j + 1 for j in range(h_len) if j not in design_h_idx]  # 1-based, fix all but CDR
    fixed_positions_dict = {name: {h_chain: fixed_pos_list}}
    for c in other_chains:
        fixed_positions_dict[name][c] = []  # visible anyway, but be explicit

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    (X, S, mask, lengths, chain_M, chain_encoding_all, _cll, _vll, _mll, _mcl, chain_M_pos,
     omit_AA_mask, residue_idx, _dm, _tpl, pssm_coef, pssm_bias, pssm_log_odds_all,
     bias_by_res_all, _tb) = tied_featurize(
        [pdb], device, chain_id_dict, fixed_positions_dict, None, None, None, None, ca_only=False)

    omit_AAs_np = np.array([aa in "X" for aa in MPNN_ALPHABET]).astype(np.float32)  # never emit X
    bias_AAs_np = np.zeros(len(MPNN_ALPHABET), dtype=np.float32)

    # global flattened indices of the CDR positions. H is the single masked chain => packed first,
    # so its heavy-chain index equals the global index. Verify against S (native tokens).
    global_cdr = {i: h_idx for i, h_idx in sample_pos}  # h_idx == global idx (H first)
    S_np = S[0].cpu().numpy()
    for i, gidx in global_cdr.items():
        tok = MPNN_ALPHABET[int(S_np[gidx])]
        if tok != native_cdr[i]:
            print(f"[{target}] WARNING: global index check failed at cdr pos {i}: "
                  f"S='{tok}' native='{native_cdr[i]}'; all-native fallback.", flush=True)
            all_native_pool("global_index_mismatch")
            return

    def sample_once(temp, seed_off):
        torch.manual_seed(args.seed + seed_off)
        randn = torch.randn(chain_M.shape, device=device)
        with torch.no_grad():
            sd = model.sample(X, randn, S, chain_M, chain_encoding_all, residue_idx, mask=mask,
                              temperature=max(temp, 1e-5), omit_AAs_np=omit_AAs_np,
                              bias_AAs_np=bias_AAs_np, chain_M_pos=chain_M_pos,
                              omit_AA_mask=omit_AA_mask, pssm_coef=pssm_coef, pssm_bias=pssm_bias,
                              pssm_multi=0.0, pssm_log_odds_flag=False,
                              pssm_log_odds_mask=None, pssm_bias_flag=False, bias_by_res=bias_by_res_all)
        return sd["S"][0].cpu().numpy()

    def assemble(S_sample):
        chars = list(native_cdr)
        for i, gidx in global_cdr.items():
            chars[i] = MPNN_ALPHABET[int(S_sample[gidx])]
        return "".join(chars)

    # Bounded dedup: try a modest number of temperature samples for unique loops, then pad with
    # (possibly duplicate) samples. The cap is small on purpose — antibody-finetuned weights (AbMPNN)
    # are near-deterministic on short/easy CDR-H3 windows, so an unbounded search would spin
    # thousands of full autoregressive model.sample() calls; duplicate designs in a matched pool are
    # already normal in this benchmark (template cluster_top1/top2 can be identical).
    max_unique_tries = 60
    seqs, seen = [], set()
    seqs.append(assemble(sample_once(1e-5, 0)))        # design 0: near-argmax (MPNN-optimal)
    seen.add(seqs[0])
    off = 1
    while len(seqs) < args.n_designs and off <= max_unique_tries:
        s = assemble(sample_once(args.sampling_temp, off)); off += 1
        if s not in seen:
            seen.add(s); seqs.append(s)
    while len(seqs) < args.n_designs:  # diversity exhausted (short/near-deterministic CDR): pad
        seqs.append(assemble(sample_once(args.sampling_temp, off))); off += 1

    designs = [{"name": nm, "cdr_seq": sq, "mpnn_mean": None, "mpnn_worst": None,
                "wt_ddg": None, "baddg_grade": None} for nm, sq in zip(names, seqs)]
    write_out(designs, align_ok=align_ok)

    print(f"\n[{target}] native  {native_cdr}")
    for d in designs:
        nmis = sum(a != b for a, b in zip(d["cdr_seq"], native_cdr))
        print(f"  {d['name']:<13} {d['cdr_seq']}  ({nmis} muts vs native)")
    print(f"\nwrote {len(designs)} ProteinMPNN designs x {len(tmpl['antigen_variants'])} "
          f"antigen variants -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
