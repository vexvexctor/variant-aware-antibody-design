#!/usr/bin/env python3
"""H3-DDG GRADING pass over the F5 antibody-design validation set.

H3-DDG is a GRADER ONLY — it never steered these designs (cf. BA-DDG, which did and is
therefore circular). It is the NeurIPS-2025 hypergraph many-body ΔΔG model
(github.com/sodaball/H3-DDG), fine-tuned FROM BA-DDG's soluble ProteinMPNN backbone
(v_48_020) on SKEMPI v2 with added Triplet/Quadruplet hypergraph layers. Architecturally
derived from BA-DDG (same ProteinMPNN backbone + thermodynamic-cycle ΔΔG head), so it is a
CORROBORATING oracle, LESS orthogonal than FoldX (which is the independent physics verdict).

This mirrors design_for_validation.py:_baddg_grade / baddg_score_designs.py, but loads the
H3-DDG checkpoint into H3-DDG's own (hypergraph) DDGPredictor and uses the model's trained
thermodynamic-cycle head output (out['ddG_pred']) — H3-DDG's headline signal.

For each (antibody_design x antigen_variant) it computes the design-vs-native binding ΔΔG on
that variant's antigen context:
  * parse the WT complex once (fixed backbone, like FoldX/BA-DDG),
  * apply the antigen variant's point mutations to BOTH the wt (aa) and mut (aa_mut) sequences
    (so the variant antigen is the shared structural context for both states),
  * apply the design's CDR substitutions to aa_mut only (native CDR stays in aa),
  * mut_flag = (aa != aa_mut) marks only the CDR positions that differ -> the thermodynamic
    cycle isolates the antibody-binding effect of native->design CDR on that variant.
ddG = mut_scores_cycle - wt_scores_cycle (NLL units). Convention: MORE NEGATIVE = design binds
TIGHTER than native on that variant (same sign as baddg_grade; <0 = win).

Folds: the checkpoint is 3-fold CV (full3fold). These F5 targets are out-of-distribution vs
SKEMPI, so by default we ENSEMBLE the 3 folds (mean ddG_pred) — the standard inference mode
for a k-fold ΔΔG model on unseen complexes. --fold N restricts to a single fold.

Output: one JSON per (target,label) in --outdir, schema
  {"target","label","steer_agg","wt_weight","native_cdr","folds",
   "h3ddg_grade": {design_name: {variant_name: ddg, ...}, ...},
   "errors": {design_name: {variant_name: msg}}}
mirroring antibody_designs[].baddg_grade so f5_aggregate_h3ddg.py can slice by split and
compute the same worst-case / best-of-pool / beat-native matrix.

Runs under system python3 (torch + biopython + easydict). GPU if available.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import torch

H3 = f"{VAAD_ROOT}/tools/H3-DDG"
sys.path.insert(0, H3)

from Bio.PDB.PDBParser import PDBParser                      # noqa: E402
from Bio.PDB.Polypeptide import one_to_index, index_to_one  # noqa: E402
from easydict import EasyDict                                # noqa: E402

from common_utils.protein.parsers import parse_biopython_structure  # noqa: E402
from common_utils.transforms import get_transform                   # noqa: E402
from dataset import MPNNPaddingCollate                              # noqa: E402
from ddg_predictor import DDGPredictor                             # noqa: E402
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

# ProteinMPNN / Bio 1-letter alphabet: A..Y = 0..19 (matches Bio's index_to_one), X = 20
# (the model's mask/unknown token). The diffusion may emit 'X' at undetermined CDR positions;
# the in-objective baddg_grade scored those as token 20, so we map them the same way.
_ALPHA = "ACDEFGHIKLMNPQRSTVWYX"
_C2I = {c: i for i, c in enumerate(_ALPHA)}


def _aa_idx(c):
    return _C2I.get(c, 20)  # unknown/non-standard -> mask token 20


CKPT = f"{VAAD_ROOT}/results/h3ddg_train/2026-06-26-03-47-43-752_full3fold/checkpoint/ddg_model.ckpt"
# matches train_config.json (architecture H3-DDG was trained with)
CFG = EasyDict({
    "ca_only": False, "hidden_dim": 128, "num_layers": 3, "backbone_noise": 0.0,
    "num_edges": 48, "loss_weight_boltzmann": 1.0, "num_tri_heads": 4,
    "use_hypergraph": True, "max_num_hyperedges": 420, "num_mut_subgraph_nodes": 20,
    "hyper_ratio": 4, "num_edges_ratio": 3.0, "edges_selection": "dynamic",
})


def load_models(device, folds):
    sd = torch.load(CKPT, map_location="cpu", weights_only=False)
    models = []
    for fi in folds:
        m = DDGPredictor(CFG).to(device)
        missing, unexpected = m.load_state_dict(sd["models"][fi], strict=False)
        if missing or unexpected:
            raise RuntimeError(f"fold {fi} state_dict mismatch: "
                               f"missing={list(missing)[:4]} unexpected={list(unexpected)[:4]}")
        m.eval()
        models.append(m)
    return models


def build_base_data(pdb_path, ab_chains, ag_chains):
    """Parse the complex once; return (data, idx_map) with idx_map[(chain,resnum,icode)]=row."""
    parser = PDBParser(QUIET=True)
    model = parser.get_structure(None, pdb_path)[0]
    data, _ = parse_biopython_structure(
        model, antibody_chain_id=list(ab_chains), antigen_chain_id=list(ag_chains))
    if data is None:
        raise RuntimeError(f"parse failed for {pdb_path}")
    data = get_transform([{"type": "select_atom", "resolution": "backbone+CB"}])(data)
    idx_map = {}
    for i, (c, rs, ic) in enumerate(zip(data["chain_id"], data["resseq"].tolist(), data["icode"])):
        idx_map[(c, int(rs), (ic or "").strip())] = i
    data["ddG"] = 0.0
    data["id"] = 0
    data["complex"] = "design"
    data["num_muts"] = 0
    data["mutstr"] = ""
    return data, idx_map


@torch.no_grad()
def score_pair(models, base_data, cdr_idx, native_cdr, design_cdr, idx_map,
               ag_mutations, collate, device):
    """Ensembled H3-DDG design-vs-native binding ΔΔG on this antigen variant.
    Returns float ddG, or {'_error': msg}. The antigen variant is applied to BOTH aa and
    aa_mut so it is the shared context; only the CDR (native->design) differs -> mut_flag."""
    data = dict(base_data)
    aa = base_data["aa"].clone()
    aa_mut = base_data["aa"].clone()

    # antigen point mutations -> applied to BOTH states (shared variant context)
    for m in ag_mutations:
        key = (m["chain"], int(m["resnum"]), (m.get("icode") or "").strip())
        res = idx_map.get(key)
        if res is None:
            return {"_error": f"ag {key} unmapped"}
        if index_to_one(int(aa[res])) != m["wt"]:
            return {"_error": f"ag {key} parsed={index_to_one(int(aa[res]))} json_wt={m['wt']}"}
        aa[res] = _aa_idx(m["mut"])
        aa_mut[res] = _aa_idx(m["mut"])

    # CDR substitutions -> design on aa_mut only (native stays on aa)
    for i, (res, n1, d1) in enumerate(zip(cdr_idx, native_cdr, design_cdr)):
        if res is None:
            return {"_error": f"cdr_pos{i} unmapped"}
        if index_to_one(int(base_data["aa"][res])) != n1:
            return {"_error": f"cdr_pos{i} parsed={index_to_one(int(base_data['aa'][res]))} json_native={n1}"}
        if d1 != n1:
            aa_mut[res] = _aa_idx(d1)

    data["aa"] = aa
    data["aa_mut"] = aa_mut
    data["mut_flag"] = (aa != aa_mut)
    if int(data["mut_flag"].sum()) == 0:
        return 0.0  # native CDR (e.g. WT antibody) -> ddG 0 by definition

    batch = collate([data])
    batch = {k: (v.to(device) if torch.is_tensor(v) else v) for k, v in batch.items()}
    vals = []
    for model in models:
        _, out, _ = model(batch)
        vals.append(float(out["ddG_pred"][0].item()))
    return sum(vals) / len(vals)


def grade_spec(spec, pdb_path, models, collate, device):
    ab_chains, ag_chains = spec["antibody_chains"], spec["antigen_chains"]
    native_cdr = spec["native_cdr"]
    cdr_res = spec["cdr_residues"]
    base_data, idx_map = build_base_data(pdb_path, ab_chains, ag_chains)
    cdr_idx = [idx_map.get((r["chain"], int(r["resnum"]), (r.get("icode") or "").strip()))
               for r in cdr_res]
    grades, errors = {}, {}
    for ab in spec["antibody_designs"]:
        g, e = {}, {}
        for ag in spec["antigen_variants"]:
            r = score_pair(models, base_data, cdr_idx, native_cdr, ab["cdr_seq"], idx_map,
                           ag.get("mutations", []), collate, device)
            if isinstance(r, dict):
                e[ag["name"]] = r["_error"]
            else:
                g[ag["name"]] = round(r, 4)
        grades[ab["name"]] = g
        if e:
            errors[ab["name"]] = e
    return grades, errors


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default=f"{VAAD_ROOT}/scratch/steer_tests/f5_rerun_manifest.tsv")
    ap.add_argument("--valdir", default=f"{VAAD_ROOT}/results/validation")
    ap.add_argument("--pdbdir", default=f"{VAAD_ROOT}/datasets/AACDB/data_zip/complex_structure")
    ap.add_argument("--outdir", default=f"{VAAD_ROOT}/scratch/steer_tests/h3ddg_out")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--fold", type=int, default=-1, help="single fold; -1 = ensemble all 3")
    ap.add_argument("--single", default=None, help="run one design JSON path (smoke test) -> stdout")
    args = ap.parse_args()

    device = torch.device(args.device)
    folds = range(3) if args.fold < 0 else [args.fold]
    print(f"loading H3-DDG folds {list(folds)} on {device} ...", flush=True)
    models = load_models(device, folds)
    collate = MPNNPaddingCollate()

    if args.single:
        spec = json.load(open(args.single))
        pdb = os.path.join(args.pdbdir, spec.get("pdb", spec["target"] + ".pdb"))
        g, e = grade_spec(spec, pdb, models, collate, device)
        print(json.dumps({"grades": g, "errors": e}, indent=2))
        return

    os.makedirs(args.outdir, exist_ok=True)
    rows = [ln.split("\t") for ln in open(args.manifest) if ln.strip()]
    print(f"{len(rows)} target x condition pairs", flush=True)
    for i, (label, target) in enumerate(rows):
        label, target = label.strip(), target.strip()
        out = os.path.join(args.outdir, f"{target}__{label}.json")
        if os.path.exists(out):
            print(f"[{i+1}/{len(rows)}] {target} {label} -> cached", flush=True)
            continue
        jpath = os.path.join(args.valdir, f"{target}_designs_{label}.json")
        if not os.path.exists(jpath):
            print(f"[{i+1}/{len(rows)}] MISSING design JSON {jpath}", flush=True)
            continue
        spec = json.load(open(jpath))
        pdb = os.path.join(args.pdbdir, spec.get("pdb", target + ".pdb"))
        t0 = time.time()
        try:
            grades, errors = grade_spec(spec, pdb, models, collate, device)
        except Exception as ex:  # noqa: BLE001
            print(f"[{i+1}/{len(rows)}] {target} {label} FAILED: {type(ex).__name__}: {ex}", flush=True)
            continue
        res = {"target": target, "label": label,
               "steer_agg": spec.get("steer_agg"), "wt_weight": spec.get("wt_weight"),
               "native_cdr": spec["native_cdr"], "folds": list(folds),
               "h3ddg_grade": grades, "errors": errors}
        tmp = out + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(res, fh)
        os.replace(tmp, out)
        nerr = sum(len(v) for v in errors.values())
        print(f"[{i+1}/{len(rows)}] {target} {label} done in {time.time()-t0:.1f}s "
              f"({len(grades)} designs, {nerr} errors)", flush=True)


if __name__ == "__main__":
    main()
