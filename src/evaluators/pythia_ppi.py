#!/usr/bin/env python3
"""Pythia out-of-family arbiter on the third-party-grader battery (51 + 23 panels).

Thin wrapper: ALL scientific logic imported unchanged from score_pythia_escape.py.
Only per-arm paths (design-JSON template, source-PDB dir, output roots) are
re-parameterized via tp_config. Metric = worst-case REGRET vs native (native == 0 by
construction), NOT the FoldX/MM-GBSA absolute metric. CPU-only. Resumable.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, f"{VAAD_ROOT}/scratch/svdd_screen")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import torch  # noqa: E402
from Bio.PDB.Polypeptide import index_to_one  # noqa: E402

from score_pythia_escape import (  # noqa: E402  reuse verbatim
    cdr_rows,
    variant_relabel,
    write_subset,
)
from score_pythia import (  # noqa: E402
    AA_TO_IDX,
    PYTHIA_ROOT,
    build_row_index,
    energies,
)
from score_foldx import renumber_pdb  # noqa: E402
from pythia.masked_ddg_scan import get_torch_model  # noqa: E402

import tp_config as C  # noqa: E402
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)


def run(arm: str, target: str, models, device: str, limit_variants: int = 0):
    jpath = C.json_path(arm, target)
    outdir = os.path.join(C.PYTHIA_OUT, arm, target)
    workdir = os.path.join(C.OUT_ROOT, "work_pythia", arm, target)
    os.makedirs(outdir, exist_ok=True)
    ckpt_path = os.path.join(outdir, "E_checkpoint.json")
    meta_path = os.path.join(outdir, "status.json")

    if not os.path.exists(jpath):
        json.dump({"status": "no_designs", "path": jpath}, open(meta_path, "w"))
        print(f"[{arm}/{target}] no design JSON {jpath}; skip")
        return
    js = json.load(open(jpath))

    ab, ag = list(js["antibody_chains"]), list(js["antigen_chains"])
    src = C.pdb_path(arm, target)
    if not os.path.exists(src):
        json.dump({"status": "no_pdb", "path": src}, open(meta_path, "w"))
        print(f"[{arm}/{target}] missing source PDB {src}")
        return

    os.makedirs(workdir, exist_ok=True)
    renum = os.path.join(workdir, "complex_renum.pdb")
    cmap = renumber_pdb(src, renum)

    apo_pdb = write_subset(renum, os.path.join(workdir, "apo.pdb"), ab)
    e_apo, protbb_apo = energies(apo_pdb, models, device)
    idx_apo = build_row_index(apo_pdb)
    rows_apo, err = cdr_rows(js, cmap, idx_apo, protbb_apo)
    if err:
        json.dump({"status": f"apo_cdr:{err}"}, open(meta_path, "w"))
        print(f"[{arm}/{target}] apo CDR mapping failed: {err}")
        return

    variants = [v for v in js["antigen_variants"] if v.get("split") == "eval"]
    if limit_variants:
        variants = variants[:limit_variants]
    designs = list(js["antibody_designs"])
    print(f"[{arm}/{target}] device={device} {len(designs)} designs x {len(variants)} "
          f"eval variants; {len(variants) + 1} forward passes", flush=True)

    E = {"WT": {}}
    errors = {}
    for v in variants:
        vname = v["name"]
        rel, verr = variant_relabel(v, cmap)
        if verr:
            errors[vname] = verr
            continue
        vpdb = write_subset(renum, os.path.join(workdir, "variant.pdb"), ab + ag, rel)
        try:
            e_c, protbb_c = energies(vpdb, models, device)
            idx_c = build_row_index(vpdb)
        except Exception as exc:  # noqa: BLE001
            errors[vname] = f"{type(exc).__name__}:{exc}"
            continue
        rows_c, cerr = cdr_rows(js, cmap, idx_c, protbb_c)
        if cerr:
            errors[vname] = cerr
            continue

        E["WT"][vname] = 0.0
        for d in designs:
            cdr = d["cdr_seq"]
            if len(cdr) != len(rows_c):
                errors[f"{d['name']}|{vname}"] = "cdr_len_mismatch"
                continue
            total, bad = 0.0, False
            for i, aa in enumerate(cdr):
                ai = AA_TO_IDX.get(aa)
                if ai is None:
                    bad = True
                    break
                total += float(e_c[rows_c[i]][ai]) - float(e_apo[rows_apo[i]][ai])
            if bad:
                errors[f"{d['name']}|{vname}"] = "nonstd_aa"
                continue
            E.setdefault(d["name"], {})[vname] = total
        json.dump(E, open(ckpt_path, "w"))

    json.dump(E, open(ckpt_path, "w"))
    done = sum(len(x) for x in E.values())
    json.dump({"status": "ok", "scored": done, "errors": errors,
               "n_designs": len(designs), "n_variants": len(variants),
               "metric": "worst_case_regret_vs_native (NOT the FoldX absolute metric)"},
              open(meta_path, "w"), indent=1)
    print(f"[{arm}/{target}] {done} values, {len(errors)} errors -> {ckpt_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, default=None)
    ap.add_argument("--arm", default=None)
    ap.add_argument("--target", default=None)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--limit-variants", type=int, default=0)
    args = ap.parse_args()

    jobs = C.jobs()
    if args.arm and args.target:
        arm, target = args.arm, args.target
    elif args.index is not None:
        if args.index >= len(jobs):
            print(f"index {args.index} >= njobs {len(jobs)}")
            return
        arm, target = jobs[args.index]
    else:
        raise SystemExit("need --index or (--arm and --target)")

    models = [
        get_torch_model(os.path.join(PYTHIA_ROOT, "pythia", "pythia-c.pt"), args.device),
        get_torch_model(os.path.join(PYTHIA_ROOT, "pythia", "pythia-p.pt"), args.device),
    ]
    run(arm, target, models, args.device, args.limit_variants)


if __name__ == "__main__":
    main()
