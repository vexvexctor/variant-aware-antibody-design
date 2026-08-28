#!/usr/bin/env python3
"""MM-GBSA on the sp1 shared-pool frozen picks, sharded by design.

All scientific logic (force field, prep, single-trajectory dG_bind, sign convention) is
imported UNCHANGED from score_mmgbsa_escape.py, the same module the published 51-panel escape
run used. This file only:
  - points at the sp1 shared-pool design JSONs (20 distinct AntiFold candidates per target,
    the union of what the 7 frozen selectors picked across 5 seeds), and
  - shards the design list, because one target is 20 endpoints x 24 eval variants = ~500
    dG_bind evaluations, far too much for a single array task.

The native antibody ("WT" endpoint) is graded in shard 0 only; margins need it exactly once.
Each shard writes its own checkpoint so shards never race on one file
(see the campaign lesson about concurrent jobs clobbering a shared target list).
"""
from __future__ import annotations
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

import argparse
import json
import os
import sys
import traceback

sys.path.insert(0, f"{VAAD_ROOT}/scratch/svdd_screen")
sys.path.insert(0, f"{VAAD_ROOT}/scratch/thirdparty_grade")

from score_mmgbsa_escape import (  # noqa: E402
    build_and_score,
    cdr_mutations,
    variant_mutations,
)
from score_mmgbsa import (  # noqa: E402
    PLATFORM_NAME,
    raw_resname_map,
    renumber_pdb,
)

SEL_DIR = f"{VAAD_ROOT}/scratch/campaign/sp1/selected"
AACDB = f"{VAAD_ROOT}/datasets/AACDB/data_zip/complex_structure"
OUT = f"{VAAD_ROOT}/scratch/sp1grade/mmgbsa_out"
WORK = f"{VAAD_ROOT}/scratch/sp1grade/work_mmgbsa"


def targets():
    return sorted(f.split("_designs_")[0] for f in os.listdir(SEL_DIR)
                  if f.endswith("_designs_sp1all.json"))


def run(target: str, shard: int, nshard: int):
    jpath = os.path.join(SEL_DIR, f"{target}_designs_sp1all.json")
    js = json.load(open(jpath))
    ab, ag = list(js["antibody_chains"]), list(js["antigen_chains"])
    src = os.path.join(AACDB, f"{target}.pdb")
    if not os.path.exists(src):
        print(f"[{target}] missing PDB {src}")
        return

    outdir = os.path.join(OUT, target)
    workdir = os.path.join(WORK, target, f"s{shard}")
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(workdir, exist_ok=True)
    ckpt_path = os.path.join(outdir, f"E_shard{shard}.json")

    renum = os.path.join(workdir, "complex_renum.pdb")
    cmap = renumber_pdb(src, renum, ab + ag)
    raw_names = raw_resname_map(renum)

    variants = [v for v in js["antigen_variants"] if v.get("split") == "eval"]
    designs = list(js["antibody_designs"])
    mine = designs[shard::nshard]                      # stride, so shards are balanced
    endpoints = [(d["name"], d["cdr_seq"]) for d in mine]
    if shard == 0:
        endpoints = [("WT", js["native_cdr"])] + endpoints

    print(f"[{target} shard {shard}/{nshard}] platform={PLATFORM_NAME} "
          f"{len(endpoints)} endpoints x {len(variants)} variants "
          f"= {len(endpoints) * len(variants)} dG_bind", flush=True)

    E = json.load(open(ckpt_path)) if os.path.exists(ckpt_path) else {}
    errors = {}
    for name, cdr_seq in endpoints:
        E.setdefault(name, {})
        cdr_muts, cdr_chain, err = cdr_mutations(js, cdr_seq, cmap, raw_names)
        if err:
            errors[name] = err
            print(f"  [{name}] SKIP {err}", flush=True)
            continue
        for v in variants:
            if v["name"] in E[name]:
                continue
            ag_by_chain, verr = variant_mutations(v, cmap, raw_names)
            if verr:
                errors[f"{name}|{v['name']}"] = verr
                continue
            try:
                E[name][v["name"]] = build_and_score(
                    renum, cdr_muts, cdr_chain, ag_by_chain, ab, ag)
            except Exception as exc:  # noqa: BLE001
                errors[f"{name}|{v['name']}"] = f"{type(exc).__name__}:{exc}"
                traceback.print_exc()
                continue
            json.dump(E, open(ckpt_path, "w"))
        print(f"  [{name}] {len(E[name])}/{len(variants)} scored", flush=True)

    json.dump(E, open(ckpt_path, "w"))
    done = sum(len(v) for v in E.values())
    json.dump({"status": "ok", "scored": done, "errors": errors,
               "n_endpoints": len(endpoints), "n_variants": len(variants)},
              open(os.path.join(outdir, f"status_shard{shard}.json"), "w"), indent=1)
    print(f"[{target} shard {shard}] {done} values, {len(errors)} errors -> {ckpt_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--nshard", type=int, default=5)
    args = ap.parse_args()
    tg = targets()
    tgt_i, shard = divmod(args.index, args.nshard)
    if tgt_i >= len(tg):
        print(f"index {args.index} beyond {len(tg)} targets x {args.nshard} shards")
        return
    run(tg[tgt_i], shard, args.nshard)


if __name__ == "__main__":
    main()
