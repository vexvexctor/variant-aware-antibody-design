#!/usr/bin/env python3
"""Edit-distance-matched RANDOM CDR-H3 control for the generator zoo.

The decisive sanity check on the grader. If random CDR-H3 edits, matched design-for-design to
the edit distance of a real arm, beat the native antibody at a similar rate, then that grader
cannot support ANY beat-native claim and the arm-vs-arm table means nothing.

One control arm is emitted per real arm so the edit distances are matched per design, per
target -- a control matched to AbMPNN's edit distance is not a valid control for a diffusion
arm that mutates twice as many positions.
"""
from __future__ import annotations
import argparse, glob, hashlib, json, os, random
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

AA = "ACDEFGHIKLMNPQRSTVWY"
NOCYS = [a for a in AA if a != "C"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="+", required=True)
    ap.add_argument("--designs-dir", default=f"{VAAD_ROOT}/scratch/genzoo/designs")
    ap.add_argument("--seed", type=int, default=20260819)
    a = ap.parse_args()

    made = 0
    for arm in a.arms:
        for src in sorted(glob.glob(f"{a.designs_dir}/*_designs_{arm}_r*.json")):
            if src.endswith(".cost.json"):
                continue
            base = os.path.basename(src)
            target, rest = base.split("_designs_", 1)
            rep = rest.rsplit("_r", 1)[1][:-5]
            out = f"{a.designs_dir}/{target}_designs_rand_{arm}_r{rep}.json"
            if os.path.exists(out):
                continue
            d = json.load(open(src))
            native = d["native_cdr"]
            # deterministic per (target, arm, rep) stream so the control is reproducible
            h = hashlib.sha256(f"{a.seed}:{target}:{arm}:{rep}".encode()).hexdigest()
            rng = random.Random(int(h[:16], 16))
            designs = []
            for x in d.get("antibody_designs", []):
                seq = x.get("cdr_seq") or native
                nmut = sum(1 for p, q in zip(native, seq) if p != q)
                chars = list(native)
                if nmut:
                    pos = rng.sample(range(len(native)), min(nmut, len(native)))
                    for i in pos:
                        # genuine substitution; never introduce an unpaired Cys
                        choices = [c for c in NOCYS if c != native[i]]
                        chars[i] = rng.choice(choices)
                designs.append({"name": x["name"], "cdr_seq": "".join(chars),
                                "mpnn_mean": None, "mpnn_worst": None,
                                "wt_ddg": None, "baddg_grade": None})
            o = dict(d)
            o.update(antibody_designs=designs, oracle=f"random_control_{arm}",
                     generator="random_edit_matched", mechanism="control", base="random",
                     matched_arm=arm, seed=a.seed,
                     note=("edit-distance-matched random CDR-H3; no oracle, no structure, no "
                           "escape variants -- pure noise at the same distance from native"))
            json.dump(o, open(out, "w"), indent=2)
            made += 1
    print(f"wrote {made} random-control design JSONs")


if __name__ == "__main__":
    main()
