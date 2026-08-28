#!/usr/bin/env python3
"""Out-of-family beat-native from Pythia, mirroring the H3-DDG protocol exactly.

Pythia shares no lineage with ProteinMPNN, so unlike H3-DDG it is not in-family for the
AbMPNN/ProteinMPNN arms. Same statistic as H3: worst-case over the HELD-OUT eval variants,
best-of-pool, versus the native antibody (which Pythia reports as 0 by construction).
Random edit-matched controls are included so we can tell whether Pythia discriminates at all
before trusting any rate it produces.
"""
from __future__ import annotations
import csv, json, os, random, statistics as st
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

G = f"{VAAD_ROOT}/scratch/genzoo"
P = f"{VAAD_ROOT}/scratch/thirdparty_grade/pythia_out"
panel = {r["target"]: r for r in csv.DictReader(open(f"{G}/cand/panel_expanded.tsv"), delimiter="\t")}


def cell(arm, t):
    f = f"{P}/gz_{arm}/{t}/E_checkpoint.json"
    dj = f"{G}/designs/{t}_designs_{arm}_r1.json"
    if not (os.path.exists(f) and os.path.exists(dj)):
        return None
    d = json.load(open(dj))
    ev = {v["name"] for v in d.get("antigen_variants", []) if v.get("split") == "eval"}  # name IS the mutation string Pythia keys by
    try:
        g = json.load(open(f))
    except Exception:
        return None
    W = []
    for name, vals in g.items():
        if name == "WT":
            continue
        vv = [vals[k] for k in vals if k in ev]
        if len(vv) == len(ev) and vv:
            W.append(max(vv))
    return min(W) if W else None


def crate(d):
    by = {}
    for t, v in d.items():
        by.setdefault(panel[t]["ag_cluster"], []).append(v)
    return 100 * st.mean([st.mean(v) for v in by.values()]) if by else float("nan")


def boot(by, n=2000, seed=0):
    rng = random.Random(seed); cs = list(by); cm = {c: st.mean(by[c]) for c in cs}
    reps = sorted(st.mean([cm[rng.choice(cs)] for _ in cs]) for _ in range(n))
    return 100 * st.mean(cm.values()), 100 * reps[int(.025 * n)], 100 * reps[int(.975 * n)]


print("=== PYTHIA (OUT-OF-FAMILY) beat-native, held-out escape split ===\n")
arms = ["base_abmpnn_ar", "ours_abmpnn_live", "base_antifold", "ours_antifold",
        "base_proteinmpnn_ar", "ours_proteinmpnn_live",
        "rand_base_abmpnn_ar", "rand_ours_abmpnn_live", "base_rfantibody", "base_diffab"]
print(f"{'arm':26s} {'beat-native':>12s} {'n':>5s} {'mean margin':>12s}")
for a in arms:
    d = {t: cell(a, t) for t in panel}
    d = {t: v for t, v in d.items() if v is not None}
    if len(d) < 5:
        continue
    print(f"{a:26s} {crate({t:(1.0 if v<0 else 0.0) for t,v in d.items()}):11.1f}% "
          f"{len(d):5d} {st.mean(d.values()):12.3f}")

print("\n=== PAIRED: generator alone vs + our objective (Pythia) ===\n")
for name, b, o in [("AbMPNN", "base_abmpnn_ar", "ours_abmpnn_live"),
                   ("AntiFold", "base_antifold", "ours_antifold"),
                   ("ProteinMPNN", "base_proteinmpnn_ar", "ours_proteinmpnn_live")]:
    rows = [t for t in panel if cell(b, t) is not None and cell(o, t) is not None]
    if len(rows) < 10:
        print(f"{name}: n={len(rows)} too few yet"); continue
    B = {t: (1.0 if cell(b, t) < 0 else 0.0) for t in rows}
    O = {t: (1.0 if cell(o, t) < 0 else 0.0) for t in rows}
    by = {}
    for t in rows:
        by.setdefault(panel[t]["ag_cluster"], []).append(O[t] - B[t])
    dl, lo, hi = boot(by)
    sig = " **" if (lo > 0 or hi < 0) else ""
    print(f"{name:14s} base {crate(B):5.1f}%  ours {crate(O):5.1f}%  "
          f"delta {dl:+6.1f}{sig} [{lo:+.0f}, {hi:+.0f}]  n={len(rows)}")
