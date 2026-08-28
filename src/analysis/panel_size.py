#!/usr/bin/env python3
"""Design-panel-size sensitivity: beat-native and cost vs M_design in {1,3,6,12}.

Every arm is graded on the SAME fixed 24 held-out variants (verified byte-identical:
the eval panel is the bottom-24 by omega, the steer panel the top-M, and every panel
target has >= 41 variants so the two never intersect for M <= 12).

Grader = H3-DDG, out-of-family relative to the BA-DDG steering objective.
Margin = design-minus-native; < 0 = design binds that variant tighter than native.
beat-native = best-of-pool worst-case margin < 0.
CI = 95% bootstrap over ANTIGEN CLUSTERS (30 distinct), not targets.
"""
import csv, json, os, random, statistics as st
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

SCR = f"{VAAD_ROOT}/scratch/panel_size"
VAL, H3 = f"{SCR}/validation", f"{SCR}/h3_out"
MS, REPS = [1, 3, 6, 12], [1, 2, 3]
TARGETS = [l.strip() for l in open(f"{SCR}/panel_targets.txt") if l.strip()]
CLU = {r["target"]: r["ag_cluster"] for r in csv.DictReader(
    open(f"{VAAD_ROOT}/antibody-robustness/results/deliverables_2026-08-19"
         "/01_vs_native/per_target_margins.csv"))}


def run_units(target, M, rep):
    """-> (per_design_worst list, cost dict) or None if incomplete."""
    dp = f"{VAL}/{target}_designs_m{M}_r{rep}.json"
    gp = f"{H3}/{target}__m{M}_r{rep}.json"
    if not (os.path.exists(dp) and os.path.exists(gp)):
        return None
    d = json.load(open(dp))
    ev = [v["name"] for v in d["antigen_variants"] if v.get("split") == "eval"]
    grade = json.load(open(gp)).get("h3ddg_grade") or {}
    worsts = []
    for name, per_var in grade.items():
        vals = [per_var[v] for v in ev if v in per_var]
        if len(vals) == len(ev) and vals:
            worsts.append(max(vals))          # worst held-out variant for this design
    if not worsts:
        return None
    cost = {}
    cp = dp[:-5] + ".cost.json"
    if os.path.exists(cp):
        c = json.load(open(cp))
        cost = {"oracle_calls": c["counters"].get("oracle_calls", 0),
                "reward_evals": c["counters"].get("reward_evals", 0),
                "wall_s": c.get("wall_s", 0), "cuda_s": c.get("cuda_s", 0)}
    return worsts, cost, len(ev)


def boot_ci(vals, clusters, n=4000, seed=0):
    """95% cluster bootstrap CI on the mean of `vals`."""
    by = {}
    for v, c in zip(vals, clusters):
        by.setdefault(c, []).append(v)
    keys, rng, out = list(by), random.Random(seed), []
    if not keys:
        return (float("nan"), float("nan"))
    for _ in range(n):
        s = [x for k in (rng.choice(keys) for _ in keys) for x in by[k]]
        out.append(sum(s) / len(s))
    out.sort()
    return out[int(.025 * n)], out[int(.975 * n)]


rows, per_unit = [], []
for M in MS:
    pooled_hit, pooled_clu, pooled_marg = [], [], []
    run_hit, run_clu, run_marg = [], [], []
    oc, we, cu, n_ev_seen = [], [], [], set()
    for t in TARGETS:
        pool = []
        for r in REPS:
            got = run_units(t, M, r)
            if not got:
                continue
            worsts, cost, n_ev = got
            n_ev_seen.add(n_ev)
            pool += worsts
            bp_run = min(worsts)
            run_hit.append(1.0 if bp_run < 0 else 0.0)
            run_marg.append(bp_run); run_clu.append(CLU.get(t, t))
            per_unit.append({"target": t, "ag_cluster": CLU.get(t, t), "M_design": M, "rep": r,
                             "n_designs": len(worsts), "best_of_pool_worst_margin": round(bp_run, 4),
                             "beat_native": int(bp_run < 0), **cost})
            if cost:
                oc.append(cost["oracle_calls"]); we.append(cost["reward_evals"]); cu.append(cost["cuda_s"])
        if pool:
            bp = min(pool)
            pooled_hit.append(1.0 if bp < 0 else 0.0)
            pooled_marg.append(bp); pooled_clu.append(CLU.get(t, t))
    if not pooled_hit:
        print(f"M={M}: no complete units yet"); continue
    lo, hi = boot_ci(pooled_hit, pooled_clu)
    rlo, rhi = boot_ci(run_hit, run_clu)
    rows.append({
        "M_design": M, "n_targets": len(pooled_hit), "n_runs": len(run_hit),
        "n_clusters": len(set(pooled_clu)), "n_eval_variants": sorted(n_ev_seen),
        "beat_native_pooled": 100 * st.mean(pooled_hit), "ci_lo": 100 * lo, "ci_hi": 100 * hi,
        "beat_native_per_run": 100 * st.mean(run_hit), "run_lo": 100 * rlo, "run_hi": 100 * rhi,
        "median_margin": st.median(pooled_marg),
        "oracle_calls": st.mean(oc) if oc else 0, "reward_evals": st.mean(we) if we else 0,
        "cuda_s": st.mean(cu) if cu else 0,
    })

with open(f"{SCR}/panel_size_per_unit.csv", "w", newline="") as fh:
    if per_unit:
        w = csv.DictWriter(fh, fieldnames=list(per_unit[0])); w.writeheader(); w.writerows(per_unit)
with open(f"{SCR}/panel_size_summary.csv", "w", newline="") as fh:
    if rows:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

print(f"{'M':>3} {'beat-native (pooled)':>26} {'per-run':>20} {'medMargin':>10} "
      f"{'oracleCalls':>12} {'rewardEvals':>12} {'cuda_s':>8} {'nT':>4} {'nRun':>5} {'nClu':>5} {'nEval':>7}")
for r in rows:
    print(f"{r['M_design']:>3} {r['beat_native_pooled']:>8.1f}% [{r['ci_lo']:>5.1f},{r['ci_hi']:>5.1f}]"
          f" {r['beat_native_per_run']:>8.1f}% [{r['run_lo']:>4.1f},{r['run_hi']:>4.1f}]"
          f" {r['median_margin']:>10.3f} {r['oracle_calls']:>12.0f} {r['reward_evals']:>12.0f}"
          f" {r['cuda_s']:>8.1f} {r['n_targets']:>4} {r['n_runs']:>5} {r['n_clusters']:>5}"
          f" {str(r['n_eval_variants']):>7}")
