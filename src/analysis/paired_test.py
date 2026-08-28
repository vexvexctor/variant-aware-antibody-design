#!/usr/bin/env python3
"""Paired M-vs-M comparison. Same targets, same fixed held-out panel, so pair per target
and bootstrap the paired difference over antigen clusters. Far more powerful than
comparing overlapping marginal CIs."""
import csv, json, os, random, statistics as st
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)
SCR=f"{VAAD_ROOT}/scratch/panel_size"
MS,REPS=[1,3,6,12],[1,2,3]
T=[l.strip() for l in open(f"{SCR}/panel_targets.txt") if l.strip()]
CLU={r["target"]:r["ag_cluster"] for r in csv.DictReader(open(
 f"{VAAD_ROOT}/antibody-robustness/results/deliverables_2026-08-19/01_vs_native/per_target_margins.csv"))}

def bp(t,M):
    """best-of-pool worst-case margin over the fixed 24 held-out variants, reps pooled."""
    pool=[]
    for r in REPS:
        dp,gp=f"{SCR}/validation/{t}_designs_m{M}_r{r}.json", f"{SCR}/h3_out/{t}__m{M}_r{r}.json"
        if not(os.path.exists(dp) and os.path.exists(gp)): continue
        ev=[v["name"] for v in json.load(open(dp))["antigen_variants"] if v.get("split")=="eval"]
        for _,pv in (json.load(open(gp)).get("h3ddg_grade") or {}).items():
            vals=[pv[v] for v in ev if v in pv]
            if len(vals)==len(ev): pool.append(max(vals))
    return min(pool) if pool else None

M2={M:{t:bp(t,M) for t in T} for M in MS}

def paired(a,b,n=10000,seed=0):
    """mean(b-a) on margins and pp difference in beat-native, cluster-bootstrapped."""
    rows=[(CLU.get(t,t), M2[a][t], M2[b][t]) for t in T if M2[a][t] is not None and M2[b][t] is not None]
    by={}
    for c,x,y in rows: by.setdefault(c,[]).append((x,y))
    keys=list(by); rng=random.Random(seed); dm=[]; dp=[]
    for _ in range(n):
        s=[p for k in (rng.choice(keys) for _ in keys) for p in by[k]]
        dm.append(st.mean(y-x for x,y in s))
        dp.append(100*(st.mean(1 if y<0 else 0 for _,y in s)-st.mean(1 if x<0 else 0 for x,_ in s)))
    dm.sort(); dp.sort()
    obs_m=st.mean(y-x for _,x,y in rows)
    obs_p=100*(st.mean(1 if y<0 else 0 for _,_,y in rows)-st.mean(1 if x<0 else 0 for _,x,_ in rows))
    return len(rows),obs_m,(dm[250],dm[9750]),obs_p,(dp[250],dp[9750])

print("Paired vs M=6 (the current choice).  margin: negative = that M is TIGHTER/better.")
print(f"{'pair':>12} {'n':>3} {'Δmargin':>9} {'95% CI':>18} {'Δbeat-native pp':>16} {'95% CI':>18} {'sig?':>5}")
for a in [1,3,12]:
    n,dm,(l,h),dp_,(pl,ph)=paired(6,a)
    sig = "yes" if (l>0 or h<0) else "no"
    print(f"{'M=6 -> M='+str(a):>12} {n:>3} {dm:>+9.3f} [{l:>+7.3f},{h:>+7.3f}] {dp_:>+16.1f} [{pl:>+7.1f},{ph:>+7.1f}] {sig:>5}")
print()
print("Paired vs M=1 (the cheapest arm).")
for a in [3,6,12]:
    n,dm,(l,h),dp_,(pl,ph)=paired(1,a)
    sig = "yes" if (l>0 or h<0) else "no"
    print(f"{'M=1 -> M='+str(a):>12} {n:>3} {dm:>+9.3f} [{l:>+7.3f},{h:>+7.3f}] {dp_:>+16.1f} [{pl:>+7.1f},{ph:>+7.1f}] {sig:>5}")
