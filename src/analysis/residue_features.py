#!/usr/bin/env python3
"""What do escape-ROBUST CDR-H3 designs have in common?

Pools every graded design across arms and asks which sequence features track the held-out
worst-case escape margin. Edit distance from native is a known confound (the grader penalises
distance monotonically), so effects are reported BOTH raw and within edit-distance strata.
"""
from __future__ import annotations
import json, glob, os, statistics as st
from collections import defaultdict
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

G = f"{VAAD_ROOT}/scratch/genzoo"
KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,
      "I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,
      "Y":-1.3,"V":4.2}
CHG = {"D":-1,"E":-1,"K":1,"R":1,"H":0.1}
AROM, TINY, FLEX = set("FWY"), set("GAS"), set("GP")

def feats(s):
    n = len(s)
    return {"len": n,
            "charge": sum(CHG.get(c,0) for c in s),
            "hydro": st.mean(KD.get(c,0) for c in s),
            "arom": sum(c in AROM for c in s)/n,
            "tiny": sum(c in TINY for c in s)/n,
            "flex": sum(c in FLEX for c in s)/n,
            "Y": s.count("Y")/n, "W": s.count("W")/n, "D": s.count("D")/n,
            "G": s.count("G")/n, "S": s.count("S")/n, "R": s.count("R")/n}

def corr(xs, ys):
    n=len(xs); mx,my=st.mean(xs),st.mean(ys)
    num=sum((a-mx)*(b-my) for a,b in zip(xs,ys))
    den=(sum((a-mx)**2 for a in xs)*sum((b-my)**2 for b in ys))**0.5
    return num/den if den else float("nan")

rows=[]
for gf in glob.glob(f"{G}/h3_out/*.json"):
    b=os.path.basename(gf)[:-5]
    if "__" not in b: continue
    t,lab=b.split("__",1)
    dp=f"{G}/designs/{t}_designs_{lab}.json"
    if not os.path.exists(dp): continue
    try:
        d=json.load(open(dp)); g=json.load(open(gf)).get("h3ddg_grade",{})
    except Exception: continue
    nat=d.get("native_cdr")
    if not nat: continue
    ev={v["name"] for v in d.get("antigen_variants",[]) if v.get("split")=="eval"}
    if not ev: continue
    arm=lab.rsplit("_r",1)[0]
    for x in d.get("antibody_designs",[]):
        vals=g.get(x["name"])
        if not vals: continue
        vv=[vals[k] for k in vals if k in ev]
        if len(vv)!=len(ev): continue
        seq=x["cdr_seq"]
        if len(seq)!=len(nat): continue
        ed=sum(1 for p,q in zip(nat,seq) if p!=q)
        rows.append((arm, seq, nat, ed, max(vv)))

print(f"pooled graded designs: {len(rows):,} across {len({r[0] for r in rows})} arms\n")
if len(rows) < 500:
    raise SystemExit("too few")

# feature correlations with the escape margin (negative margin = more robust)
F=[feats(s) for _a,s,_n,_e,_m in rows]
M=[m for *_x,m in rows]
E=[e for *_x,e,_m in rows]
keys=[k for k in F[0] if k!="len"]
print("Spearman-ish (Pearson) of feature vs worst-case margin  (NEGATIVE r = feature associated with MORE robust):")
print(f"  {'feature':10s} {'all designs':>12s} {'edits<=median':>14s} {'edits>median':>13s}")
med=st.median(E)
lo_i=[i for i,e in enumerate(E) if e<=med]; hi_i=[i for i,e in enumerate(E) if e>med]
out=[]
for k in keys:
    r_all=corr([f[k] for f in F], M)
    r_lo=corr([F[i][k] for i in lo_i],[M[i] for i in lo_i])
    r_hi=corr([F[i][k] for i in hi_i],[M[i] for i in hi_i])
    out.append((abs(r_all),k,r_all,r_lo,r_hi))
for _a,k,r,rl,rh in sorted(out, reverse=True):
    print(f"  {k:10s} {r:12.3f} {rl:14.3f} {rh:13.3f}")

# top vs bottom decile by margin, within edit-distance strata to kill the confound
print("\nTop vs bottom decile by robustness, WITHIN matched edit distance:")
byed=defaultdict(list)
for i,(a,s,n,e,m) in enumerate(rows): byed[e].append(i)
agg=defaultdict(lambda:[[],[]])
for e,idx in byed.items():
    if len(idx)<60: continue
    idx=sorted(idx,key=lambda i:M[i])
    k=max(5,len(idx)//10)
    for i in idx[:k]:
        for kk in keys: agg[kk][0].append(F[i][kk])
    for i in idx[-k:]:
        for kk in keys: agg[kk][1].append(F[i][kk])
print(f"  {'feature':10s} {'robust':>9s} {'fragile':>9s} {'delta':>9s}")
deltas=[]
for kk in keys:
    a,b=agg[kk]
    if not a: continue
    deltas.append((abs(st.mean(a)-st.mean(b)), kk, st.mean(a), st.mean(b)))
for _d,kk,a,b in sorted(deltas, reverse=True):
    print(f"  {kk:10s} {a:9.3f} {b:9.3f} {a-b:+9.3f}")
