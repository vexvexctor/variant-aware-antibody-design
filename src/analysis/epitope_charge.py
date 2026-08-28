#!/usr/bin/env python3
"""Is the charge effect absolute, or COMPLEMENTARITY to the epitope?

Robust designs carry more negative net charge at matched edit distance. That could be a generic
bias of the grader, or real electrostatic complementarity. Complementarity predicts the SIGN of
the effect should FLIP with the epitope's own charge: against an acidic epitope, positively
charged CDRs should win. Absolute bias predicts no flip.
"""
from __future__ import annotations
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)
import json, glob, os, statistics as st
from collections import defaultdict

G=f"{VAAD_ROOT}/scratch/genzoo"
PDB=f"{VAAD_ROOT}/datasets/AACDB/data_zip/complex_structure"
CHG={"D":-1,"E":-1,"K":1,"R":1}
AA3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H",
     "ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
     "TYR":"Y","VAL":"V"}

def epitope_charge(target, tmpl):
    """Net charge of antigen residues whose CB (or CA) lies within 10A of any antibody CB/CA.

    CB-level rather than all-atom: the all-atom version is O(1e7) pairs per target and does not
    finish in reasonable time across 200 targets, and the residue-level contact set is
    essentially unchanged at a slightly larger cutoff.
    """
    import numpy as np
    p=os.path.join(PDB, f"{target}.pdb")
    if not os.path.exists(p): return None
    ab=set(tmpl["antibody_chains"]); ag=set(tmpl["antigen_chains"])
    A=[]; Bc=[]; Ba=[]
    for line in open(p):
        if not line.startswith("ATOM"): continue
        an=line[12:16].strip()
        if an not in ("CB","CA"): continue
        ch=line[21]
        xyz=(float(line[30:38]),float(line[38:46]),float(line[46:54]))
        if ch in ab: A.append(xyz)
        elif ch in ag:
            Bc.append(xyz); Ba.append((int(line[22:26]), AA3.get(line[17:20].strip(),"X")))
    if not A or not Bc: return None
    A=np.asarray(A,dtype=np.float32); B=np.asarray(Bc,dtype=np.float32)
    d2=((B[:,None,:]-A[None,:,:])**2).sum(-1).min(1)
    seen={}
    for i,ok in enumerate(d2<=100.0):
        if ok and Ba[i][0] not in seen: seen[Ba[i][0]]=Ba[i][1]
    return sum(CHG.get(a,0) for a in seen.values()), len(seen)

cache={}
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
    nat=d.get("native_cdr"); 
    if not nat: continue
    if t not in cache:
        cache[t]=epitope_charge(t,d)
    ec=cache[t]
    if ec is None: continue
    ev={v["name"] for v in d.get("antigen_variants",[]) if v.get("split")=="eval"}
    if not ev: continue
    for x in d.get("antibody_designs",[]):
        vals=g.get(x["name"])
        if not vals: continue
        vv=[vals[k] for k in vals if k in ev]
        if len(vv)!=len(ev): continue
        s=x["cdr_seq"]
        if len(s)!=len(nat): continue
        rows.append((ec[0], sum(CHG.get(c,0) for c in s),
                     sum(1 for p,q in zip(nat,s) if p!=q), max(vv)))

def corr(xs,ys):
    n=len(xs); mx,my=st.mean(xs),st.mean(ys)
    num=sum((a-mx)*(b-my) for a,b in zip(xs,ys))
    den=(sum((a-mx)**2 for a in xs)*sum((b-my)**2 for b in ys))**0.5
    return num/den if den else float("nan")

print(f"designs with an epitope-charge estimate: {len(rows):,}")
eps=sorted({r[0] for r in rows})
print(f"epitope net charge ranges {min(eps)} .. {max(eps)}\n")
print("Correlation of DESIGN charge vs escape margin, split by EPITOPE charge")
print("(negative r = more positive CDR charge is MORE robust)\n")
print(f"  {'epitope charge':>16s} {'n designs':>10s} {'r(design charge, margin)':>26s}")
for lo,hi,lab in [(-99,-3,"acidic (<= -3)"),(-2,2,"neutral (-2..+2)"),(3,99,"basic (>= +3)")]:
    sub=[r for r in rows if lo<=r[0]<=hi]
    if len(sub)<200: 
        print(f"  {lab:>16s} {len(sub):10d}  (too few)"); continue
    print(f"  {lab:>16s} {len(sub):10d} {corr([s[1] for s in sub],[s[3] for s in sub]):26.3f}")
