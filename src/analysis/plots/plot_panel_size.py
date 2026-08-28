#!/usr/bin/env python3
"""Beat-native and held-out margin vs design-panel size M, and vs reward-evaluation cost."""
import csv, json, os, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCR="results/controls"
rows=list(csv.DictReader(open(f"{SCR}/panel_size_summary.csv")))
M=[int(r["M_design"]) for r in rows]
bn=[float(r["beat_native_pooled"]) for r in rows]
lo=[float(r["ci_lo"]) for r in rows]; hi=[float(r["ci_hi"]) for r in rows]
oc=[float(r["oracle_calls"]) for r in rows]
err=[[b-l for b,l in zip(bn,lo)],[h-b for b,h in zip(bn,hi)]]

# mean held-out margin per M (the continuous metric, where the signal is)
per=list(csv.DictReader(open(f"{SCR}/panel_size_per_unit.csv")))
marg={m:[] for m in M}
for r in per: marg[int(r["M_design"])].append(float(r["best_of_pool_worst_margin"]))
mm=[st.mean(marg[m]) for m in M]

INK,ACC,ACC2,GRID="#1a1a1a","#2f6f9f","#b1543a","#d8d8d8"
fig,ax=plt.subplots(1,3,figsize=(15,4.6))

ax[0].errorbar(M,bn,yerr=err,marker="o",ms=7,lw=2,capsize=4,color=ACC,ecolor=ACC,mfc="white",mew=2)
ax[0].set_ylabel("beat-native (%)"); ax[0].set_title("Beat-native vs panel size — flat",loc="left",fontsize=11)
for m,b in zip(M,bn): ax[0].annotate(f"{b:.0f}%",(m,b),textcoords="offset points",xytext=(0,11),ha="center",fontsize=9)

ax[1].plot(M,mm,marker="o",ms=7,lw=2,color=ACC2,mfc="white",mew=2)
ax[1].set_ylabel("mean held-out worst-case margin\n(lower = tighter than native)")
ax[1].set_title("Margin vs panel size — M=12 significantly better",loc="left",fontsize=11)
ax[1].set_ylim(min(mm)-0.012, max(mm)+0.016)
for m,v in zip(M,mm): ax[1].annotate(f"{v:.3f}",(m,v),textcoords="offset points",xytext=(0,11),ha="center",fontsize=9)

for a in ax[:2]:
    a.set_xscale("log"); a.set_xticks(M); a.set_xticklabels([str(m) for m in M])
    a.set_xlabel("design panel size  $M_{design}$")
    a.axvline(6,color="#b0b0b0",ls="--",lw=1,zorder=0)

ax[2].errorbar(oc,bn,yerr=err,marker="o",ms=7,lw=2,capsize=4,color=ACC,ecolor=ACC,mfc="white",mew=2)
for m,x,b in zip(M,oc,bn):
    ax[2].annotate(f"M={m}",(x,b),textcoords="offset points",xytext=(0,11),ha="center",fontsize=9)
ax[2].set_xlabel("reward-evaluation cost (oracle calls / design run)")
ax[2].set_ylabel("beat-native (%)")
ax[2].set_title("Performance vs cost — no return on spend",loc="left",fontsize=11)

for a in ax:
    a.grid(True,color=GRID,lw=.6,alpha=.7); a.set_axisbelow(True)
    for s in ("top","right"): a.spines[s].set_visible(False)
    for s in ("left","bottom"): a.spines[s].set_color("#999")
    a.tick_params(colors="#555",labelsize=9)

fig.suptitle("Design-panel-size sensitivity — 51 targets x 3 reps, H3-DDG (out-of-family), "
             "identical 24 held-out variants at every M; CI = cluster bootstrap (30 antigen clusters)",
             fontsize=10,color=INK,x=.01,ha="left",y=.99)
fig.tight_layout(rect=[0,0,1,.93])
fig.savefig("figures/fig_panel_size.png",dpi=200,facecolor="white")
print("wrote panel_size_sensitivity.png")
for m,b,v,x in zip(M,bn,mm,oc): print(f"M={m:>2}  beat-native={b:>5.1f}%  mean margin={v:>+7.3f}  oracle_calls={x:>6.0f}")
