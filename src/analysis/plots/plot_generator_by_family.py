#!/usr/bin/env python3
"""Diverging heatmap: paired improvement (ours - its own base) in percentage points,
generator x protein family. Polarity is the data's job -> diverging scale,
blue (improvement) / gray 0 / red (decline), per the reference palette."""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np

BLUE, GRAY, RED = "#2a78d6", "#f0efec", "#e34948"
INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
cmap = LinearSegmentedColormap.from_list("div", [RED, GRAY, BLUE])

rows = list(csv.DictReader(open("results/headline/generator_by_protein_family.csv")))
GENS = ["AntiFold", "AbMPNN", "ProteinMPNN", "ESM-IF1", "Masked diffusion", "RFantibody (de-novo)"]
FAMS = ["Coronavirus spike", "Flavivirus envelope (E)", "Influenza hemagglutinin",
        "— all viral glycoprotein", "— all 190 targets"]
idx = {(r["generator"], r["family"]): r for r in rows}

Z = np.full((len(GENS), len(FAMS)), np.nan)
for i, g in enumerate(GENS):
    for j, f in enumerate(FAMS):
        if (g, f) in idx: Z[i, j] = float(idx[(g, f)]["dpp"])

fig, ax = plt.subplots(figsize=(11.4, 5.0))
lim = np.nanmax(np.abs(Z))
im = ax.imshow(Z, cmap=cmap, norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim), aspect="auto")

for i, g in enumerate(GENS):
    for j, f in enumerate(FAMS):
        r = idx.get((g, f))
        if not r: continue
        d, n = float(r["dpp"]), int(r["n"])
        sig = float(r["lo"]) > 0 or float(r["hi"]) < 0
        # ink tokens for text, never the series/cell color
        ax.text(j, i - .13, f"{d:+.1f}", ha="center", va="center", fontsize=11.5,
                color=INK, fontweight="bold" if sig else "normal")
        ax.text(j, i + .22, f"n={n}" + ("  ✓" if sig else ""), ha="center", va="center",
                fontsize=8, color=INK2)
        if sig:  # 2px surface ring marks the resolved cells
            ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                       edgecolor=INK, lw=1.6, zorder=3))

ax.set_xticks(range(len(FAMS)))
ax.set_xticklabels(["Coronavirus\nspike", "Flavivirus\nenvelope (E)", "Influenza\nhemagglutinin",
                    "all viral\nglycoprotein", "all 190\ntargets"], fontsize=9.5, color=INK)
ax.set_yticks(range(len(GENS))); ax.set_yticklabels(GENS, fontsize=10, color=INK)
ax.axvline(2.5, color=SURF, lw=4)          # 2px+ surface gap: families | aggregates
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
ax.set_xticks(np.arange(-.5, len(FAMS), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(GENS), 1), minor=True)
ax.grid(which="minor", color=SURF, lw=2.5); ax.tick_params(which="minor", length=0)

cb = fig.colorbar(im, ax=ax, pad=.02, fraction=.033)
cb.set_label("paired Δ beat-native (percentage points), ours − its own base", fontsize=9, color=INK2)
cb.ax.tick_params(labelsize=8, colors=INK2); cb.outline.set_visible(False)

ax.set_title("Improvement from our objective, by generator and antigen protein family",
             loc="left", fontsize=12.5, color=INK, pad=34)
fig.text(.005, .875, "H3-DDG grader (out-of-family), held-out escape panel. Bold + ring = 95% "
         "bootstrap CI excludes zero. Family cells are small-n; the two right-hand columns are "
         "aggregates, not families.", fontsize=8.6, color=INK2)
fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
fig.tight_layout(rect=[0, 0, 1, .97])
fig.savefig("figures/fig_generator_by_family.png", dpi=200, facecolor=SURF)
print("wrote family_breakdown.png")
