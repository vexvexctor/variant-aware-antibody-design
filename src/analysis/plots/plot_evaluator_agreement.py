#!/usr/bin/env python3
"""Grader agreement: score-level (Spearman) and decision-level (Cohen's kappa).
Both encode magnitude of agreement -> sequential single-hue ramp, light->dark."""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

INK, INK2, SURF, GRIDC = "#0b0b0b", "#52514e", "#fcfcfb", "#e8e7e3"
seq = LinearSegmentedColormap.from_list("seq", ["#f4f8fd", "#2a78d6", "#12365f"])
ORDER = ["BA-DDG", "H3", "Pythia", "FoldX", "MM-GBSA"]
rows = list(csv.DictReader(open("results/evaluator_agreement/pairwise_long.csv")))
M = {(r["evaluator_a"], r["evaluator_b"]): r for r in rows}
def g(a, b, k):
    r = M.get((a, b)) or M.get((b, a))
    return float(r[k]) if r else np.nan

fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4))
for ax, key, title, vmax in [
    (axes[0], "spearman", "Score agreement  (Spearman of held-out margins)", 0.70),
    (axes[1], "cohen_kappa",    "Decision agreement  (Cohen's κ on the beat-native call)", 0.70)]:
    Z = np.full((5, 5), np.nan)
    for i, a in enumerate(ORDER):
        for j, b in enumerate(ORDER):
            if i != j: Z[i, j] = g(a, b, key)
    im = ax.imshow(Z, cmap=seq, vmin=0, vmax=vmax, aspect="equal")
    for i in range(5):
        for j in range(5):
            if i == j:
                ax.add_patch(plt.Rectangle((j-.5, i-.5), 1, 1, facecolor=GRIDC, edgecolor=SURF, lw=2))
                ax.text(j, i, "—", ha="center", va="center", fontsize=11, color=INK2); continue
            v = Z[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=11.5,
                    color="#ffffff" if v > .42 else INK,
                    fontweight="bold" if v >= .30 else "normal")
    ax.set_xticks(range(5)); ax.set_yticks(range(5))
    ax.set_xticklabels(ORDER, fontsize=9.5, color=INK, rotation=20, ha="right")
    ax.set_yticklabels(ORDER, fontsize=9.5, color=INK)
    ax.set_xticks(np.arange(-.5, 5, 1), minor=True); ax.set_yticks(np.arange(-.5, 5, 1), minor=True)
    ax.grid(which="minor", color=SURF, lw=2.5); ax.tick_params(which="minor", length=0)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title(title, loc="left", fontsize=11, color=INK, pad=10)
    cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.03)
    cb.ax.tick_params(labelsize=8, colors=INK2); cb.outline.set_visible(False)

fig.suptitle("Do our evaluators agree?  1,785 identical designs, every evaluator scoring the same units",
             fontsize=12.5, color=INK, x=.008, ha="left", y=.985)
fig.text(.008, .915, "Bold ≥ 0.30. κ is chance-corrected: FoldX↔MM-GBSA agree on 82% of calls but "
         "κ = −0.01 — they agree only by both saying “no”.",
         fontsize=8.8, color=INK2)
fig.patch.set_facecolor(SURF)
for ax in axes: ax.set_facecolor(SURF)
fig.tight_layout(rect=[0, 0, 1, .90])
fig.savefig("figures/fig_evaluator_agreement.png", dpi=200, facecolor=SURF)
print("wrote grader_agreement.png")
