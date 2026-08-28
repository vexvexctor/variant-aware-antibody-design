# Paper

**The manuscript source is not in this archive.**

No `.tex`, `.bib` or `main.pdf` for this paper exists on the machine this repository was
assembled from. Rather than reconstruct it, this directory is left empty by design.

To complete the repository for release, copy in:

```
paper/main.tex  main.pdf  appendix.tex  refs.bib  figs/  tables/
```

Consequences while it is missing:

* The manuscript's **Figure 1** (method schematic) is unavailable. `README.md` documents the
  method with pseudocode instead.
* Appendix figures that were produced outside this archive — trajectory summary, biochemical/
  residue figure, target-level native figure — are not present. The **result files** behind them
  are, under `results/`, so they can be redrawn.
* Figures that *are* regenerated from released CSVs live in `figures/`.

`REPRODUCIBILITY.md` maps every manuscript claim to the code and data that produce it, including
the ones whose figures are missing.
