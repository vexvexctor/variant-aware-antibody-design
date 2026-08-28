# Experiment 5 — Wild-type affinity retention (analysis half)

**Question.** Does the escape-robustness steering of the headline arm come at the
cost of wild-type (WT) binding affinity? I.e. do we buy held-out-variant
robustness by giving up potency on the original antigen?

**Headline arm.** `svdd_af_baddg` (SVDD on the AntiFold prior, BA-DDG reward),
graded across all **51 targets**.

**Grader / data-source decision.** `grades_baddg_long.csv` does **not** contain the
headline SVDD arm — BA-DDG is that arm's *steering objective*, so grading on it
would be in-loop (Goodhart). Per the project grader rule (*steer BA → grade H3*),
the independent per-variant grader for `svdd_af_baddg` is **H3-DDG**, taken from
`grades_h3_long.csv`. Its `h3_ddg` values are already the **design-minus-native
margin** per antigen variant (**< 0 = design tighter than native; lower = better**).
The WT antigen is the `WT` variant; held-out escape variants are the `split == "eval"`
rows (24 held-out variants per design·rep).

**Consistency check.** The per-target *worst-case* (max-variant) margin improves on
58.8% of targets — exactly reproducing the published H3-DDG beat-native headline
(58.8%), confirming the selection logic below matches the pipeline's beat-native
formula.

---

## 1. Per-design metrics (held-out variants)

For each design·rep unit I compute four design-minus-native margins:

| metric | definition | lower = better means |
|---|---|---|
| **WT margin** | `h3_ddg` on the `WT` antigen | design tighter than native on the original antigen |
| **mean-variant** | mean `h3_ddg` over the 24 held-out eval variants | design more robust on average |
| **CVaR20-variant** | mean of the worst 20% (highest ddG, α=0.2 → 5 variants) of eval margins | design more robust in the bad tail |
| **max-variant** | max `h3_ddg` over eval variants (worst case) | design robust even in the single worst escape |

Full per-design·rep table: `results/report_data/wt_retention_per_design.csv`.

**Per-target collapse (shipped design).** For each target I select the design·rep
that would actually be shipped — the **best-of-pool** unit, i.e. the one minimizing
the worst-case (max) held-out-variant margin. This matches the published
beat-native selection (min-over-pool of max-over-eval). The selected unit's four
margins are the per-target values used below.

Per-target medians of the shipped design (N=51):

| WT | mean-variant | CVaR20-variant | max-variant |
|---|---|---|---|
| **-0.852** | -0.760 | -0.520 | -0.419 |

All four medians are negative → the shipped robust design is, on the typical target,
*simultaneously* tighter than native on the WT antigen **and** more robust across
held-out escape variants.

Improvement rates across the 51 targets:
- WT margin improved (< 0): **78.4%**
- CVaR20 robustness improved (< 0): **62.7%**
- worst-case (max) robustness improved (< 0): **58.8%** (= published headline)

---

## 2. Quadrant scatter — the WT ↔ robustness trade-off

x = WT margin, y = held-out CVaR20 margin, one point per target (shipped design),
`svdd_af_baddg`, N=51. A target lands in a quadrant by the signs of the two margins
(negative = improved vs native):

| quadrant | meaning | targets | **%** |
|---|---|---|---|
| **Q1** | WT improved **and** robustness improved | 32 | **62.7%** |
| **Q2** | WT worsened, robustness improved | 0 | **0.0%** |
| **Q3** | WT improved, robustness worsened | 8 | **15.7%** |
| **Q4** | both worsened | 11 | **21.6%** |

**Trade-off correlation.** Spearman between WT margin and CVaR20 margin across the
51 targets: **ρ = +0.928, p = 1.1e-22** (strongly positive).

**Interpretation — no trade-off.** A robustness-vs-WT *cost* would show up as a
*negative* correlation (better robustness ↔ worse WT) and a populated Q2 (robustness
bought by sacrificing WT). Instead the correlation is strongly **positive** and
**Q2 is empty (0%)**: WT binding and escape robustness move together. The headline
arm essentially never gains robustness by giving up WT affinity — when it improves
one it improves the other (Q1, 62.7%), and its failures are joint (Q4, 21.6%, where
the design beats native on neither axis). The 15.7% Q3 targets keep WT gains while
the worst-tail robustness slips slightly; there is no symmetric Q2 population.

**Bottom line:** robustness does **not** come at the cost of WT binding for
`svdd_af_baddg`. The two objectives are aligned (ρ=+0.93), so escape-robustness
steering retains — and on 78% of targets improves — wild-type affinity.

---

## 3. Outputs

- `results/report_data/wt_retention_quadrants.csv` — tidy per-target:
  `target, arm, wt_margin, mean_margin, cvar20_margin, max_margin, quadrant`
- `results/report_data/wt_retention_per_design.csv` — per design·rep unit (supporting)
- `results/figures/exp5_wt_quadrants.png` — labeled 4-quadrant scatter

Analysis script: `scratchpad/exp5.py` (pandas + scipy; matplotlib for the figure).
Grader = H3-DDG (out-of-family vs the BA-DDG steering objective).
