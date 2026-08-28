# EXP4 — CVaR cross-objective matrix (antibody escape-robustness)

Design-time objective = `--steer-agg {mean,worst,cvar} --cvar-alpha X` on the SAME generator (SVDD, base=antifold, oracle=baddg_ddg, --svdd-k 6 --svdd-temp 1.0) so the matrix isolates the *objective*, not the vehicle. 51-target panel, 3 unseeded reps.

**FROZEN protocol (per Exp2):** one deploy design per (target,rep) = the pool member minimising BA-DDG CVaR20 over the OPT/steer variants; then grade THAT design's HELD-OUT (eval-split) variants with H3-DDG (out-of-family grader — correct, since steering is BA-DDG). Margin<0 = design binds the held-out escape variant tighter than native.

**Flag mapping (objective -> arm):**

| objective | flags | arm | new? |
|---|---|---|---|
| WT-only | steer set = WT antigen alone (the `single` design; steered on `[(base,1.0)]`, objective-independent) | scr_svdd_af_baddg `single` | reuse |
| mean(tau=1.0) | `--steer-agg mean` (uniform mean over steer panel; SVDD has no separate tau knob — the mean is the full-panel mean) | exp4_mean | NEW |
| CVaR50 | `--steer-agg cvar --cvar-alpha 0.5` | exp4_cvar50 | NEW |
| CVaR20 | `--steer-agg cvar --cvar-alpha 0.2` (= headline) | scr_svdd_af_baddg | reuse |
| CVaR10 | `--steer-agg cvar --cvar-alpha 0.1` | exp4_cvar10 | NEW |
| strict-max(k=1) | `--steer-agg worst` | exp4_worst | NEW |

## Held-out H3 margin matrix
row = design-time objective; cols = eval metric on held-out H3 margins (median over targets of rep-mean). Lower = tighter/more-robust; <0 = beats native. `beatN` = fraction beating native worst-case (targ = rep-averaged; run = per (target,rep)). `WTmarg` = median held-out margin on the WT antigen. `editN` = mean CDR-H3 edits from native.

| objective | worst | cvar10 | cvar20 | cvar50 | mean | beatN(targ) | beatN(run) | WTmarg | editN | nT | nRuns |
|---|---|---|---|---|---|---|---|---|---|---|---|
| WT-only | 0.863 | 0.703 | 0.650 | 0.493 | 0.318 | 0.098 | 0.216 | 0.223 | 6.58 | 51 | 153 |
| mean(tau=1.0) | 0.146 | 0.046 | 0.006 | -0.120 | -0.370 | 0.431 | 0.399 | -0.324 | 6.46 | 51 | 153 |
| CVaR50(a=0.5) | 0.274 | 0.191 | 0.113 | -0.046 | -0.233 | 0.373 | 0.327 | -0.184 | 6.38 | 51 | 153 |
| CVaR20(a=0.2) | 0.258 | 0.085 | 0.011 | -0.143 | -0.322 | 0.451 | 0.386 | -0.304 | 6.45 | 51 | 153 |
| CVaR10(a=0.1) | 0.132 | 0.069 | -0.012 | -0.146 | -0.355 | 0.431 | 0.444 | -0.262 | 6.35 | 51 | 153 |
| strict-max(k1) | 0.213 | 0.082 | 0.041 | -0.036 | -0.204 | 0.314 | 0.366 | -0.236 | 6.25 | 51 | 153 |

## Sequence stability across objectives
Pairwise at matched (target, first shared rep): mean CDR-H3 edit distance, fraction of targets choosing an IDENTICAL frozen sequence, mean Jaccard of mutation sets (vs native).

| obj A | obj B | n | mean CDR-H3 edit | frac identical | mut-set Jaccard |
|---|---|---|---|---|---|
| WT-only | mean(tau=1.0) | 51 | 6.55 | 0.00 | 0.62 |
| WT-only | CVaR50(a=0.5) | 51 | 6.33 | 0.00 | 0.64 |
| WT-only | CVaR20(a=0.2) | 51 | 5.24 | 0.16 | 0.69 |
| WT-only | CVaR10(a=0.1) | 51 | 6.27 | 0.00 | 0.65 |
| WT-only | strict-max(k1) | 51 | 6.12 | 0.00 | 0.65 |
| mean(tau=1.0) | CVaR50(a=0.5) | 51 | 5.78 | 0.00 | 0.65 |
| mean(tau=1.0) | CVaR20(a=0.2) | 51 | 5.88 | 0.00 | 0.65 |
| mean(tau=1.0) | CVaR10(a=0.1) | 51 | 5.75 | 0.02 | 0.65 |
| mean(tau=1.0) | strict-max(k1) | 51 | 5.82 | 0.00 | 0.64 |
| CVaR50(a=0.5) | CVaR20(a=0.2) | 51 | 5.63 | 0.00 | 0.69 |
| CVaR50(a=0.5) | CVaR10(a=0.1) | 51 | 6.14 | 0.00 | 0.63 |
| CVaR50(a=0.5) | strict-max(k1) | 51 | 6.14 | 0.00 | 0.65 |
| CVaR20(a=0.2) | CVaR10(a=0.1) | 51 | 5.73 | 0.00 | 0.70 |
| CVaR20(a=0.2) | strict-max(k1) | 51 | 5.98 | 0.04 | 0.68 |
| CVaR10(a=0.1) | strict-max(k1) | 51 | 6.02 | 0.00 | 0.65 |

CSVs: `results/report_data/exp4_cross_objective_matrix.csv`, `..._seq_stability.csv`.
