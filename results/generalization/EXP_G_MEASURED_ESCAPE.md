# EXP_G — Adversarial selection vs REAL MEASURED escape (AbAgym SARS-CoV-2 RBD)

**Question.** Do designs picked by a per-candidate **adversarial** objective retain binding against the **top DMS-measured escape** mutations (ground truth) better than designs picked by the **fixed-panel** (in-house omega escape panel) or **in-family** (potency-only, `wt_only`) baselines?

**Ground truth.** AbAgym anti-RBD deep-mutational-scanning; per (site,substitution) the MAX MinMax-normalized DMS escape across the anti-RBD DMS panel; top-18 measured escapes restricted to each antibody's structural epitope (<=5A heavy-atom contacts). WT residue identity read from each target PDB (Wuhan RBD numbering).

**Designs.** campaign2/B1 objective-selected CDR-H3 picks (adv_worst/adv_cvar vs omega_worst/omega_cvar vs wt_only) x 3 seeds, graded on the AACDB antibody-RBD complex.

**Metric.** retention_loss(d,e)=E[d][e]-E[d][WT_RBD] (FoldX interaction ddG; **lower = binding retained = more robust**). Contrast Delta = adversarial - baseline; **Delta<0 => adversarial more robust**. 95% CI = cluster bootstrap over targets.

Systems graded: 7/7 — 7MZK_HLA, 7OR9_HLE, 7TBF_HLA, 7VYR_HLR, 7WPH_HLB, 7YCL_HLD, 8IDN_HLA

## Pooled result

| contrast | metric | pooled Delta (kcal/mol) | 95% CI | adv win-rate | n |
|---|---|---|---|---|---|
| adv_worst_vs_fixed(omega_worst) | worst | +1.261 | [+0.259, +2.722] | 14% | 21 |
| adv_worst_vs_fixed(omega_worst) | mean | +0.138 | [-0.118, +0.354] | 19% | 21 |
| adv_worst_vs_infamily(wt_only) | worst | +1.065 | [+0.143, +2.104] | 33% | 21 |
| adv_worst_vs_infamily(wt_only) | mean | +0.114 | [+0.009, +0.196] | 33% | 21 |
| adv_cvar_vs_fixed(omega_cvar) | worst | +0.983 | [+0.146, +2.008] | 14% | 21 |
| adv_cvar_vs_fixed(omega_cvar) | mean | +0.093 | [-0.194, +0.369] | 19% | 21 |
| adv_cvar_vs_infamily(wt_only) | worst | +1.046 | [+0.296, +1.841] | 24% | 21 |
| adv_cvar_vs_infamily(wt_only) | mean | +0.122 | [-0.061, +0.370] | 29% | 21 |

## Per-system worst-case retention_loss (mean over 3 seeds, lower=robust)

| target | n_esc | native | adv_worst | omega_worst | wt_only | adv-omega | adv-wt |
|---|---|---|---|---|---|---|---|
| 7MZK_HLA | 18 | +9.05 | +5.67 | +5.29 | +5.39 | +0.38 | +0.28 |
| 7OR9_HLE | 18 | +16.05 | +14.39 | +12.44 | +13.09 | +1.94 | +1.29 |
| 7TBF_HLA | 18 | +6.87 | +6.60 | +1.38 | +3.06 | +5.22 | +3.53 |
| 7VYR_HLR | 18 | +8.76 | +4.21 | +3.64 | +2.48 | +0.57 | +1.73 |
| 7WPH_HLB | 18 | +7.29 | +7.60 | +7.34 | +6.13 | +0.26 | +1.47 |
| 7YCL_HLD | 18 | +14.65 | +4.37 | +4.60 | +4.73 | -0.23 | -0.37 |
| 8IDN_HLA | 18 | +2.54 | +2.16 | +1.47 | +2.64 | +0.69 | -0.48 |

_Delta<0 (adv columns lower than baseline) = adversarial-selected design loses less binding against the real measured escapes._
