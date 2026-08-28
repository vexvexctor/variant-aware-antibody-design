# Results

Tidy aggregates only. Raw per-design output (~77 GB) is not released; see `REPRODUCIBILITY.md`.

| directory | contents |
|---|---|
| `headline/` | main generator comparison, per-target margins, generator × protein family, WT affinity retention |
| `controls/` | temperature-matched, edit-distance null, objective ablation, panel-size sensitivity, budget-matched search, antibody-cluster split |
| `generalization/` | Pythia-PPI out-of-family, post-cutoff antigens, natural SARS-CoV-2 variants, measured-DMS transfer, AACDB coverage |
| `evaluator_agreement/` | Spearman, Cohen's κ, percent agreement, base rates on the shared 1,785-design pool |
| `trajectories/` | 121 raw guided/unguided decoding traces. **Descriptive only** — partial capture |
| `biology/` | CDR-H3 residue composition, charge complementarity |
| `compute/` | per-arm wall/CUDA seconds, reward-evaluation counters |

Figures in `figures/` are regenerated from these files by `scripts/08_make_figures.sh`. No result
value is typed into a plotting script.

**Before quoting any number**, read `headline/HEADLINE_HONEST.md` (the four controls that
narrowed the claim) and `evaluator_agreement/EVALUATOR_AGREEMENT.md` (which comparisons are
legal given evaluator lineage).
