# Computational overhead by arm

All runs on one H100 per task, identical accounting via `profile_util.Profiler`. TRAINING COST IS ZERO FOR EVERY ARM HERE: no arm fine-tunes, trains a reward model, or distills anything -- the priors and the BA-DDG reward are frozen external checkpoints. So the whole comparison is inference-time cost.


## Marginal cost per run (one target, one rep, pool of 6 designs)

| arm | n runs | wall s | wall s / design | CUDA s | reward evals | reward evals / design | base forwards | peak VRAM GB |
|---|---|---|---|---|---|---|---|---|
| `base_abmpnn_npz` | 600 | 0.03 | 0.005 | 0.00 | 0 | 0 | 0 | 0.00 |
| `base_proteinmpnn` | 600 | 0.03 | 0.006 | 0.00 | 0 | 0 | 0 | 0.00 |
| `base_antifold` | 600 | 0.04 | 0.006 | 0.00 | 0 | 0 | 0 | 0.00 |
| `base_esmif` | 600 | 0.05 | 0.008 | 0.00 | 0 | 0 | 0 | 0.00 |
| `base_antifold_t1` | 600 | 1.24 | 0.207 | 0.02 | 0 | 0 | 0 | 0.00 |
| `base_iglm` | 475 | 2.27 | 0.379 | 0.00 | 0 | 0 | 0 | nan |
| `base_proteinmpnn_ar` | 600 | 7.26 | 1.209 | 6.46 | 0 | 0 | 6 | 0.48 |
| `base_esmif_ar_t1` | 600 | 7.97 | 1.328 | 0.00 | 0 | 0 | 82 | nan |
| `base_proteinmpnn_ar_t1` | 600 | 8.30 | 1.383 | 7.74 | 0 | 0 | 6 | 0.48 |
| `base_abmpnn_ar_t1` | 600 | 8.35 | 1.392 | 7.64 | 0 | 0 | 6 | 0.48 |
| `base_esmif_ar` | 600 | 8.49 | 1.415 | 0.00 | 0 | 0 | 82 | nan |
| `base_abmpnn_ar` | 600 | 8.61 | 1.434 | 7.83 | 0 | 0 | 8 | 0.48 |
| `base_diffab` | 552 | 14.74 | 2.457 | 0.00 | 0 | 0 | 0 | nan |
| `obj_wt` | 567 | 21.30 | 3.549 | 12.77 | 456 | 76 | 95 | 0.93 |
| `base_diffusion_s` | 325 | 24.96 | 4.160 | 16.03 | 95 | 16 | 96 | 0.94 |
| `grad_diffusion` | 600 | 26.15 | 4.359 | 15.90 | 0 | 0 | 0 | 2.79 |
| `base_diffusion` | 600 | 26.32 | 4.387 | 15.98 | 0 | 0 | 0 | 2.79 |
| `ours_antifold` | 600 | 87.65 | 14.608 | 76.57 | 457 | 76 | 95 | 0.94 |
| `ours_abmpnn_npz` | 600 | 88.12 | 14.687 | 77.16 | 459 | 76 | 95 | 0.94 |
| `ours_abmpnn_live` | 600 | 88.99 | 14.832 | 77.87 | 460 | 77 | 95 | 0.95 |
| `ours_esmif` | 600 | 89.51 | 14.918 | 78.56 | 473 | 79 | 95 | 0.94 |
| `ours_proteinmpnn` | 600 | 93.79 | 15.632 | 82.74 | 498 | 83 | 95 | 0.94 |
| `ours_diffusion` | 600 | 93.98 | 15.664 | 82.91 | 501 | 84 | 95 | 0.94 |
| `ours_proteinmpnn_live` | 600 | 95.50 | 15.917 | 84.11 | 498 | 83 | 95 | 0.95 |
| `ours_esmif_live` | 570 | 106.79 | 17.798 | 87.95 | 477 | 80 | 95 | 0.94 |
| `obj_mint` | 504 | 155.21 | 25.868 | 138.55 | 465 | 78 | 96 | 6.70 |
| `base_rfantibody` | 53 | 506.09 | 84.348 | 0.00 | 0 | 0 | 0 | nan |

## One-off amortised cost: per-target logit-table export

Paid once per target, then reused by every arm and every rep that uses that base. For the table-based bases (AntiFold, and the frozen-npz ablation arms) this export IS the generator's real inference cost -- the 0.0x s 'baseline' below is only the cost of drawing samples from an already-computed table, so read the two together.

| export | n targets | wall s | CUDA s |
|---|---|---|---|
| `export_abmpnn_npz` | 200 | 2.28 | 0.44 |
| `export_proteinmpnn` | 200 | 2.24 | 0.42 |

## Overhead of adding our objective to a generator

| generator | baseline wall s | + our objective wall s | multiplier | added reward evals |
|---|---|---|---|---|
| antifold | 0.04 | 87.65 | 2391x | 457 |
| ProteinMPNN | 7.26 | 95.50 | 13x | 498 |
| AbMPNN | 8.61 | 88.99 | 10x | 460 |
| ESM-IF1 | 8.49 | 89.51 | 11x | 473 |
| our diffusion | 26.32 | 93.98 | 4x | 501 |
