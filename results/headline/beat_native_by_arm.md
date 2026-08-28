## Beat-native by arm

Rows flagged ⚠️ are de-novo generators: they design a new CDR-H3 backbone (RFantibody also re-docks the Fv), so grading their SEQUENCE on the native complex understates them by construction and their row should not be read as their performance. Their fair number comes from `grade_own_structure.py`.

| arm | beat-native (best-of-pool) | 95% CI | n targets | mean margin | one-shot beat-native | wall s/run | reward evals/run |
|---|---|---|---|---|---|---|---|
| `base_abmpnn_ar` | 77.1% | [70, 83] | 190 | -0.755 | 33.2% | 8 | 0 |
| `base_abmpnn_ar_t1` | 65.3% | [58, 72] | 190 | -0.346 | 5.3% | 8 | 0 |
| `base_abmpnn_npz` | 69.0% | [62, 76] | 190 | -0.547 | 37.9% | 0 | 0 |
| `base_antifold` | 50.1% | [43, 58] | 190 | +0.224 | 17.9% | 0 | 0 |
| `base_antifold_t1` | 46.9% | [40, 55] | 190 | +0.231 | 4.7% | 1 | 0 |
| `base_diffab` ⚠️ de-novo backbone: see own-structure protocol | 15.8% | [8, 25] | 75 | +2.239 | 2.7% | 15 | 0 |
| `base_diffusion` | 1.9% | [0, 4] | 190 | +4.773 | 0.0% | 26 | 0 |
| `base_diffusion_s` | 4.2% | [0, 9] | 66 | +3.135 | 0.0% | 24 | 89 |
| `base_esmif` | 39.0% | [32, 46] | 190 | +0.208 | 17.4% | 0 | 0 |
| `base_esmif_ar` | 41.5% | [34, 49] | 190 | +0.208 | 15.3% | 9 | 0 |
| `base_esmif_ar_t1` | 40.6% | [33, 48] | 189 | +0.177 | 2.6% | 8 | 0 |
| `base_iglm` | 0.0% | [0, 0] | 74 | +6.446 | 0.0% | 2 | 0 |
| `base_proteinmpnn` | 51.6% | [44, 59] | 190 | -0.222 | 28.9% | 0 | 0 |
| `base_proteinmpnn_ar` | 53.5% | [45, 61] | 190 | -0.380 | 23.2% | 7 | 0 |
| `base_proteinmpnn_ar_t1` | 38.4% | [31, 46] | 188 | +0.373 | 2.7% | 8 | 0 |
| `base_rfantibody` ⚠️ de-novo backbone: see own-structure protocol | 0.0% | [0, 0] | 53 | +4.236 | 0.0% | 506 | 0 |
| `grad_diffusion` | 44.3% | [37, 52] | 190 | +0.128 | 15.3% | 26 | 0 |
| `obj_mint` | 49.3% | [37, 62] | 67 | +0.224 | 4.5% | 129 | 411 |
| `obj_wt` | 56.2% | [44, 69] | 66 | -0.147 | 7.6% | 19 | 410 |
| `ours_abmpnn_live` | 83.4% | [78, 89] | 190 | -0.907 | 22.1% | 89 | 463 |
| `ours_abmpnn_npz` | 79.9% | [74, 86] | 190 | -0.835 | 27.9% | 89 | 463 |
| `ours_antifold` | 64.3% | [57, 72] | 190 | -0.466 | 14.2% | 88 | 462 |
| `ours_diffusion` | 11.9% | [7, 17] | 190 | +2.106 | 0.5% | 95 | 505 |
| `ours_esmif` | 61.6% | [54, 69] | 190 | -0.411 | 9.5% | 90 | 478 |
| `ours_proteinmpnn` | 63.5% | [56, 71] | 190 | -0.536 | 13.7% | 94 | 502 |
| `ours_proteinmpnn_live` | 62.6% | [55, 70] | 190 | -0.434 | 11.6% | 96 | 501 |
| `rand_base_abmpnn_ar` | 34.0% | [27, 41] | 190 | +0.527 | 0.0% | nan | nan |
| `rand_base_antifold` | 29.8% | [23, 37] | 190 | +0.719 | 0.0% | nan | nan |
| `rand_base_proteinmpnn_ar` | 7.9% | [4, 12] | 190 | +1.571 | 0.0% | nan | nan |
| `rand_ours_abmpnn_live` | 20.4% | [15, 27] | 190 | +0.883 | 0.0% | nan | nan |
| `rand_ours_antifold` | 21.9% | [16, 28] | 190 | +0.932 | 0.0% | nan | nan |
| `rand_ours_proteinmpnn_live` | 7.2% | [3, 11] | 190 | +1.921 | 0.0% | nan | nan |

## Paired: generator + our objective vs generator alone
