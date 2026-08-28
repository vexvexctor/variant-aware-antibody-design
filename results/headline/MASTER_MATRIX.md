# Generator zoo — master matrix (H3-DDG grader, held-out escape panel)

Reference = the NATIVE affinity-matured antibody (H3-DDG margin 0 by construction). `beat-native` = best-of-pool worst-case margin < 0 on variants never steered against. CI = 95% bootstrap over ANTIGEN CLUSTERS (mmseqs-70), so N is independent antigens, not PDB entries.

Panel: 200 targets, 172 antigen clusters. Graded so far: 190 targets.


## Beat-native by arm

Rows flagged ⚠️ are de-novo generators: they design a new CDR-H3 backbone (RFantibody also re-docks the Fv), so grading their SEQUENCE on the native complex understates them by construction and their row should not be read as their performance. Their fair number comes from `grade_own_structure.py`.

| arm | beat-native (best-of-pool) | 95% CI | n targets | mean margin | one-shot beat-native | wall s/run | reward evals/run |
|---|---|---|---|---|---|---|---|
| `base_abmpnn_ar` | 77.1% | [70, 83] | 190 | -0.755 | 33.2% | 8 | 0 |
| `base_abmpnn_npz` | 69.0% | [62, 76] | 190 | -0.547 | 37.9% | 0 | 0 |
| `base_antifold` | 50.1% | [43, 58] | 190 | +0.224 | 17.9% | 0 | 0 |
| `base_diffab` ⚠️ de-novo backbone: see own-structure protocol | 0.0% | [0, 0] | 1 | +1.690 | 0.0% | 33 | 0 |
| `base_diffusion` | 1.9% | [0, 4] | 190 | +4.867 | 0.0% | 26 | 0 |
| `base_esmif` | 39.2% | [32, 47] | 189 | +0.212 | 17.5% | 0 | 0 |
| `base_esmif_ar` | 41.5% | [34, 49] | 190 | +0.208 | 15.3% | 9 | 0 |
| `base_proteinmpnn` | 51.3% | [44, 59] | 187 | -0.208 | 28.3% | 0 | 0 |
| `base_proteinmpnn_ar` | 53.5% | [45, 61] | 190 | -0.380 | 23.2% | 7 | 0 |
| `base_rfantibody` ⚠️ de-novo backbone: see own-structure protocol | 0.0% | [0, 0] | 53 | +4.236 | 0.0% | 506 | 0 |
| `grad_diffusion` | 43.6% | [36, 51] | 188 | +0.168 | 14.4% | 26 | 0 |
| `obj_mint` | 100.0% | [100, 100] | 1 | -1.477 | 100.0% | 48 | 179 |
| `obj_wt` | 100.0% | [100, 100] | 1 | -1.500 | 100.0% | 12 | 173 |
| `ours_abmpnn_live` | 83.4% | [78, 89] | 190 | -0.907 | 22.1% | 89 | 463 |
| `ours_abmpnn_npz` | 79.8% | [73, 86] | 189 | -0.812 | 27.5% | 89 | 463 |
| `ours_antifold` | 64.3% | [57, 72] | 190 | -0.466 | 14.2% | 88 | 462 |
| `ours_diffusion` | 11.4% | [7, 16] | 189 | +2.123 | 0.5% | 95 | 504 |
| `ours_esmif` | 61.6% | [54, 69] | 190 | -0.409 | 9.5% | 90 | 478 |
| `ours_proteinmpnn` | 63.3% | [56, 71] | 189 | -0.538 | 13.8% | 95 | 503 |
| `ours_proteinmpnn_live` | 62.6% | [55, 70] | 190 | -0.434 | 11.6% | 96 | 501 |
| `rand_base_abmpnn_ar` | 34.0% | [27, 41] | 190 | +0.527 | 0.0% | nan | nan |
| `rand_base_antifold` | 29.8% | [23, 37] | 190 | +0.719 | 0.0% | nan | nan |
| `rand_base_proteinmpnn_ar` | 7.9% | [4, 12] | 190 | +1.571 | 0.0% | nan | nan |
| `rand_ours_abmpnn_live` | 20.4% | [15, 27] | 190 | +0.883 | 0.0% | nan | nan |
| `rand_ours_antifold` | 21.9% | [16, 28] | 190 | +0.932 | 0.0% | nan | nan |
| `rand_ours_proteinmpnn_live` | 7.2% | [3, 11] | 190 | +1.921 | 0.0% | nan | nan |

## Paired: generator + our objective vs generator alone

Same base distribution, same pool size, same panel — the only difference is the objective. Positive pp = our objective helps that generator.

| generator | base alone | + our objective | Δ pp | Δ CI | n paired | Δ mean margin | extra wall s | extra reward evals |
|---|---|---|---|---|---|---|---|---|
| AntiFold | 50.1% | 64.3% | +14.2 ** | [+8, +21] | 190 | -0.690 | +88 | +462 |
| AbMPNN | 77.1% | 83.4% | +6.3 ** | [+0, +12] | 190 | -0.152 | +81 | +463 |
| ProteinMPNN | 53.5% | 62.0% | +8.5 ** | [+2, +15] | 190 | -0.054 | +89 | +501 |
| our masked diffusion | 1.9% | 10.8% | +8.9 ** | [+4, +14] | 189 | -2.867 | +68 | +504 |

### Ablation: same models with a FROZEN prior on both sides

Isolates what re-conditioning the generator each denoising step buys. For the autoregressive models this is NOT the headline comparison.

| generator | base alone | + our objective | Δ pp | Δ CI | n paired | Δ mean margin | extra wall s | extra reward evals |
|---|---|---|---|---|---|---|---|---|
| AbMPNN (frozen) | 68.8% | 79.8% | +11.0 ** | [+5, +18] | 189 | -0.278 | +89 | +463 |
| ProteinMPNN (frozen) | 51.7% | 62.6% | +10.9 ** | [+5, +18] | 186 | -0.303 | +94 | +502 |
| ESM-IF1 (frozen) | 39.2% | 62.0% | +22.7 ** | [+16, +30] | 189 | -0.617 | +90 | +477 |

## Lift over each arm's OWN edit-distance-matched random control

Beat-native alone confounds design quality with how far the design moved from the native CDR: more edits score worse under the grader regardless of who proposed them. Each arm is therefore compared to a random control that mutates the SAME number of positions per design. LIFT = arm rate minus its own null rate, on the targets where both are graded; it is the part of the score not explained by edit distance.

Restricted to targets where ALL FOUR cells are graded (base, ours, and both nulls), so the two lifts are computed on the SAME targets and are comparable. Per-arm lift on differing subsets is not.

| generator | arm | beat-native | its random null | LIFT pp | Δ LIFT (ours-base) | n |
|---|---|---|---|---|---|---|
| AntiFold | `base_antifold` | 50.1% | 29.8% | **+20.4** | | 190 |
| AntiFold | `ours_antifold` | 64.3% | 21.9% | **+42.4** | **+22.0** ** [+12, +33] | 190 |
| AbMPNN | `base_abmpnn_ar` | 77.1% | 34.0% | **+43.1** | | 190 |
| AbMPNN | `ours_abmpnn_live` | 83.4% | 20.4% | **+62.9** | **+19.8** ** [+9, +30] | 190 |
| ProteinMPNN | `base_proteinmpnn_ar` | 53.5% | 7.9% | **+45.6** | | 190 |
| ProteinMPNN | `ours_proteinmpnn_live` | 62.6% | 7.2% | **+55.4** | **+9.8** ** [+3, +18] | 190 |

## Old-51 subset vs the new targets (is the expansion changing the story?)

| arm | beat-native on old 51 | beat-native on the 149 new | n old | n new |
|---|---|---|---|---|
| `base_abmpnn_ar` | 84.3% | 75.5% | 51 | 139 |
| `base_abmpnn_npz` | 72.5% | 69.1% | 51 | 139 |
| `base_antifold` | 52.9% | 48.2% | 51 | 139 |
| `base_diffab` | 0.0% | nan% | 1 | 0 |
| `base_diffusion` | 0.0% | 2.2% | 51 | 139 |
| `base_esmif` | 45.1% | 39.1% | 51 | 138 |
| `base_esmif_ar` | 43.1% | 43.2% | 51 | 139 |
| `base_proteinmpnn` | 62.7% | 50.7% | 51 | 136 |
| `base_proteinmpnn_ar` | 66.7% | 51.8% | 51 | 139 |
| `base_rfantibody` | 0.0% | 0.0% | 26 | 27 |
| `grad_diffusion` | 39.2% | 44.5% | 51 | 137 |
| `obj_mint` | 100.0% | nan% | 1 | 0 |
| `obj_wt` | 100.0% | nan% | 1 | 0 |
| `ours_abmpnn_live` | 82.4% | 84.9% | 51 | 139 |
| `ours_abmpnn_npz` | 82.4% | 79.0% | 51 | 138 |
| `ours_antifold` | 72.5% | 64.0% | 51 | 139 |
| `ours_diffusion` | 13.7% | 10.1% | 51 | 138 |
| `ours_esmif` | 68.6% | 61.9% | 51 | 139 |
| `ours_proteinmpnn` | 68.6% | 63.8% | 51 | 138 |
| `ours_proteinmpnn_live` | 74.5% | 60.4% | 51 | 139 |
| `rand_base_abmpnn_ar` | 47.1% | 30.9% | 51 | 139 |
| `rand_base_antifold` | 25.5% | 30.2% | 51 | 139 |
| `rand_base_proteinmpnn_ar` | 15.7% | 6.5% | 51 | 139 |
| `rand_ours_abmpnn_live` | 29.4% | 18.0% | 51 | 139 |
| `rand_ours_antifold` | 19.6% | 21.6% | 51 | 139 |
| `rand_ours_proteinmpnn_live` | 5.9% | 7.2% | 51 | 139 |
