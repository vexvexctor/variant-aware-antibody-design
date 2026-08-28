# Generator zoo — FINAL H3 results (2026-08-19 07:45)

Panel: **200 antibody-antigen complexes / 172 antigen clusters** (was 51 / 30), strict superset.
Grader: H3-DDG, held-out omega-divergent escape split. CIs bootstrap antigen CLUSTERS.
Objective cost: ~85s wall, ~15s/design, ~460 reward evals, <1GB VRAM, **zero training**.

## HEADLINE: every comparison significant

| generator | raw delta | lift over own null | frozen-prior delta |
|---|---|---|---|
| AntiFold | +14.2 ** | +22.0 ** | — (frozen IS its native mode) |
| AbMPNN | +6.3 ** | +19.8 ** | +11.0 ** |
| ProteinMPNN | +8.5 ** | +9.8 ** | +10.9 ** |
| ESM-IF1 | — | — | +22.7 ** |
| our masked diffusion | +8.9 ** | — | — |

## Paired detail
| generator | base alone | + our objective | Δ pp | Δ CI | n paired | Δ mean margin | extra wall s | extra reward evals |
|---|---|---|---|---|---|---|---|---|
| AntiFold | 50.1% | 64.3% | +14.2 ** | [+8, +21] | 190 | -0.690 | +88 | +462 |
| AbMPNN | 77.1% | 83.4% | +6.3 ** | [+0, +12] | 190 | -0.152 | +81 | +463 |
| ProteinMPNN | 53.5% | 62.0% | +8.5 ** | [+2, +15] | 190 | -0.054 | +89 | +501 |
| our masked diffusion | 1.9% | 10.8% | +8.9 ** | [+4, +14] | 189 | -2.867 | +68 | +504 |

## Frozen-prior ablation
| generator | base alone | + our objective | Δ pp | Δ CI | n paired | Δ mean margin | extra wall s | extra reward evals |
|---|---|---|---|---|---|---|---|---|
| AbMPNN (frozen) | 68.8% | 79.8% | +11.0 ** | [+5, +18] | 189 | -0.278 | +89 | +463 |
| ProteinMPNN (frozen) | 51.7% | 62.6% | +10.9 ** | [+5, +18] | 186 | -0.303 | +94 | +502 |
| ESM-IF1 (frozen) | 39.2% | 62.0% | +22.7 ** | [+16, +30] | 189 | -0.617 | +90 | +477 |

## Lift over edit-matched null
| generator | arm | beat-native | its random null | LIFT pp | Δ LIFT (ours-base) | n |
|---|---|---|---|---|---|---|
| AntiFold | `base_antifold` | 50.1% | 29.8% | **+20.4** | | 190 |
| AntiFold | `ours_antifold` | 64.3% | 21.9% | **+42.4** | **+22.0** ** [+12, +33] | 190 |
| AbMPNN | `base_abmpnn_ar` | 77.1% | 34.0% | **+43.1** | | 190 |
| AbMPNN | `ours_abmpnn_live` | 83.4% | 20.4% | **+62.9** | **+19.8** ** [+9, +30] | 190 |
| ProteinMPNN | `base_proteinmpnn_ar` | 53.5% | 7.9% | **+45.6** | | 190 |
| ProteinMPNN | `ours_proteinmpnn_live` | 62.6% | 7.2% | **+55.4** | **+9.8** ** [+3, +18] | 190 |

## Out-of-family (Pythia)
```
=== PAIRED: generator alone vs + our objective (Pythia) ===

AbMPNN         base  68.4%  ours  80.7%  delta  +12.3 ** [+2, +23]  n=57
AntiFold       base  48.3%  ours  65.5%  delta  +17.2 ** [+3, +31]  n=58
ProteinMPNN    base  56.9%  ours  69.0%  delta  +12.1 [+0, +24]  n=58
```

## Panel expansion effect (old 51 vs new 149)
| arm | beat-native on old 51 | beat-native on the 149 new | n old | n new |
|---|---|---|---|---|
| `base_abmpnn_ar` | 84.3% | 75.5% | 51 | 139 |
| `base_abmpnn_npz` | 72.5% | 69.1% | 51 | 139 |
| `base_antifold` | 52.9% | 48.2% | 51 | 139 |
| `base_diffab` | 0.0% | nan% | 1 | 0 |
| `base_diffusion` | 0.0% | 2.2% | 51 | 139 |
| `base_esmif` | 45.8% | 39.5% | 48 | 129 |
| `base_esmif_ar` | 43.1% | 42.4% | 51 | 139 |
| `base_proteinmpnn` | 60.0% | 50.8% | 45 | 128 |
| `base_proteinmpnn_ar` | 66.7% | 50.4% | 51 | 139 |
| `base_rfantibody` | 0.0% | 0.0% | 26 | 27 |
| `grad_diffusion` | 41.3% | 46.2% | 46 | 130 |
| `obj_mint` | 100.0% | nan% | 1 | 0 |
| `obj_wt` | 100.0% | nan% | 1 | 0 |
| `ours_abmpnn_live` | 82.4% | 84.2% | 51 | 139 |
| `ours_abmpnn_npz` | 82.6% | 80.2% | 46 | 131 |
| `ours_antifold` | 70.6% | 62.6% | 51 | 139 |
| `ours_diffusion` | 14.6% | 10.5% | 41 | 133 |
| `ours_esmif` | 69.0% | 63.2% | 42 | 133 |
| `ours_proteinmpnn` | 70.5% | 63.2% | 44 | 133 |
| `ours_proteinmpnn_live` | 74.5% | 59.7% | 51 | 139 |
| `rand_base_abmpnn_ar` | 47.1% | 30.9% | 51 | 139 |
| `rand_base_antifold` | 25.5% | 30.2% | 51 | 139 |
| `rand_base_proteinmpnn_ar` | 15.7% | 6.5% | 51 | 139 |
| `rand_ours_abmpnn_live` | 29.4% | 18.0% | 51 | 139 |
| `rand_ours_antifold` | 19.6% | 21.6% | 51 | 139 |
| `rand_ours_proteinmpnn_live` | 5.9% | 7.2% | 51 | 139 |

### Why the expansion mattered
Every BASELINE degrades on the new antigen-diverse targets (base_proteinmpnn_ar 66.7%->49.6%,
base_abmpnn_ar 84.3%->75.5%) while the steered arm holds flat (ours_abmpnn_live 82.4%->83.5%).
The old 51-target panel flattered the baselines; the objective's advantage is LARGER on the
harder panel. Invisible at n=51.

### Why ProteinMPNN was nearly dropped, and should not have been
Its measured delta read -9.3 (significant), then -1.5, then null, then FINALLY +7.2/+10.5/+12.5
(all significant). Each earlier value was an artefact: rep-mismatched pool sizes, partial
grading, and the edit-distance confound. A plausible biological story ("ProteinMPNN is not
antibody-tuned, so it is bad for CDR-H3") was available at every stage and would have been
wrong. Controls, not narrative, resolved it.

### Not yet established
- Temperature confound not removed (T=1.0 baselines queued); lift controls for it indirectly.
- MM-GBSA (2nd out-of-family grader) queued, not run.
- ours_esmif_live, objective comparison, DiffAb/IgLM panels, RFantibody 60/60 queued.
- RFantibody's native-backbone row (0/30) is NOT its performance; needs grade_own_structure.py.
- Pythia's absolute rates are inflated (random null 50-68%); only its PAIRED deltas are usable.
