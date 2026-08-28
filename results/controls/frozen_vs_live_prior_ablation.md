### Ablation: same models with a FROZEN prior on both sides

Isolates what re-conditioning the generator each denoising step buys. For the autoregressive models this is NOT the headline comparison.

| generator | base alone | + our objective | Δ pp | Δ CI | n paired | Δ mean margin | extra wall s | extra reward evals |
|---|---|---|---|---|---|---|---|---|
| AbMPNN (frozen) | 69.0% | 79.9% | +10.9 ** | [+5, +17] | 190 | -0.287 | +89 | +463 |
| ProteinMPNN (frozen) | 51.6% | 63.5% | +11.9 ** | [+6, +19] | 190 | -0.314 | +94 | +502 |
| ESM-IF1 (frozen) | 39.0% | 61.6% | +22.6 ** | [+16, +30] | 190 | -0.619 | +90 | +478 |

## Lift over each arm's OWN edit-distance-matched random control
