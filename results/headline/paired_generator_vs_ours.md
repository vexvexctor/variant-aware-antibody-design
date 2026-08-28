## Paired: generator + our objective vs generator alone

Same base distribution, same pool size, same panel — the only difference is the objective. Positive pp = our objective helps that generator.

| generator | base alone | + our objective | Δ pp | Δ CI | n paired | Δ mean margin | extra wall s | extra reward evals |
|---|---|---|---|---|---|---|---|---|
| AntiFold | 50.1% | 64.3% | +14.2 ** | [+8, +21] | 190 | -0.689 | +88 | +462 |
| AbMPNN | 77.1% | 83.4% | +6.3 ** | [+0, +12] | 190 | -0.152 | +81 | +463 |
| ProteinMPNN | 53.5% | 62.6% | +9.1 ** | [+3, +16] | 190 | -0.054 | +89 | +501 |
| our masked diffusion | 1.9% | 11.9% | +10.1 ** | [+5, +15] | 190 | -2.729 | +68 | +505 |

### Ablation: same models with a FROZEN prior on both sides
