# Temperature-matched comparison — FINAL, n=190 targets / 172 antigen clusters

## The confound
Baselines sampled CDR-H3 at `--sampling-temp 0.2` (near-argmax); SVDD draws its candidates from
the generator at **T = 1.0**. So "generator alone vs + our objective" was also, silently,
"T=0.2 vs T=1.0". The `*_t1` arms are the identical generator sampled at T=1.0.

## Result (H3-DDG, held-out escape split, cluster-bootstrap CI)

| generator | base T=0.2 | base T=1.0 | + our objective | delta vs T=1.0 | delta vs T=0.2 |
|---|---|---|---|---|---|
| AbMPNN | 66.6% | 41.7% | 66.4% | **+24.7 ** [+17,+32]** | **-0.2** |
| AntiFold | 41.3% | 32.8% | 51.7% | **+18.8 ** [+13,+25]** | **+10.4 *** |
| ProteinMPNN | 47.7% | 24.2% | 41.8% | **+17.6 ** [+11,+25]** | **-5.9** |
| ESM-IF1 | — | — | — | insufficient overlap | — |

## Reading — BOTH columns matter and they say different things

**vs T=1.0 (mechanistic):** the objective adds +17.6 to +24.7pp over the distribution it
actually samples from. All significant at full panel size.

**vs T=0.2 (deployment-realistic):** AntiFold +10.4** , AbMPNN -0.2, ProteinMPNN -5.9.
Against a well-tuned low-temperature baseline the objective **beats AntiFold, ties AbMPNN and
slightly trails ProteinMPNN**.

## The honest single sentence
The objective recovers essentially all the quality that low-temperature sampling buys — from a
much noisier starting distribution — while additionally optimising escape robustness; but on
raw beat-native it does not clearly beat a well-tuned T=0.2 baseline except on AntiFold.

Quoting only "+24.7" would overstate the result. Quoting only "-0.2" would hide that the
baseline's advantage is a sampling-temperature effect the objective fully compensates for.

## Open follow-up (RUNNING, job 73135)
`ours_*_sharp`: SVDD with the base sharpened to T=0.2 before candidate sampling. Tests whether
the objective STACKS on low-temperature sampling or merely substitutes for it. If it lands
above the T=0.2 baselines, the deployment claim holds on all three generators.
