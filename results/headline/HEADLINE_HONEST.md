# Where the headline actually stands (2026-08-19, end of day)

The claim narrowed over the course of the campaign. Every narrowing came from a CONTROL we
added, not from new data. Recording the sequence because the intermediate numbers were all
circulated and each looked solid at the time.

| stage | claim | what killed / changed it |
|---|---|---|
| initial | AbMPNN -9.3pp, "significant" | **rep-mismatched pools** — baselines had 3 graded reps (~18 designs), steered arms ~1 (~6). best-of-pool is monotone in pool size. |
| rep-matched | +7 to +18pp, significant on 5 generators | **edit-distance confound** — the grader penalises distance-from-native monotonically and our objective adds ~1.4 edits |
| + lift-over-null | +8 to +22pp, significant | **temperature confound** — baselines sampled T=0.2, SVDD samples T=1.0 |
| + temperature-matched | +17.6 to +24.7 vs T=1.0 (significant); **+10.4 / -0.2 / -5.9 vs T=0.2** | matched comparison is the mechanistic claim; the T=0.2 column is the deployment claim |
| + base-temp sweep, matched targets | see below | **subset artifact** — deltas computed over different target subsets are not comparable |

## Final, on matched targets

| generator | best config | delta vs its OWN tuned (T=0.2) baseline | n |
|---|---|---|---|
| **AntiFold** | raw prior (T=1.0) | **+12.5 ** [+3,+22]** | 93 |
| AbMPNN | flat across T=0.1/0.2/0.5/1.0 | **-2.0 .. +0.0** (no benefit at any temperature) | 73 |
| ProteinMPNN | undetermined | insufficient shared targets | 23 |

## The honest claim
The objective **clearly helps AntiFold** (+12.5pp over a tuned baseline, significant) and shows
**no measurable benefit on AbMPNN at any base temperature**. The plausible mechanism: AntiFold's
prior is weak and already sharp, leaving room for reward-guided selection; a strong
autoregressive prior sampled near-argmax is already close to its own optimum, so selection has
nothing to add.

What remains solidly true regardless:
- vs the distribution SVDD actually samples from (T=1.0): +17.6 to +24.7pp, significant, n=190.
  The objective genuinely compensates for a noisy base distribution.
- Out-of-family (Pythia) paired deltas: +12 to +17pp, significant on AntiFold and AbMPNN.
- Zero training cost; ~15s and ~77 reward evals per design.
- Panel expanded 51 -> 200 targets, 30 -> 172 antigen clusters; every baseline degrades on the
  harder new targets while the steered arm holds flat.

## WHY the negatives may be partly a BENCHMARK limitation (the most useful thing learned today)

`WT_ESCAPE_COUPLING.txt`: held-out escape worst-case correlates with WT binding at
**r = +0.720 pooled (r^2 = 0.52)**, median within-target r = +0.659, and on **86/190 targets
(45%) r > 0.7**.

So roughly half the variance in "survives the escape panel" is simply "binds the wild-type
antigen tightly". An escape-SPECIFIC objective can only compete for the remainder -- which is
exactly why `ours` ties `WT-affinity only` (-1.4pp, n.s.) and why the gaps between objectives
and between generators are compressed. It also predicts the temperature result: low-temperature
sampling yields high-affinity designs, which then score well on escape BECAUSE of the coupling.

The panel is 86% triple mutants proposed by ProteinMPNN around the same epitope and filtered by
ESM-1v plausibility -- a process that tends to produce variants whose binding tracks the WT
interaction. The benchmark does not isolate escape.

**Fix (cheap, re-slices existing data, no new generation):**
1. Single-mutant-enriched eval panel — 7,058 singles already in the library; single
   substitutions are where real escape decouples from affinity.
2. Decoupled-variant selection — stratify the held-out panel on LOW WT-binding correlation
   instead of on omega. Least-coupled targets are listed in WT_ESCAPE_COUPLING.txt.

This should be done BEFORE adding more graders: FoldX/MM-GBSA add measurement precision to a
quantity that is confounded upstream.

## What would strengthen it
1. Finish the sweep (T=0.5/0.7 and ProteinMPNN overlap) — expect flat, not peaked.
2. MM-GBSA as the second out-of-family grader (running).
3. A single-mutant-enriched eval panel: the held-out set is 86% triple mutants and only
   47/200 targets contain any single mutant, which is not the clinically relevant regime.
