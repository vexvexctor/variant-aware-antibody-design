# Do our evaluators agree? Score-level and decision-level

_Run 2026-08-27._

**One shared frozen candidate pool: 1,785 identical (target, seed, selector, design) units,
every evaluator scoring exactly the same designs.** No unit is scored by one evaluator and not
another, so every pair below is exactly matched — unlike the earlier partial matrix in
`CLAIM_AUDIT.md`, which pooled across arms with differing coverage (884 pairs for BA↔H3 there
vs 1,785 here; that is why BA↔H3 reads +0.575 there and +0.661 here).

Margin `w` = held-out worst-case design-minus-native; **`w < 0` = "beats native"** is the
decision every headline rests on.

## Score agreement — Spearman of the margins

|  | BA-DDG | H3 | Pythia | FoldX | MM-GBSA |
|---|---:|---:|---:|---:|---:|
| **BA-DDG** | — | **0.661** | 0.283 | 0.044 | 0.121 |
| **H3** | 0.661 | — | 0.469 | 0.207 | 0.246 |
| **Pythia** | 0.283 | 0.469 | — | 0.475 | 0.345 |
| **FoldX** | 0.044 | 0.207 | 0.475 | — | 0.294 |
| **MM-GBSA** | 0.121 | 0.246 | 0.345 | 0.294 | — |

## Decision agreement — Cohen's κ on the beat-native call

**This is the number to quote.** Raw percent agreement is inflated by base rates and is
misleading here (see below).

|  | BA-DDG | H3 | Pythia | FoldX | MM-GBSA |
|---|---:|---:|---:|---:|---:|
| **BA-DDG** | — | **0.429** | 0.170 | 0.067 | 0.004 |
| **H3** | 0.429 | — | **0.325** | 0.081 | −0.000 |
| **Pythia** | 0.170 | 0.325 | — | **0.317** | 0.023 |
| **FoldX** | 0.067 | 0.081 | 0.317 | — | −0.012 |
| **MM-GBSA** | 0.004 | −0.000 | 0.023 | −0.012 | — |

Landis–Koch bands: <0.20 slight, 0.21–0.40 fair, 0.41–0.60 moderate.

## Raw percent agreement — and why not to use it

|  | BA-DDG | H3 | Pythia | FoldX | MM-GBSA |
|---|---:|---:|---:|---:|---:|
| **BA-DDG** | — | 73.4% | 62.1% | 63.7% | 66.3% |
| **H3** | 73.4% | — | 68.1% | 60.4% | 60.1% |
| **Pythia** | 62.1% | 68.1% | — | 71.8% | 63.5% |
| **FoldX** | 63.7% | 60.4% | 71.8% | — | **82.4%** |
| **MM-GBSA** | 66.3% | 60.1% | 63.5% | 82.4% | — |

**FoldX ↔ MM-GBSA is the trap.** It is the highest percent agreement in the table (82.4%) and
the *lowest* κ (−0.012). They agree because they both almost always say "no", not because they
see the same thing. Base rates on the identical 1,785 designs:

| evaluator | beat-native rate |
|---|---:|
| H3 | 39.7% |
| Pythia | 37.0% |
| BA-DDG | 33.4% |
| FoldX | 16.5% |
| MM-GBSA | **1.3%** |

## Reading

1. **BA-DDG ↔ H3 is the only moderate pair** (ρ=0.661, κ=0.429). Both import ProteinMPNN; they
   are the same family. This is why H3 must be described as a **held-out evaluator, not an
   independent one**, and why cross-generator absolute rates on H3 are invalid for
   ProteinMPNN/AbMPNN-derived arms.
2. **Pythia is the genuine bridge.** It is the only evaluator with fair agreement to *both*
   families — H3 (κ=0.325) and FoldX (κ=0.317) — while sitting outside the ProteinMPNN lineage.
   That is what makes it the out-of-family confirmation of choice.
3. **BA-DDG ↔ FoldX is essentially unrelated** (ρ=0.044, κ=0.067). Near-zero correlation is
   *not* evidence of clean orthogonality: combined with the random-edit control (random CDR-H3
   edits score 76.5% under FoldX), it is evidence FoldX **barely discriminates** on this task.
4. **MM-GBSA agrees with nothing** (κ from −0.012 to 0.023 against all four). At a 1.3%
   beat-native rate it is close to a constant "no", which is why the effect does not replicate
   on it. Its disagreement is not independent confirmation of a negative — it is a
   near-degenerate rater.
5. **The evaluators are not redundant.** Outside the BA/H3 pair, chance-corrected agreement is
   fair-to-slight everywhere. A claim resting on any single evaluator is weakly supported;
   the defensible pattern is a within-generator paired delta that holds on H3 **and** Pythia.

## For the paper

Put the κ matrix next to the central table. It pre-empts the reviewer who computes it
themselves, converts the "your graders are entangled" objection into a demonstration of rigour,
and it is the evidence for the sentence *"H3-DDG is held out from optimization but is not
independent of BA-DDG; Pythia provides the out-of-family check."*

## Files

| file | what |
|---|---|
| `GRADER_AGREEMENT.md` | this write-up |
| `grader_agreement.png` | the two matrices side by side |
| `grader_agreement.csv` | every pair: n, Spearman, % agreement, κ, per-evaluator beat-native rates |
| `grader_agreement.py` | the computation (stdlib only) |
| `plot_grader_matrix.py` | the figure |

Source: `scratch/sp1grade/central_table_per_unit.csv` (8,925 rows = 1,785 units × 5 evaluators).
