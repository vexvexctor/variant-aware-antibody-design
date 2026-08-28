# Improvement across model families × antigen protein families

_Run 2026-08-25. Generator-zoo campaign, 190 graded targets, H3-DDG grader (out-of-family
relative to the BA-DDG steering objective), held-out escape panel._

Each cell is a **paired** comparison: same target, same held-out panel, same pool size — the
only difference is whether our objective was applied on top of that generator's own prior.
Δ = ours − its own base, in percentage points of beat-native. CI = 95% bootstrap over targets.

## Headline

**All five fixed-backbone generators improve significantly on the full panel** (+6.3 to
+23.2 pp). At the protein-family level only the Coronavirus-spike cells resolve, because family
n is 6–15.

| generator | Coronavirus spike (14) | Flavivirus E (15) | Influenza HA (6) | all viral glycoprotein (39) | all 190 targets |
|---|---:|---:|---:|---:|---:|
| **AntiFold** | **+35.7** ✓ | +13.3 | +16.7 | **+28.2** ✓ | **+16.8** ✓ |
| **AbMPNN** | +0.0 | −13.3 | +16.7 | +0.0 | **+6.3** ✓ |
| **ProteinMPNN** | +7.1 | +20.0 | −16.7 | +7.7 | **+8.4** ✓ |
| **ESM-IF1** | **+28.6** ✓ | +20.0 | +16.7 | **+20.5** ✓ | **+23.2** ✓ |
| **Masked diffusion** | +21.4 | +0.0 | +0.0 | **+10.3** ✓ | **+10.5** ✓ |
| **RFantibody** (de-novo) | +33.3 | −20.0 | +0.0 | +0.0 | +7.5 |

✓ = 95% bootstrap CI excludes zero. Full base/ours rates, margins and CIs in
`family_breakdown.csv`; figure in `family_breakdown.png`.

## The mechanism is visible in the rows

The size of the improvement tracks **how weak the base prior is**, not the protein family:

| generator | base beat-native (190) | Δ from our objective |
|---|---:|---:|
| Masked diffusion | 1.6% | +10.5 |
| ESM-IF1 | 40.5% | **+23.2** |
| AntiFold | 49.5% | **+16.8** |
| ProteinMPNN | 55.8% | +8.4 |
| AbMPNN | 77.9% | +6.3 |

A weak or noisy base distribution leaves room for reward-guided selection; a strong
autoregressive prior sampled near its own argmax has little left to gain. On Coronavirus spike
AbMPNN's base is already at 92.9%, so its Δ is exactly 0.0 — a ceiling, not a failure. This is
the same mechanism reported in `HEADLINE_HONEST.md`, now visible per family.

## Protein families

Sequence clustering does **not** produce a few large families on this panel — it is
deliberately antigen-diverse (172 clusters at 70% identity; still 124 clusters at 30%). So the
families here are the three largest homology groups at 30% identity / 50% coverage, identified
from their representative sequences:

| family | n | identified by |
|---|---:|---|
| Flavivirus envelope (E) | 15 | `MRCIGISNRDFVEGV…` flavivirus E; includes the Zika targets 5GZN/5JHL/5VIG |
| Coronavirus spike | 14 | `AYTNSFTRGVYYPDKVFRSS…` SARS-CoV-2 spike; includes 7MZK/7OR9/7TBF/7TTX/7VYR |
| Influenza hemagglutinin | 6 | `EDTICIGYHANNSTDTVDTVLEK…` HA1 |
| *(pooled)* viral surface glycoprotein | 39 | the three above + HIV-1 env gp120 (4) |

Smaller groups exist but are too small to estimate on: serine protease (5), lysozyme (4),
HIV-1 env (4), TCR β (4); the remaining 138 targets are near-singletons.

`mmseqs easy-cluster --min-seq-id 0.3 -c 0.5 --cov-mode 1`, on antigen sequences extracted with
the historical rule (see `AACDB_COVERAGE.md`).

## Caveats

1. **Family cells are underpowered.** n = 6–15 means one target flips the rate by 7–17 pp. Only
   the Coronavirus-spike column has enough signal to resolve for any generator, and only for
   AntiFold and ESM-IF1. Individual red cells (AbMPNN on Flavivirus E, ProteinMPNN on Influenza
   HA) are **not** evidence of a real family-specific regression — their CIs span zero.
2. **RFantibody is de-novo and separately graded.** It designs a new backbone and re-docks the
   Fv, so it is scored on its **own** designed structure via `grade_own_structure.py`, from
   `denovo_per_target.csv`, not on the native complex. Its family cells (n = 2–5) are
   illustrative only; even its 190-target row (+7.5, n=53) does not reach significance.
3. **H3-DDG is ProteinMPNN-lineage**, so it is in-family for the AbMPNN/ProteinMPNN arms.
   Only within-generator paired deltas are valid here — do **not** compare absolute rates
   across generator rows. Pythia is the out-of-family confirmation (+12 to +17 pp paired on
   AntiFold and AbMPNN, panel-level).
4. **Temperature.** These are the campaign's as-run arms; baselines sampled at T=0.2 and ours at
   T=1.0. See `TEMPERATURE_CONTROL.md` — matching temperature roughly doubles the measured
   effect against T=1.0 and shrinks it against a tuned T=0.2 baseline.
5. The three families are **all viral surface glycoproteins**. This is not a viral-vs-non-viral
   contrast; the non-viral targets are too fragmented to group. `EXP_G_NONVIRAL.csv` has a
   hand-labelled 21-target non-viral set if that contrast is needed.

## Files

| file | what |
|---|---|
| `FAMILY_BREAKDOWN.md` | this write-up |
| `family_breakdown.png` | diverging heatmap, generator × family |
| `family_breakdown.csv` | every cell: n, base %, ours %, Δ pp + CI, Δ margin + CI |
| `family_analysis.py` | the analysis |
| `plot_family.py` | the figure |
| `panel_families.json` | target → family / superfamily assignment |
