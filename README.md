# Variant-aware CDR-H3 Design

Antibody sequence generators normally condition on a single fixed antigen structure. If the
antigen then varies — as viral surface proteins routinely do — nothing in the generation
procedure has asked whether the designed antibody still binds. This work asks whether exposing a
**frozen** generator to several plausible antigen variants *during inference* can favour
antibodies whose predicted binding persists across antigen variation, without retraining the
generator or touching its weights.


## Method in brief

We build **RobustBench**, 200 antibody–antigen complexes each carrying 30 model-generated
antigen variants split into **6 design-time** and **24 held-out**. During CDR-H3 decoding a
frozen generator proposes candidate residues; each candidate completion is scored against the
6 design variants, the resulting design-minus-native margins are aggregated with
**CVaR at α = 0.2** over the worst tail, and that reward reweights the generator's own proposal
distribution before a residue is committed. The generator is never updated and no gradients are
taken through it. Designs are then scored on the **24 held-out variants**, which no part of the
design procedure ever saw.

```text
for step in CDR-H3 decoding order:
    props   = frozen_generator.propose(prefix)          # generator is never updated
    for each candidate residue c in top-k props:        # k = 6
        x_c      = complete_cdr(prefix + c)
        margins  = [score(x_c, v) - score(native, v) for v in design_variants]   # 6 variants
        R[c]     = -CVaR(margins, alpha=0.2)            # worst tail; lower margin is better
    commit(sample(props reweighted by R))               # derivative-free
evaluate(final_design, heldout_variants)                # 24 variants, never seen
```

## Key results

Paired, same base distribution, same pool size, same held-out panel; H3-DDG evaluator (held out
from optimization); 190 graded targets; intervals bootstrap antigen clusters.

| generator | unguided | + variant-aware | Δ pp |
|---|---:|---:|---:|
| ESM-IF1 | 39.2% | 62.0% | +22.7 [+16, +30] |
| AntiFold | 50.1% | 64.3% | +14.2 [+8, +21] |
| ProteinMPNN | 53.5% | 62.0% | +8.5 [+2, +15] |
| AbMPNN | 77.1% | 83.4% | +6.3 [+0, +12] |
| masked diffusion | 1.9% | 10.8% | +8.9 [+4, +14] |

Reading these correctly requires four things, all documented in
`results/headline/HEADLINE_HONEST.md`: replicate counts must be matched (best-of-pool is
monotone in pool size); the evaluator penalises distance from native, so effects are also
reported against an edit-matched null; baselines sampled at T = 0.2 while steered runs sampled
at T = 1.0, and both framings are reported; and held-out escape correlates with wild-type
binding at r = +0.72, so this benchmark does not cleanly separate an escape-specific objective
from an affinity one.

Evaluators are **not** mutually independent — see §15.

## Repository structure

```
configs/      as-run settings: benchmark, steering, generators, evaluators
data/         manifests, the 6/24 splits, sequences  (no raw structures)
src/
  benchmark/    RobustBench construction: extraction, clustering, panels
  generators/   frozen generator wrappers
  steering/     the objective (cvar.py) and the decoding loop (variant_aware.py)
  evaluators/   BA-DDG, H3-DDG, Pythia-PPI, FoldX, MM-GBSA
  analysis/     metrics, cluster bootstrap, controls, agreement, plots
scripts/      00-08 numbered pipeline + reproduce_all.sh
results/      headline, controls, generalization, evaluator_agreement, trajectories, biology, compute
figures/      regenerated from results/ by scripts/08_make_figures.sh
tests/        CVaR sign/tail, split disjointness, metric conventions
paper/        manuscript (see paper/README.md)
```

## Installation

```bash
git clone <repo-url> && cd variant-aware-antibody-design
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 -m unittest discover -s tests      # 27 tests, no GPU or data needed
```

The analysis and figure path needs only numpy/pandas/scipy/matplotlib/PyYAML.

## Downloading models and data

No third-party weights, structures or binaries are redistributed here.

```bash
export VAAD_ROOT=/path/to/data-root
export VAAD_TOOLS=$VAAD_ROOT/tools
bash scripts/00_download_dependencies.sh    # prints the required layout and sources
```

You need AACDB structures, ProteinMPNN, AntiFold, AbMPNN, ESM-IF1/ESM-1v, Pythia-PPI,
RFantibody and mmseqs2. **FoldX is proprietary** — obtain your own academic licence. Nothing in
this repository hard-codes a machine path; everything resolves through `VAAD_ROOT`
(`src/utils/paths.py`).

## Quick start

Rebuild every paper statistic and figure from the released result files. No GPU, no downloads,
about a minute:

```bash
bash scripts/reproduce_all.sh --from-precomputed
```

To design one target end to end (needs the staged dependencies above and a GPU):

```bash
bash scripts/03_run_variant_aware_design.sh 1IC7_HLY
```

## Reproducing RobustBench

```bash
bash scripts/01_build_robustbench.sh
```

Extracts antigen sequences (4.5 Å epitope / 8.0 Å shell contact rules), clusters with mmseqs2 at
`--min-seq-id 0.7 -c 0.8 --cov-mode 0 --cluster-mode 0`, scans structurally eligible candidates,
selects the 200 targets and builds the 6/24 panels. Released manifests are already in `data/`,
so this is only needed to rebuild from scratch. See `data/README.md` for the **171 vs 172**
cluster-count distinction — both numbers are correct and they answer different questions.

## Running unguided generation

```bash
bash scripts/02_generate_baselines.sh
```

Each frozen generator alone, no reward. Per-generator settings are in `configs/generators.yaml`.

## Running variant-aware generation

```bash
bash scripts/03_run_variant_aware_design.sh <TARGET_ID>
```

Identical generator and sampling settings, with the CVaR objective active during decoding. The
objective is `src/steering/cvar.py`; the decoding loop is `src/steering/variant_aware.py`.

## Held-out evaluation

```bash
bash scripts/04_evaluate_heldout.sh <manifest.tsv> <design_dir> <out_dir>
```

Scores each design on its 24 held-out variants. A design **beats native** only when its *worst*
held-out margin is below zero. Rates are paired, intersected across conditions, and
bootstrapped over antigen clusters using the fixed benchmark labels.

## Reproducing paper results

```bash
bash scripts/reproduce_all.sh --from-precomputed   # statistics + figures from released CSVs
bash scripts/reproduce_all.sh --full               # full generation + scoring; multi-GPU, days
```

`REPRODUCIBILITY.md` maps every manuscript claim to its script, input and output, and states
plainly which artifacts are missing from this archive.

## Controls and ablations

```bash
bash scripts/05_run_controls.sh
```

Temperature-matched control, edit-distance-matched null, objective ablation (CVaR vs mean vs
worst-case), design-panel-size sensitivity (M = 1/3/6/12), budget-matched generate-and-rank, and
RFantibody kept separate because it generates its own backbone. Results in `results/controls/`.

Two of these do not favour the method and are reported as found: a budget-matched genetic
algorithm beats the steered arm on the in-lineage evaluator, and beat-native is flat in M over
1–12, so six design variants is a pre-registered default rather than a tuned optimum.

## Evaluator agreement

```bash
cd src && python3 -m analysis.evaluator_agreement \
    --pool ../data/processed/shared_1785_design_pool.csv \
    --outdir ../results/evaluator_agreement
```

On one shared pool of **1,785 identical designs** scored by all five evaluators:

| | BA-DDG | H3 | Pythia | FoldX | MM-GBSA |
|---|---:|---:|---:|---:|---:|
| **Spearman / κ vs BA-DDG** | — | .661 / .429 | .283 / .170 | .044 / .067 | .121 / .004 |

Full matrices in `results/evaluator_agreement/`. Three consequences:

- **BA-DDG and H3-DDG share ProteinMPNN lineage** (κ = 0.429). H3-DDG is a *held-out* evaluator,
  not an independent one; use within-generator paired deltas, not cross-generator absolute rates.
- **Pythia-PPI is the cleaner out-of-family check** — the only evaluator with fair agreement to
  both families while sharing lineage with neither.
- **MM-GBSA calls beat-native on 1.3% of the pool** and has κ ≈ 0 against everything, behaving
  close to a constant-negative classifier here. Failure to replicate under MM-GBSA is therefore
  weak evidence, not strong independent evidence against an effect.

Raw percent agreement is *not* the headline statistic: FoldX and MM-GBSA agree on 82.4% of calls
at κ = −0.012, because they mostly agree on saying "no".

## Data availability

`data/` ships the benchmark definition — target manifest, antigen/VH sequences, cluster labels
and the 6/24 splits (6,200 variant rows). The full candidate library behind it (≈3.04 M variants,
1.57 GB) is not in git; `data/README.md` documents the archive and the code that regenerates it.
Structures come from AACDB and are not redistributed.

## Limitations

**All evidence here is computational.** Every number is produced by a learned or empirical
scoring function, not by a binding assay, and none of it substitutes for experimental
measurement. Specifically: the antigen variants are model-generated rather than measured, and
86% of the held-out set are triple mutants; agreement between our panel and measured
deep-mutational-scanning escape is moderate (Spearman ≈ 0.33); held-out performance correlates
with wild-type binding at r = +0.72, so the benchmark cannot fully isolate variant-specific
robustness; no result measures escape of the designed antibodies themselves, only their scores
against variants measured on other antibodies; and the evaluators are correlated, so agreement
among them is not independent confirmation. Wet-lab validation is required before any claim
about real binding.

## Citation

See `CITATION.cff`. Cite AACDB alongside this work, and cite the upstream generators and
evaluators you use.
