# Exp1 — Budget-matched AntiFold search baselines: FROZEN, out-of-family (extension)

Extends `scratch/budget_baselines/BUDGET_MATCHED_SUMMARY.md` (the pooled best-of-8 core run: svdd/bestofn/random/greedy/beam/ga all matched at B≈379 BA-DDG reward-evals, same AntiFold base, same reward, same held-out omega eval split). This doc adds the four gaps: (1) three more out-of-family graders, (2) best-of-N temperature sweep, (3) anytime curve, (4) >=5 seeds.

**FROZEN primary (per Exp2).** Every table below selects ONE design per (arm,target,rep) — the top-1 by BA-DDG CVaR20 on the OPT/steer variants — and grades THAT single frozen design out-of-loop on the held-out eval variants. No oracle-selection over a pool (the pooled best-of-8 numbers in BUDGET_MATCHED_SUMMARY.md are the selection-inflated upper bound).

Beat-native margin<0 = frozen design's worst held-out variant strictly tighter than the native antibody's (for cofold, worst-case ipTM exceeds native's). Each (arm,target,rep) is one independent frozen sample; rate = fraction of samples beating native.

## 1. Frozen beat-native, all six graders

| arm | BA-DDG<br>(in-loop) | H3<br>(in-family (MPNN)) | FoldX<br>(OUT-of-family) | MMGBSA<br>(OUT-of-family) | Pythia<br>(OUT-of-family) | cofold<br>(OUT-of-family) |
|---|---|---|---|---|---|---|
| SVDD (value-guided) | 30.2% (n=255, +0.34) | 41.2% (n=255, +0.25) | 55.9% (n=102, +0.05) | 33.3% (n=102, +129.22) | 39.2% (n=102, +1.79) | 49.0% (n=51, +0.01) |
| best-of-N | 45.9% (n=255, -0.17) | 51.4% (n=255, -0.01) | 57.8% (n=102, -0.21) | 47.1% (n=102, +12.71) | 44.1% (n=102, +0.90) | 43.1% (n=51, +0.01) |
| random-search | 15.7% (n=255, +1.68) | 5.1% (n=255, +2.85) | 38.0% (n=100, +2.30) | 32.4% (n=102, +130241.71) | 18.6% (n=102, +6.24) | 58.8% (n=51, +0.00) |
| greedy asc. | 40.4% (n=255, -0.09) | 52.2% (n=255, +0.00) | 48.0% (n=102, +0.61) | 31.4% (n=102, +82.54) | 40.2% (n=102, +1.47) | 51.0% (n=51, -0.00) |
| beam search | 58.8% (n=255, -1.00) | 49.4% (n=255, -0.15) | 43.6% (n=94, +0.91) | 35.3% (n=102, -7.27) | 44.1% (n=102, +1.21) | 52.9% (n=51, +0.01) |
| genetic algorithm | 62.4% (n=255, -1.05) | 61.6% (n=255, -0.52) | 53.9% (n=102, +1.04) | 37.2% (n=102, +4.31) | 40.2% (n=102, +0.86) | 58.8% (n=51, -0.00) |

BA-DDG/H3 use all available reps; FoldX reps 1-2; MMGBSA/Pythia/cofold = frozen top1 of reps 1-2. **FINAL — all grader arrays complete (MMGBSA 306/306, Pythia 306/306, cofold 306/306).**

**CAVEAT — MMGBSA mean-margin column is polluted by outliers** (e.g. random-search +130241, SVDD +129): occasional MM-GBSA clash/unconverged energies blow up the implicit-solvent ΔG. The beat-native *rate* is rank-based and robust; **treat the MMGBSA rate as valid but ignore its mean-margin annotation** (use median or winsorize to report a margin). FoldX/H3/BA-DDG/cofold margins are clean.

### Headline question — does SVDD's edge live in the OUT-of-family arbiters?
- **BA-DDG** (in-loop): SVDD 30.2% (mean +0.34) vs best competitor genetic algorithm 62.4% (-1.05) → SVDD trails.
- **H3** (in-family (MPNN)): SVDD 41.2% (mean +0.25) vs best competitor genetic algorithm 61.6% (-0.52) → SVDD trails.
- **FoldX** (OUT-of-family): SVDD 55.9% (mean +0.05) vs best competitor best-of-N 57.8% (-0.21) → SVDD trails.
- **MMGBSA** (OUT-of-family): SVDD 33.3% (mean +129.22) vs best competitor best-of-N 47.1% (+12.71) → SVDD trails.
- **Pythia** (OUT-of-family): SVDD 39.2% (mean +1.79) vs best competitor best-of-N 44.1% (+0.90) → SVDD trails.
- **cofold** (OUT-of-family): SVDD 49.0% (mean +0.01) vs best competitor genetic algorithm 58.8% (-0.00) → SVDD trails.

## 2. Seed robustness (>=5 reps)

| grader | arm | n_seeds | mean | best | worst | std | pooled | frac seeds >50% |
|---|---|---|---|---|---|---|---|---|
| BA-DDG | SVDD (value-guided) | 5 | 30.2% | 35.3% | 25.5% | 0.032 | 30.2% | 0% |
| BA-DDG | best-of-N | 5 | 45.9% | 51.0% | 39.2% | 0.040 | 45.9% | 20% |
| BA-DDG | random-search | 5 | 15.7% | 17.6% | 13.7% | 0.012 | 15.7% | 0% |
| BA-DDG | greedy asc. | 5 | 40.4% | 45.1% | 37.2% | 0.029 | 40.4% | 0% |
| BA-DDG | beam search | 5 | 58.8% | 68.6% | 51.0% | 0.073 | 58.8% | 100% |
| BA-DDG | genetic algorithm | 5 | 62.4% | 66.7% | 54.9% | 0.042 | 62.4% | 100% |
| H3 | SVDD (value-guided) | 5 | 41.2% | 43.1% | 37.2% | 0.021 | 41.2% | 0% |
| H3 | best-of-N | 5 | 51.4% | 56.9% | 39.2% | 0.062 | 51.4% | 80% |
| H3 | random-search | 5 | 5.1% | 5.9% | 3.9% | 0.010 | 5.1% | 0% |
| H3 | greedy asc. | 5 | 52.2% | 56.9% | 47.1% | 0.032 | 52.2% | 80% |
| H3 | beam search | 5 | 49.4% | 54.9% | 45.1% | 0.038 | 49.4% | 40% |
| H3 | genetic algorithm | 5 | 61.6% | 70.6% | 54.9% | 0.056 | 61.6% | 100% |

## 3. Best-of-N AntiFold sampling-temperature sweep

temp=1.0 == the existing default best-of-N arm; 0.1/0.5 are the new prespecified grid points (temp NOT tuned on eval variants).

| grader | temp | n | beat-native | mean margin |
|---|---|---|---|---|
| BA-DDG | 1.0 | 255 | 45.9% | -0.165 |
| BA-DDG | 0.1 | 102 | 18.6% | +1.036 |
| BA-DDG | 0.5 | 102 | 39.2% | +0.121 |
| H3 | 1.0 | 255 | 51.4% | -0.008 |
| H3 | 0.1 | 102 | 36.3% | +0.580 |
| H3 | 0.5 | 102 | 52.0% | -0.157 |

## 4. Anytime curve

Best-so-far incumbent's HELD-OUT BA-DDG CVaR20 (higher=tighter) at budget = 10/25/50/75/100% of matched B, mean over targets. Figure: `results/figures/exp1_anytime.png`. H3-graded incumbents (stricter, out-of-family) in `h3_any_out/` once job completes.

(n=51 targets)

| arm | 10% | 25% | 50% | 75% | 100% |
|---|---|---|---|---|---|
| SVDD (value-guided) | -0.16 | +0.01 | +0.31 | +0.36 | +0.52 |
| best-of-N | -0.01 | +0.36 | +0.51 | +0.58 | +0.60 |
| random-search | -2.00 | -1.64 | -1.44 | -1.33 | -1.32 |
| greedy asc. | -0.75 | -0.28 | +0.09 | +0.25 | +0.50 |
| beam search | -0.85 | -0.38 | +0.30 | +1.07 | +1.63 |
| genetic algorithm | +0.13 | +0.61 | +1.26 | +1.37 | +1.44 |

## Data files
- `results/report_data/exp1_frozen_summary.csv` (arm×oracle), `exp1_frozen_margins.csv` (per target/rep)
- `results/report_data/exp1_seed_stats.csv`, `exp1_seed_per_seed.csv`
- `results/report_data/exp1_tempsweep.csv`
- `results/report_data/exp1_anytime.csv` + `results/figures/exp1_anytime.png`

Regenerate: `bash scratch/budget_baselines/run_all_exp1.sh`
