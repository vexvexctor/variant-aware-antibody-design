# Design-panel-size sensitivity — does M_design = 6 pay for itself?

_Run 2026-08-25. 51-target panel × M_design ∈ {1, 3, 6, 12} × 3 unseeded reps = 612 design runs,
612 H3-DDG gradings. Zero failures._

**Short answer: no, and the curve does not saturate at six.** Beat-native is flat from M=1 to
M=12 (no pairwise difference is significant). The continuous margin *does* improve with M, but
only M=12 separates from M=6 — and it costs 2× the oracle calls to get there. **M=6 is the one
setting that is hard to defend on this evidence**: it is indistinguishable from M=1 at 5× the
cost, and significantly worse than M=12 on margin.

## Protocol

Everything except `--max-cluster` is pinned to the headline cell `svdd_af_baddg`:
base = AntiFold (frozen npz prior), oracle = BA-DDG, `--steer-agg cvar --cvar-alpha 0.2`,
`--n-designs 6`, `--svdd-k 6`, `--holdout-by omega`, `--eval-cluster 24`.

Grader = **H3-DDG**, out-of-family relative to the BA-DDG steering objective (project grader
rule: steer BA → grade H3). Margin = design-minus-native; **< 0 = design binds that held-out
variant tighter than the native affinity-matured antibody**. beat-native = best-of-pool
worst-case margin < 0. All reps pooled at a **matched 3 reps per arm** (best-of-pool is monotone
in pool size, so this matters).

CI = 95% bootstrap over **antigen clusters** (30 distinct), not targets.

### Validity: the held-out panel really is identical at every M

This is the whole basis of the comparison, so it was checked rather than assumed. The steer
panel is the **top-M by omega**; the held-out panel is the **bottom-24 by omega**. Every target
on the panel has ≥ 41 variants (min 41, median 4,715), so with M ≤ 12 the two sets cannot
intersect. Verified on the generated JSONs:

- **51/51 targets have a byte-identical held-out eval set across M = 1, 3, 6, 12** (md5 of the
  ordered variant-name list).
- **0 steer/eval overlaps** at any M.
- Every eval panel is exactly **24** variants.

## Results

| M_design | beat-native (pooled) | 95% CI | per-run | mean held-out margin | oracle calls / run | CUDA s / run |
|---:|---:|---|---:|---:|---:|---:|
| **1** | **70.6%** | [56.7, 84.4] | 55.6% | −0.056 | **441** | 12.0 |
| 3 | 66.7% | [54.4, 78.8] | 54.2% | −0.108 | 1,198 | 34.2 |
| **6** (current) | 66.7% | [51.2, 78.7] | 51.0% | −0.102 | 2,334 | 67.7 |
| 12 | 64.7% | [48.9, 78.4] | 54.9% | **−0.132** | 4,599 | 134.1 |

n = 51 targets / 153 runs / 30 antigen clusters per arm. `reward_evals` is ~441 at every M; the
cost that scales is **oracle calls** (= reward evals × M), which is the hardware-independent
budget axis.

### Paired tests (same targets, cluster-bootstrapped)

Marginal CIs overlap heavily, which proves nothing either way. Pairing per target is the
powered test:

| pair | Δ margin | 95% CI | Δ beat-native | 95% CI | significant? |
|---|---:|---|---:|---|:--:|
| M=6 → M=1 | −0.025 | [−0.182, +0.100] | +3.9 pp | [−4.6, +15.6] | **no** |
| M=6 → M=3 | −0.028 | [−0.122, +0.079] | +0.0 pp | [−6.9, +8.7] | **no** |
| M=6 → M=12 | **−0.154** | **[−0.282, −0.040]** | −2.0 pp | [−8.7, +4.9] | **margin only** |
| M=1 → M=3 | −0.003 | [−0.137, +0.164] | −3.9 pp | [−14.5, +5.4] | no |
| M=1 → M=6 | +0.025 | [−0.100, +0.182] | −3.9 pp | [−15.6, +4.6] | no |
| M=1 → M=12 | **−0.129** | **[−0.232, −0.010]** | −5.9 pp | [−18.4, +2.7] | **margin only** |

Negative Δ margin = the second arm is tighter/better.

## Reading

1. **Beat-native does not respond to M at all.** No pairwise contrast is significant, and the
   point estimate mildly *decreases* with M. Binarising at zero throws away the signal.
2. **The continuous margin does respond, but only at M=12.** M=12 is significantly tighter than
   both M=6 (−0.154) and M=1 (−0.129). M=3 and M=6 are indistinguishable from M=1 and from each
   other. So the gain is real but small, and it arrives late — there is **no knee at six**.
3. **Almost all of the variant-aware benefit is already present at M=1.** One escape variant in
   the objective (plus WT) buys what six do. For reference, EXP4 measured WT-only steering at
   9.8% beat-native on the same panel — so the jump from "no variants" to "one variant" is
   enormous, and everything after that is nearly flat.
4. **Cost scales linearly and buys nothing on the headline metric.** M=6 spends 5.3× the oracle
   calls of M=1 and 5.6× the CUDA seconds for +0.0 to −3.9 pp (n.s.).

## What this means for the paper

The hoped-for defence — "performance saturates at six, so six is the natural choice" — **is not
what the data shows**, and claiming it would not survive a reviewer running this sweep.

Two honest options:

- **Defend M=1 on cost.** Same beat-native, 5× cheaper, and it strengthens the training-free /
  cheap-inference story. The cost is a slightly weaker mean margin (−0.056 vs −0.102).
- **Defend M=12 on margin.** The only arm with a statistically real improvement, at 2× the cost
  of the current setting.

**M=6 is the weakest of the three to argue for.** If it stays, the honest wording is that the
method is *insensitive* to design-panel size over 1–12, and six was a default rather than a
tuned optimum — which is a perfectly publishable robustness statement, just not a
performance one.

One caveat linking to the known benchmark limitation: held-out escape correlates with plain WT
binding at r = +0.72 on this panel. If most of what the held-out metric measures is affinity,
then adding escape variants to the design objective *should* have limited returns — this result
is consistent with that coupling, and a depth-1 (single-mutant) panel could plausibly restore a
real dependence on M.

## Files

| file | what |
|---|---|
| `PANEL_SIZE_SENSITIVITY.md` | this write-up |
| `panel_size_sensitivity.png` | beat-native vs M, margin vs M, performance vs cost |
| `panel_size_summary.csv` | one row per M: rates, CIs, costs |
| `panel_size_per_unit.csv` | 612 rows — per (target, M, rep) margin, beat-native flag, cost counters |
| `design_array.sh` | the design sweep (only `--max-cluster` varies) |
| `grade_h3_array.sh`, `make_shards.sh` | H3-DDG grading |
| `aggregate_panel_size.py`, `paired_test.py`, `plot_panel_size.py` | analysis |

Reproduce: `sbatch --array=0-611%48 design_array.sh` → `./make_shards.sh 48` →
`sbatch --array=0-47%48 grade_h3_array.sh` → `python3 aggregate_panel_size.py && python3 paired_test.py && python3 plot_panel_size.py`

_Note: the H3 grader needs Bio + torch + yaml; use `${VAAD_ROOT}/tools/boltz-env/bin/python3`
(the default python lacks Bio, `scratch/lyra_env` lacks yaml)._
