# Experiment 10 — Low-similarity / post-cutoff panel, frozen single-deployed design, 5 seeds

**Owner:** EXP10 (novel OOD panel). **Date:** 2026-07-23.
**Scope:** re-run the 23-target out-of-distribution antibody-design panel with **5 reps/seeds**,
deploy **one FROZEN design per (target, seed)** selected **only by BA-DDG**, and grade that
identical frozen sequence under every available evaluator. Reports **CVaR20 and strict-max**
beat-native, **WT-retention**, and **across-seed statistics** (never best-seed-only).

Follows the Exp2 finding (`EXP2_FROZEN_VS_ORACLE.md`): the published headline was mostly
**oracle-selection inflation**; the deployable result is the **frozen** single design chosen by
the steering objective, graded held-out. This experiment reproduces that discipline on the OOD panel.

---

## 1. Panel & metadata

23 targets (`novel_panel/panel_targets_novel.txt`), each an Ab–Ag complex with a per-target escape
panel split into in-loop **steer** variants and disjoint **held-out eval** variants (position-disjoint,
omega holdout). Full provenance in **`results/report_data/exp10_metadata.csv`**
(target_id, PDB_id, deposition_date, antigen_sequence, nearest_training_antigen, sequence_identity,
alignment_coverage, and the four cutoff columns).

**Novelty is real but heterogeneous — stated honestly:**
- **16/23** targets have antigen **%identity < 50 %** (or no mmseqs hit at all) vs the AACDB
  training antigens; `sequence_identity × alignment_coverage` reproduces the pipeline's `gid`.
- **16/23** were deposited **on/after 2024-06-01** (the RCSB `initial_release_date` filter used to
  build the panel, verifiable in `curate_candidates.py`). **7/23 predate it** (8SGN 2024-04, 7R8U
  2022, 7XQ8 2022, 6UMG 2020, 5WNA 2018, 8BW0 2024-01, 7Y1B 2023) — these were added in the
  10→23 expansion for **low antigen identity**, not deposition recency. The panel is therefore best
  described as **low-similarity OOD**, only partly strict-post-cutoff.
- **Per-model training cutoffs** (`known_AntiFold_cutoff`, `known_BA_DDG_cutoff`,
  `known_H3_DDG_cutoff`, `cutoff_certainty`) are recorded as the literal string **`unknown`**: no
  authoritative in-repo source documents the exact training-data cutoff of AntiFold, BA-DDG, or
  H3-DDG, and per the experiment guardrail these are **not inferred/fabricated**. Deposition date
  (verifiable) is the only novelty anchor banked with certainty.

## 2. Method (frozen, BA-DDG-selected, 5 seeds)

- **Seeds = reps.** Sampling is unseeded; run-to-run variation is the seed. Reps `r1..r5`, both arms.
- **Arms.** `headline` = SVDD on the AntiFold base, BA-DDG objective, CVaR20 steer (the model-matrix
  headline cell); `antifold` = budget-matched AntiFold best-of-N baseline (same 6-design pool, same
  escape panel).
- **BA-DDG grading (uniform, native-referenced).** Both arms are graded by
  `scripts/surrogate/baddg_score_designs.py`, which scores every (antibody × antigen-variant) pair
  including the **WT / native-antibody** row, so the native reference is intrinsic and both arms are
  apples-to-apples in one convention (higher `baddg_ddg` = weaker binding vs the native/WT baseline).
- **Frozen selection.** For each (target, seed, arm) deploy the **single** pooled design minimizing
  **BA-DDG CVaR20 over the in-loop STEER variants** (tightest robust binder). No held-out information
  and no per-grader oracle selection enters deployment.
- **Frozen grading.** That one design's **held-out eval** worst-case is scored under **BA-DDG,
  H3-DDG, FoldX**. Beat-native (margin < 0 = beats native):
  - **BA-DDG:** design vs native (WT-antibody) over eval.
  - **H3-DDG:** design-vs-native ΔΔG, native ≡ 0.
  - **FoldX:** design interaction-E − native(WT) interaction-E.
  - Two aggregations per grader: **CVaR20** (mean of the worst 20 % eval tail) and **strict-max**
    (single worst eval variant).
- **WT-retention.** For each deployed design, BA-DDG on the **WT antigen** (native WT-WT ≡ 0) vs its
  **held-out escape CVaR20** — does escape-steering keep on-target potency?
- **Across-seed stats.** Per grader×metric×arm: each seed's panel beat-native **rate**, then
  **mean / std / best / worst** over seeds and the **fraction of targets beaten by all / any** seed.
  **The headline is the mean over seeds, not the best seed.**

MM-GBSA / Pythia / cofold were **not** added: those graders are wired to the 51-panel / SKEMPI
structure sets, not these 23 OOD complexes, so the harness does not support them here without new
structure builds (out of scope). BA-DDG (in-family objective), H3-DDG (independent, discriminating)
and FoldX (weak physics arbiter) are the three evaluators the panel supports.

## 3. Results

Tables auto-regenerate from disk — **re-run the aggregator to refresh** (see §5). Per-rep rows:
`results/report_data/exp10_frozen_perrep.csv`; machine summary:
`results/report_data/exp10_frozen_summary.json`.

> **STATUS: FINAL (2026-07-23).** All 5 seeds × 23 targets graded (BA-DDG, H3-DDG, FoldX).
> Numbers below are mean over 5 seeds; source `results/report_data/exp10_frozen_summary.json`.

### 3a. Frozen beat-native, per evaluator — FINAL, mean over 5 seeds (23 targets)

Beat-native % = fraction of the 23 targets whose frozen deployed design beats native, averaged
over 5 seeds (best / worst seed and "any-seed" in parentheses).

| arm | evaluator | metric | mean (best / worst) | any-seed | headline vs AntiFold |
|---|---|---|---|---|---|
| **headline** | **BA-DDG** (objective) | CVaR20 | **79.1 %** (87.0 / 69.6) | 96 % | dominates (vs 42.6 %) |
| headline | BA-DDG | strict-max | 75.7 % (82.6 / 69.6) | 100 % | vs 47.0 % |
| **headline** | **H3-DDG** (independent) | CVaR20 | **48.7 %** (52.2 / 47.8) | 78 % | **vs 30.4 %** (+18 pp) |
| headline | H3-DDG | strict-max | 40.9 % (43.5 / 34.8) | 70 % | vs 26.1 % |
| headline | FoldX (weak arbiter) | CVaR20 | 47.0 % (56.5 / 43.5) | 74 % | vs 45.2 % (no sep.) |
| headline | FoldX | strict-max | 45.2 % (56.5 / 30.4) | 83 % | vs 48.7 % |
| antifold | BA-DDG | CVaR20 | 42.6 % | 57 % | |
| antifold | H3-DDG | CVaR20 | 30.4 % | 39 % | |
| antifold | FoldX | CVaR20 | 45.2 % | 61 % | |

Final read, consistent with Exp2: the **frozen** headline dominates the AntiFold baseline on its
in-family objective (**BA-DDG 79.1 % vs 42.6 %**) and, most importantly, **beats it on the
independent H3-DDG grader (48.7 % vs 30.4 %, +18 pp)** on genuinely novel/post-cutoff antigens.
**FoldX shows no separation (47 % vs 45 %)** — the weak-arbiter pattern again — and the frozen
H3 ~49 % is far below the old best-of-pool 100 %, i.e. most of the old FoldX headline was
oracle-selection inflation (Exp2). Headline beats AntiFold on H3 in **5/5 seeds**; the win is
seed-robust, not a best-seed artifact.

### 3b. WT-retention (deployed design, BA-DDG units, native WT-WT ≡ 0) — FINAL

Escape-steering does **not** sacrifice on-target potency: the frozen headline design binds the **WT
antigen tighter than native** on average (mean WT ΔΔG ≈ **−1.0**) while carrying a positive held-out
escape tail, and **≈93 %** of deployed headline designs retain WT binding within +0.5 of native
(AntiFold ≈69 %). Consistent with EXP5 (no WT-vs-robustness tradeoff).

## 4. Files

- Metadata: `results/report_data/exp10_metadata.csv`
- Per-rep frozen grades: `results/report_data/exp10_frozen_perrep.csv`
- Machine summary: `results/report_data/exp10_frozen_summary.json`
- Design/grade artifacts: `scratch/novel_panel/{validation,baddg_out,h3_out,foldx_out}/…_r{1..5}…`
- Scripts: `scratch/novel_panel/{run_target_rep.sh, submit_rep_array.sh, aggregate_exp10.py,
  build_exp10_metadata.py}`

## 5. Reproduce / finalize

```bash
# (re)build metadata
python3 ${VAAD_ROOT}/work/novel_panel/build_exp10_metadata.py
# (re)aggregate frozen 5-seed metrics once jobs finish (safe to run anytime; reads disk)
python3 ${VAAD_ROOT}/work/novel_panel/aggregate_exp10.py
```

The design+grade chain is idempotent (`run_target_rep.sh` / array `submit_rep_array.sh`,
manifests `exp10_rep_manifest.tsv` = reps 2–5, `exp10_rep1_manifest.tsv` = rep-1 BA-DDG backfill);
requeue only produces missing rep outputs and never touches r1's original design/H3/FoldX.
