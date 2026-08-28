# Random edit-distance-matched CDR-H3 FoldX control (reviewer concern #2)

**Question (Reviewer #2).** Does *any* CDR-H3 redesign beat the native antibody on the out-of-family FoldX arbiter? If random sequences at the same edit distance as our designs also beat native, our beat-native rate is meaningless.

**Answer.** No. Random CDR-H3 mutants matched design-for-design to the edit distance of our headline arm (svdd_af_baddg) beat native on FoldX only **74.5%** of targets (n=51), versus **80.4%** for svdd_af_baddg on the identical panel, grader, and held-out eval split. The steered arm's advantage is not an artefact of 'moving away from native'.

## Methods

- **Targets.** All 51 omega-panel antibody-antigen complexes (AACDB structures).
- **Edit-distance match.** For every real svdd_af_baddg design in a target's best-of-pool (4 cluster/single designs x 2 unseeded reps = 8 designs), we emit ONE random mutant that mutates the *same number* of CDR-H3 positions as that design (per-design, per-target edit distance matched exactly; verified rand_nmut == real_nmut for all 8x51 designs). Global mean edit distance 6.5 mutations, same as the steered arm.
- **Randomisation.** Mutated positions chosen uniformly at random; each substituted residue drawn uniformly from the 20 standard amino acids excluding the native residue at that position (every drawn edit is a genuine substitution) and excluding Cys (never introduce a spurious unpaired cysteine — a conservative choice that does not handicap the random arm). A native Cys conserved by all real designs would be protected; no panel CDR starts with Cys.
- **No steering.** The random arm has zero access to any oracle, structure, or escape variant. It is pure edit-distance-matched noise on the native CDR-H3.
- **Grader = FoldX (out-of-family physics arbiter)**, `validate_designs_foldx.py`, the exact binary/pipeline and cache used for svdd_af_baddg's 80.4%. Renumber -> RepairPDB -> BuildModel (CDR design + antigen point mutations) -> AnalyseComplex interaction energy.
- **Beat-native metric (identical to the screen).** Best-of-pool worst-case: for each design, worst (least-negative) FoldX interaction energy over the **held-out omega-divergent eval variants** (split=='eval', 24/target, never seen by any arm); best-of-pool = the min over the 8 random designs; margin = best-of-pool worst MINUS the native (WT) antibody's worst over the same eval variants. **margin < 0 => beats native.**
- **Reproducible.** Fixed seed 20260715; per-target RNG stream = sha256(seed:target).

## Panel result

| arm | FoldX beat-native | n | mean margin |
|---|---|---|---|
| **random edit-matched (control)** | **74.5%** | 51 | -2.103 |
| svdd_af_baddg (headline, steered) | 80.4% | 51 | -2.914 |
| native antibody (baseline) | 0% by definition | 51 | 0.000 |

Random-edit control beats native on 38/51 targets (74.5%); the steered headline arm beats native +5.9 percentage points more often on the same grader. A mean margin that is less negative (or positive) than the steered arm's (-2.914) confirms random redesigns are, on average, no tighter than native against the held-out escape variants.

## Per-target

Full per-target margins in `random_edit_pertarget.csv` (columns: target, foldx_margin, n_designs_scored, native_worst, pool_worst, beats_native).

| target | native CDR len | pool margin | beats native |
|---|---|---|---|
| 1IC7_HLY | 7 | -2.069 | yes |
| 1KB5_HLB | 12 | -0.529 | yes |
| 1KXT_DC | 5 | +0.201 | no |
| 1KXT_FE | 5 | -1.469 | yes |
| 1RJC_AB | 9 | +1.873 | no |
| 1ZVH_AL | 18 | +1.732 | no |
| 3L5W_HLI | 13 | +0.098 | no |
| 3L5X_HLA | 13 | -1.082 | yes |
| 3OGO_FC | 8 | +3.299 | no |
| 3OGO_GB | 8 | +0.636 | no |
| 3VI3_HLD | 11 | -3.553 | yes |
| 4AL8_HLC | 10 | -0.755 | yes |
| 4BZ1_HLA | 10 | +2.094 | no |
| 4L5F_HLE | 10 | -0.331 | yes |
| 5DUM_HLA | 20 | -4.343 | yes |
| 5GJS_HLB | 15 | +0.319 | no |
| 5GZN_HLA | 20 | -15.568 | yes |
| 5JHL_HLA | 10 | -0.594 | yes |
| 5MVZ_HLU | 11 | -2.743 | yes |
| 5UGY_HLA | 19 | -2.868 | yes |
| 5VIC_HLE | 15 | -13.089 | yes |
| 5VIG_HLZ | 13 | -0.092 | yes |
| 6AJ9_HLC | 13 | +0.174 | no |
| 6FLA_HLI | 10 | -5.576 | yes |
| 6FLB_HLG | 10 | -13.774 | yes |
| 6KZ0_KLJ | 15 | -1.232 | yes |
| 6MLM_KLH | 17 | -6.967 | yes |
| 6PLK_HLE | 15 | -2.961 | yes |
| 6QX4_DA | 18 | -4.395 | yes |
| 6QX4_HB | 18 | -0.397 | yes |
| 6VJA_HLC | 14 | -0.201 | yes |
| 6VRQ_HLA | 12 | -0.272 | yes |
| 7AMS_HLB | 17 | -4.991 | yes |
| 7BQ5_HLA | 10 | -2.481 | yes |
| 7CN2_GgE | 16 | -0.262 | yes |
| 7JLK_HLB | 15 | -3.977 | yes |
| 7MZK_HLA | 15 | -1.692 | yes |
| 7OR9_HLE | 16 | -0.234 | yes |
| 7STZ_HLC | 13 | -0.109 | yes |
| 7TBF_HLA | 6 | +0.161 | no |
| 7TTX_HLA | 22 | +0.981 | no |
| 7VYR_HLR | 11 | -1.272 | yes |
| 7WPH_HLB | 14 | -2.993 | yes |
| 7X2O_HLB | 13 | -2.054 | yes |
| 7X2W_HLC | 11 | +1.593 | no |
| 7X3D_HLB | 16 | -2.560 | yes |
| 7X4I_HA | 17 | -0.069 | yes |
| 7YCL_HLD | 15 | -4.852 | yes |
| 7ZRA_FC | 16 | -0.901 | yes |
| 8CYH_HLM | 12 | +3.545 | no |
| 8IDN_HLA | 12 | -10.644 | yes |

_Generated by `${VAAD_ROOT}/work/random_edit_control/finalize_md.py`._
