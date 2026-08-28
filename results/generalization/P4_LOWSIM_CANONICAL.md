# P4 — Canonical low-similarity (novel / post-cutoff) results

_Re-analysis only. Builds ONE canonical long file `report_data/low_sim_canonical.csv` for the **23-target novel post-cutoff antigen panel** (antigens with ~0–49% identity to train, post-2024-06), then regenerates every previously-divergent low-sim number from it. Convention: `margin < 0` = beats native._

## The canonical file

`low_sim_canonical.csv` — one row per **(target, seed, method, grader, agg_metric)**, 1,955 rows, 23 targets, 5 seeds. Columns: `target, seed, method, grader, agg_metric, frozen_selected_margin, oracle_best_margin, pool_size, selection_objective, split_id, n_eval`.

- **method** ∈ {`bestofn`, `svdd` (=the SVDD "headline" arm), `antifold` (baseline)}

- **grader** ∈ {`baddg` (in-family), `h3` (independent, ProteinMPNN-lineage), `foldx` (weak arbiter)}; **agg_metric** ∈ {`cvar20`, `max`} is the worst-case aggregation over held-out eval variants (kept as its own column so cvar20/max never collapse — the single biggest source of the divergent numbers).

- **frozen_selected_margin** = the ONE deployed (frozen-selected) design's held-out margin for that seed. **oracle_best_margin** = best (min) margin over the seed pool for that (method,target,grader,metric) — the cheating best-of-pool upper bound (target-level, repeated across the 5 seed rows). **selection_objective** = `baddg_cvar20` (the deploy rule for all arms). **split_id** = `novel_postcutoff_23`. **pool_size** = 5.

Sources merged: `exp10_frozen_perrep.csv` (svdd + antifold, all 3 graders × 2 metrics), `tier1_2b_novel_bestofn_perrep.csv` (bestofn: baddg-cvar20, h3-cvar20, h3-max) and `tier1_2b_novel_bestofn_foldx.csv` (bestofn: foldx cvar20+max).

## Two distinct beat-native definitions (the root of the divergence)

Every low-sim number is one of two things, and files that quoted them rarely said which:

1. **FROZEN (deployable):** does the single deployed design beat native? Reported as the **5-seed mean** of the per-seed 23-target rate (≈ per-rep rate over the 115 (target,seed) cells).

2. **ORACLE best-of-pool (upper bound, NOT deployable):** does the best-of-5-seeds design beat native, per target (n=23)? This is systematically higher and is where the old inflated headlines came from.

### Regenerated FROZEN 5-seed-mean beat-native (per-rep, of 115)

| method | grader | cvar20 | strict-max |
|---|---|--:|--:|
| bestofn | baddg | 66.1% (76/115) | — |
| bestofn | h3 | 59.1% (68/115) | 54.8% (63/115) |
| bestofn | foldx | 49.6% (57/115) | 45.2% (52/115) |
| svdd | baddg | 79.1% (91/115) | 75.7% (87/115) |
| svdd | h3 | 48.7% (56/115) | 40.9% (47/115) |
| svdd | foldx | 47.0% (54/115) | 45.2% (52/115) |
| antifold | baddg | 42.6% (49/115) | 47.0% (54/115) |
| antifold | h3 | 30.4% (35/115) | 26.1% (30/115) |
| antifold | foldx | 45.2% (52/115) | 48.7% (56/115) |

### Regenerated ORACLE best-of-pool beat-native (per target, of 23)

| method | grader | cvar20 | strict-max |
|---|---|--:|--:|
| bestofn | baddg | 82.6% (19/23) | — |
| bestofn | h3 | 82.6% (19/23) | 82.6% (19/23) |
| bestofn | foldx | 73.9% (17/23) | 78.3% (18/23) |
| svdd | baddg | 95.7% (22/23) | 100.0% (23/23) |
| svdd | h3 | 78.3% (18/23) | 69.6% (16/23) |
| svdd | foldx | 73.9% (17/23) | 82.6% (19/23) |
| antifold | baddg | 56.5% (13/23) | 60.9% (14/23) |
| antifold | h3 | 39.1% (9/23) | 26.1% (6/23) |
| antifold | foldx | 60.9% (14/23) | 56.5% (13/23) |

## Reconciliation ledger — every divergent number traced

| reported number | file / context | what it actually is | canonical value | reconciles? |
|---|---|---|---|:--:|
| **48.7% vs 30.4%** | EXP10_POSTCUTOFF, REFRAME_EVIDENCE, TIER1 panels | frozen 5-seed-mean, **H3 cvar20**: svdd(headline) vs antifold | 48.7% (56/115) vs 30.4% (35/115) | ✅ exact |
| **52 vs 30** ("another table") | EXP10 `48.7% (52.2/47.8)` parenthetical | best-seed / worst-seed of the SAME svdd H3 cvar20 cell (not a new metric); 30 = antifold | best-seed 52.2% (12/23), worst 47.8% (11/23); antifold 30.4% | ✅ (it was a best-seed, not a separate estimate) |
| **47.0 vs 45.2** | EXP10 / TIER1_2B_FOLDX_CONFIRM ("no separation") | frozen 5-seed-mean, **FoldX**: cvar20 vs strict-max of svdd (and ≈ antifold/bestofn — FoldX doesn't discriminate) | svdd FoldX cvar20 47.0% (54/115), max 45.2% (52/115); antifold foldx max 48.7%, cvar20 45.2% | ✅ exact |
| **oracle 78.3%** | (implied best-of-pool) | ORACLE best-of-pool, **H3 cvar20, svdd** (18/23) | 78.3% (18/23) | ✅ exact |
| **oracle 91.3%** (21/23) | REVIEWER2_RESPONSE.md, citing `NOVEL_OOD_RESULTS_23.md` | **that file DOES NOT EXIST** (only the 10-target `NOVEL_OOD_RESULTS.md` does) → phantom cell. Intended cell = svdd FoldX oracle-best-of-pool on the 23-panel = **19/23 (82.6%)**; AntiFold = **15/23 (65.2%)** | ✅ **RESOLVED — regenerated from canonical; 21/23 & 91.3% cut** |
| **59.1%** | TIER1_2B panels (best-of-N H3) | frozen 5-seed-mean, **H3 cvar20, bestofn** | 59.1% (68/115) | ✅ exact |

### The one number that does not reconcile

**`91.3%` (21/23) RESOLVED — it was a phantom.** It appears only in `REVIEWER2_RESPONSE.md`, citing
`NOVEL_OOD_RESULTS_23.md` — a file that **does not exist** (the only novel-OOD file is the 10-target
`NOVEL_OOD_RESULTS.md`). Regenerated from the canonical 23-panel: the intended headline cell (svdd
FoldX oracle-best-of-pool) is **19/23 (82.6%)**, and the AntiFold comparator is **15/23 (65.2%)**, not
18/23. Both `91.3%`/`21/23` and the `18/23` AntiFold figure are **cut**; source fixed in
`REVIEWER2_RESPONSE.md`. (For reference, H3 on this panel: svdd 18/23=78.3%, AntiFold 9/23=39.1%.)

### Why the numbers looked contradictory

They were not contradictory — they were **different (arm × grader × cvar20/max × frozen-vs-oracle) cells quoted without their qualifiers.** The canonical file makes all of them a single lookup: `52` was the best seed of the same cell whose mean is `48.7`; `47.0/45.2` are FoldX cvar-vs-max of one arm; `78.3` is the oracle upper bound of the `48.7` frozen cell. Only `91.3` has no cell.
