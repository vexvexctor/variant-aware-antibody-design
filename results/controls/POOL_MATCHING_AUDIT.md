# Pool-matching audit — 2026-08-22

Bookkeeping check over existing results. No GPU. Question: are the arms in each reported
comparison actually matched on rep count and pool size, or only assumed to be?

## Step 1-2. Rep counts and mismatches

Counted graded reps per (target, arm) for every arm in the Table 3 and Table 5 comparisons.

| comparison | targets, all arms graded | missing an arm | **rep-mismatched** |
|---|---|---|---|
| T3 AntiFold | 190 | 10 | **0** |
| T3 AbMPNN | 190 | 10 | **0** |
| T3 ProteinMPNN | 190 | 10 | **0** |
| T3 ESM-IF1 | 190 | 10 | **0** |
| T3 diffusion | 190 | 10 | **0** |
| T5 compute-matched | 190 | 10 | **0** |

**Zero mismatches.** The reason is structural, not luck: `pool_worsts` requires FULL held-out
variant coverage and drops any (target, rep) that lacks it, and the arms with gaps are missing
all three reps of the affected target rather than one. Such a target is excluded wholesale
from every arm. Every one of the 190 comparison targets carries exactly 3 reps x 6 designs
= 18 designs in each 6-design arm.

## Step 3. The missing bb_bestofn cells

Expected 600 (target, rep); graded 570; **30 missing** — not 24. Breakdown:

- design JSON absent: **0**
- design present but ungraded: **24**  <- the 24
- graded but dropped for partial held-out-variant coverage: **6**

All 30 fall on just **10 targets, all 3 reps each**. Those are exactly the 10 "missing an arm"
targets above, so they are excluded from the comparison rather than contributing a short pool.
The 24 ungraded cells are the tail of the grading chain, which stops after 6 rounds.

## Step 4-5. Recompute under rep-truncation, diff vs published

Recomputed every rate restricting all arms to the INTERSECTION of their rep ids, and diffed.

| comparison | as published | rep-matched | diff |
|---|---|---|---|
| AntiFold: ours / base / delta | 64.3 / 50.1 / +14.2 | 64.3 / 50.1 / +14.2 | **0.00** |
| AbMPNN | 83.4 / 77.1 / +6.3 | 83.4 / 77.1 / +6.3 | **0.00** |
| ProteinMPNN | 62.6 / 53.5 / +9.1 | 62.6 / 53.5 / +9.1 | **0.00** |
| ESM-IF1 | 61.6 / 41.5 / +20.1 | 61.6 / 41.5 / +20.1 | **0.00** |
| diffusion | 12.6 / 1.9 / +10.7 | 12.6 / 1.9 / +10.7 | **0.00** |
| compute-matched: ours - bestofn | -16.9 | -16.9 | **0.00** |
| compute-matched: ours - p36 | -0.4 | -0.4 | **0.00** |

**Every number is identical to two decimals. No table changes.**

Pool matching is now VERIFIED rather than inferred, on both the aggregator path (which already
intersected reps in `emit_pairs`) and the ad-hoc compute-matched analysis (which did not
intersect, and turns out not to have needed to, because the target sets are uniform).

## The one real caveat: p36 is NOT pool-matched

Rep counts match, but design counts do not:

| arm | designs per target (H3, 3 reps) | designs per target (Pythia, r1) |
|---|---|---|
| ours_antifold | 18 | 6 |
| bb_bestofn | 18 | 6 |
| base_antifold | 18 | 6 |
| **base_antifold_p36** | **108** | **36** |

p36 carries a deliberate **6x pool advantage** — it is the "does more sampling alone close the
gap" control, and the advantage is what makes it conservative. But its deltas
(`ours - p36 = -0.4 pp` on H3, `-1.7 pp` on Pythia) must NOT be described as pool-matched.
Correct phrasing: *"ours ties a baseline given six times the sampling budget and six times the
selection pool."* The bb_bestofn comparison IS pool-matched (18 vs 18) and reward-budget matched,
so it is the one that carries the compute-matched claim.

## Verdict

Nothing to update in Tables 3 or 5. State pool-matching as verified. Fix the p36 wording in the
archival version.
