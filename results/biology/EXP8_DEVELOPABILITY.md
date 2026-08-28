# Experiment 8 — Developability and sequence-shift diagnostics

_Generated 2026-07-23. Panel = the SVDD validation targets. All CDR-H3 physicochemical/liability
metrics are computed on each design's `cdr_seq` with the identical code as
`report_data/extract_developability.py` (native CDR-H3 as the in-distribution reference). The
three "extra" axes (AbLang2, AGGRESCAN, Hopp-Woods) are whole-VH proxies from
`DEVELOPABILITY_EXTRA.md` / `developability_extra_perdesign.csv` (native heavy chain + design
CDR-H3). **CamSol is ABSENT** (web-server only, not runnable on this host); AGGRESCAN a3v
(aggregation) and Hopp-Woods hydrophilicity (solubility) are used as the sequence-only proxies in
its place. **AbNatiV: installed (`abnativ==2.0.8`) but not runnable** — its VH scorer requires
ANARCI/HMMER numbering, which is not present, so a per-arm AbNatiV humanness score could not be
produced; the AbLang2 VH pseudo-log-likelihood below is the humanness proxy._

**Arms.** `native` (WT CDR-H3), `antifold-baseline` (AntiFold SOTA one-shot), `best-of-N`
(AntiFold prior, best-of-N budget baseline), `GA` (genetic algorithm budget baseline), `SVDD`
(headline `svdd_af_baddg`: SVDD + AntiFold base + BA-DDG escape-robustness steering). GA and
best-of-N sequences are read from the budget-baseline design JSONs
(`scratch/budget_baselines/validation/*_designs_bb_{ga,bestofn}_r*.json`, 51 targets, 408 designs
each); native/antifold/SVDD from `report_data/developability.csv`.

---

## Part 1 — Developability distributions across arms

### Per-arm central tendency (mean over designs; median in Part 1 table)

| metric (↑worse unless noted) | native | antifold | best-of-N | GA | SVDD |
|---|---|---|---|---|---|
| gravy (hydrophobicity) | −0.795 | −0.605 | −0.732 | **−0.985** | −0.582 |
| net_charge | −0.334 | −0.582 | −0.425 | −0.277 | −0.549 |
| aromatic_frac | 0.248 | 0.267 | 0.267 | **0.290** | 0.252 |
| n_cys (count) | 0.006 | 0.000 | 0.015 | 0.012 | 0.012 |
| n_met | 0.318 | 0.357 | 0.297 | 0.306 | 0.343 |
| n_trp | 0.347 | 0.214 | 0.306 | 0.368 | 0.336 |
| nglyc (N-X-S/T) | 0.000 | 0.002 | 0.012 | 0.005 | 0.015 |
| deamid (NG/NS) | 0.068 | 0.029 | 0.034 | 0.034 | 0.059 |
| isomer (DG/DS/DT/DH) | 0.227 | 0.139 | 0.221 | 0.282 | 0.243 |
| cleavage (DP) | 0.011 | 0.011 | 0.005 | 0.007 | 0.007 |
| n_liabilities (sum) | 0.307 | 0.181 | 0.272 | 0.328 | 0.324 |
| designs (n) | 176 | 552 | 408 | 408 | 408 |
| frac ≥1 liability | 0.267 | 0.176 | 0.250 | 0.301 | 0.279 |

### Paired target-level tests vs native (Wilcoxon signed-rank on per-target means)

`mean_delta` = mean over targets of (arm − native); positive = arm has more of the metric.
`flag` = WORSE only when the shift is significant (p<0.05) **and** in the liability-worsening
direction (higher gravy/aromatic/liability-count = worse; lower = better). Metrics that are
all-zero within pairs give no test (—).

| arm | metric | n_pairs | median arm | median native | mean_delta | Wilcoxon p | flag |
|---|---|---|---|---|---|---|---|
| antifold | gravy | 92 | −0.627 | −0.819 | +0.174 | 0.0070 | **WORSE** |
| antifold | net_charge | 92 | 0.00 | 0.00 | −0.301 | 0.0066 | (more negative) |
| antifold | aromatic_frac | 92 | 0.267 | 0.250 | +0.016 | 0.059 | — |
| antifold | n_trp | 92 | 0 | 0 | −0.112 | 0.018 | better |
| antifold | isomer | 92 | 0 | 0 | −0.089 | 0.058 | — |
| antifold | n_liabilities | 92 | 0 | 0 | −0.123 | 0.025 | better |
| best-of-N | aromatic_frac | 51 | 0.267 | 0.250 | +0.038 | 0.0019 | **WORSE** |
| best-of-N | n_cys | 51 | 0 | 0 | +0.015 | 0.031 | **WORSE** (counts) |
| best-of-N | nglyc | 51 | 0 | 0 | +0.012 | 0.063 | — |
| best-of-N | gravy | 51 | −0.782 | −0.819 | +0.006 | 0.93 | — |
| GA | aromatic_frac | 51 | 0.286 | 0.250 | +0.060 | 0.00012 | **WORSE** |
| GA | gravy | 51 | −1.020 | −0.819 | −0.248 | 0.013 | better (more hydrophilic) |
| GA | n_trp | 51 | 0 | 0 | +0.015 | 0.76 | — |
| GA | isomer | 51 | 0 | 0 | −0.012 | 0.87 | — |
| SVDD | aromatic_frac | 51 | 0.250 | 0.250 | +0.023 | 0.022 | **WORSE** |
| SVDD | nglyc | 51 | 0 | 0 | +0.015 | 0.031 | **WORSE** (counts) |
| SVDD | net_charge | 51 | 0.00 | 0.00 | −0.257 | 0.073 | (more negative) |
| SVDD | gravy | 51 | −0.610 | −0.819 | +0.156 | 0.20 | — |

(Metrics not listed for an arm were non-significant, p>0.05.)

**Flags (significant liability worsening vs native):**
- **AntiFold**: only `gravy` worsens (+0.17, more hydrophobic loops); it simultaneously *reduces*
  Trp and total liabilities (both better). Net the most developable arm (frac ≥1 liability 0.18).
- **best-of-N**: `aromatic_frac` ↑ (+0.038) and `n_cys` ↑ — 6 Cys residues introduced across 408
  designs (native ≈ 1). Reported as counts, not "unpaired" (no structural pairing check).
- **GA**: `aromatic_frac` ↑ the most (+0.060) and the highest liability load of any arm
  (n_liabilities 0.328, frac ≥1 liability 0.30); it is the most hydrophilic (gravy −0.985) but
  packs aromatics/isomerization sites.
- **SVDD (ours)**: `aromatic_frac` ↑ (+0.023, smallest of the design arms) and `nglyc` ↑ — 6
  N-glycosylation sequons introduced across 408 designs (native 0). Charge trends more negative
  (p=0.073, n.s.). No significant worsening of Trp, Met, deamidation, isomerization, cleavage, or
  total liability count vs native.

**Cys note (all arms):** total Cys introduced is tiny — native 1, antifold 0, best-of-N 6, GA 5,
SVDD 5 across 408–552 designs. Counts only; no structural confirmation of pairing state, so none
are labelled "unpaired."

### Extra whole-VH proxy axes (AbLang2 / AGGRESCAN / Hopp-Woods)

Paired Wilcoxon vs native on per-target means (n=51 each). `SVDD` here = `svdd_af_baddg_merged`,
`unsteered` = `mx_unsteer`. GA / best-of-N were not scored on these axes (would require an
AbLang2 rerun; no SLURM per scope).

| arm | metric | median arm | median native | mean_delta | Wilcoxon p | flag |
|---|---|---|---|---|---|---|
| antifold | AbLang2 pll (↑better) | −2.425 | −2.420 | −0.004 | 0.15 | — |
| antifold | AGGRESCAN a3v (↓better) | 0.031 | 0.017 | +0.012 | 1.3e-4 | **WORSE** |
| antifold | Hopp-Woods sol (↑better) | −0.213 | −0.205 | −0.017 | 0.0011 | **WORSE** |
| SVDD | AbLang2 pll | −2.433 | −2.420 | −0.022 | 1.7e-5 | **WORSE** |
| SVDD | AGGRESCAN a3v | 0.026 | 0.017 | +0.008 | 0.0054 | **WORSE** |
| SVDD | Hopp-Woods sol | −0.209 | −0.205 | −0.016 | 0.0084 | **WORSE** |
| unsteered | AbLang2 pll | −2.473 | −2.420 | −0.063 | 4.9e-9 | **WORSE** |
| unsteered | AGGRESCAN a3v | 0.019 | 0.017 | +0.004 | 0.031 | **WORSE** |
| unsteered | Hopp-Woods sol | −0.203 | −0.205 | −0.007 | 0.036 | **WORSE** |

**Reading:** every design arm is marginally less human-like (AbLang2), more aggregation-prone
(AGGRESCAN), and less soluble (Hopp-Woods) than the WT loop — expected, since all arms swap a
hydrophilic native loop for a higher-affinity one — but the magnitudes are small. On AbLang2
humanness the ordering is unsteered (−0.063) < SVDD (−0.022) < antifold (−0.004, n.s.): escape
steering costs a little naturalness, but *less* than running the generator unsteered, and SVDD is
statistically indistinguishable from AntiFold on aggregation and solubility (see
`DEVELOPABILITY_EXTRA.md`).

---

## Part 2 — Reward-exploitation test (is BA-DDG being gamed?)

**Design:** per-design BA-DDG reward = worst-case (max) `baddg_ddg` over held-out eval-split
antigen variants (more-negative = stronger/more-robust binder = better). Pooled over the BA-DDG-
steered arms `{mx_baguo, f5_worst_wt1, f5_worst_wt0, f5_cvar_wt1, f5_mean_wt1, one_worst_wt0}`
(n = 10,558 designs, from `grades_baddg_long.csv`). Correlated (Spearman) against per-design
`gravy`, `n_trp`, `nmut_vs_native` (edit distance from developability.csv), and `mpnn_mean` (the
base-model score in `designs_meta.csv`). Because reward level is dominated by which target a
design belongs to, the decisive column is **target-controlled** (both reward and predictor
mean-centered within each target, then pooled) — this isolates "within a target, do higher-Trp /
more-hydrophobic / more-edited designs earn better reward?", the actual gaming question.

_On "AntiFold log-prob": no AntiFold sequence log-likelihood column exists in the released
`designs_meta`/`overlay` tables for these steered arms. `mpnn_mean` is the closest available
prior/base score; it correlates −0.53 with `wt_ddg` (WT binding), so it behaves as a base
**binding** oracle rather than a pure naturalness likelihood — reported below with that caveat._

| predictor | Spearman (raw, pooled) | p (raw) | Spearman (target-controlled) | p (ctrl) |
|---|---|---|---|---|
| hydrophobicity (gravy) | +0.139 | 7e-47 | **+0.011** | 0.24 |
| n_Trp | +0.027 | 5e-3 | **+0.009** | 0.34 |
| edit distance (nmut) | +0.179 | 1e-76 | **+0.013** | 0.17 |
| base-model prior (mpnn_mean) | −0.614 | ~0 | −0.376 | ~0 |

**Interpretation — the reward is NOT gamed via obvious sequence artifacts:**
- **Hydrophobicity, Trp, edits:** all three target-controlled correlations are ≈0 and
  non-significant (ρ = 0.009–0.013, p ≥ 0.17). Within a target, adding Trp, raising GRAVY, or
  making more edits buys *no* reward. The raw pooled positive values (gravy +0.14, nmut +0.18)
  are cross-target confounds and, tellingly, point the *wrong way for gaming*: sign is positive,
  i.e. more hydrophobic / more-edited designs get slightly **worse** worst-case BA-DDG. There is
  no Trp-stacking, hydrophobicity, or large-edit exploit signature.
- **Base-model prior (mpnn_mean):** strongly *aligned* with reward (higher prior score → more
  negative/better BA-DDG, ρ = −0.61 raw, −0.38 target-controlled). Whether read as prior
  plausibility or simply as a second binding oracle, this is the opposite of off-manifold gaming
  (gaming would show reward improving as the prior score *falls*). Better BA-DDG comes from
  more-prior-favored sequences, not adversarial ones.

This is consistent with Part 1: the steered SVDD loops carry roughly native-level liabilities and
only a small AbLang2/aggregation penalty — not the profile of an oracle-gamed adversarial sequence.

---

## Methods (brief)
- CDR-H3 metrics: Kyte-Doolittle GRAVY, net charge (K+R−D−E+0.1·H), aromatic fraction (F+W+Y),
  and regex liability motifs (N-glyc N[^P][ST]; deamidation N[GS]; isomerization D[GSTH]; cleavage
  DP), identical to `extract_developability.py`. GA/best-of-N computed the same way from budget
  JSON `cdr_seq`.
- Paired tests: per-target mean per metric, Wilcoxon signed-rank (scipy) vs native on the target
  intersection; zero-difference-only pairs → no test.
- Extra axes: AbLang2 `ablang2-paired` VH pseudo-log-likelihood, AGGRESCAN a3v mean, Hopp-Woods
  mean, from the precomputed per-design CSV; paired Wilcoxon vs native.
- Reward test: worst-case eval-split `baddg_ddg` per design; Spearman raw and within-target-
  centered.
- Tooling notes: **CamSol** unavailable (web server) → AGGRESCAN + Hopp-Woods proxies. **AbNatiV**
  installed (2.0.8) but unrunnable (needs ANARCI numbering; absent) → AbLang2 pll is the humanness
  proxy.
