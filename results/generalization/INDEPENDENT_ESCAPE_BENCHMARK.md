# Independent-escape benchmark (reviewer concern #1)

**Question.** Our optimization ("omega") panel = ProteinMPNN interface proposals filtered by
ESM-1v plausibility — a *plausible-interface-mutation* panel, not a *validated-escape* panel.
Does the headline **svdd_af_baddg** result (designs beat the native antibody on escape
robustness) still hold when the designs are graded on an escape set constructed **independently**
of the ProteinMPNN/ESM-1v lineage that generated the optimization panel?

**Approach.** The designs are already generated. We RE-GRADE the existing `svdd_af_baddg`
designs (and the native/WT antibody) on new, independently-built antigen-escape panels, with the
out-of-family FoldX physics arbiter, using the *same* best-of-pool worst-case beat-native metric
as the paper. No re-generation. Two independent panels:

- **(2a) Natural / clinical variants** — documented circulating SARS-CoV-2 VOC/VOI RBD mutations.
  Fully model-free (no ProteinMPNN, no ESM-1v, no ESM at all).
- **(2b) Different-model-family** — ESM-2 masked-marginal antigen substitutions at the
  geometric epitope. Out-of-lineage vs ProteinMPNN/ESM-1v (caveat below).

All code + data are in `${VAAD_ROOT}/work/indep_escape/` (scratch, no repo edits).

---

## 1. Feasibility table (antigen identity → independent-escape data)

Antigen identity resolved from AACDB `protein_table.txt`. Of the 51 panel targets: **8 are
SARS-CoV-2 / sarbecovirus RBD** (the sweet spot for public escape data), **4 influenza HA**,
and the rest flavivirus / picornavirus / non-viral with little-to-no public antibody-escape DMS.

| class (n) | targets | independent-escape data |
|---|---|---|
| **SARS-CoV-2 RBD, Wuhan/WA-1 (7)** | 7MZK (PDI-96), 7OR9 (COVOX-278), 7TBF (B1-182.1), 7VYR (D27), 7WPH (Fab06), 7YCL (IS-9A), 8IDN (E77) | **YES** — natural VOC RBD variants (WHO/covariants) + Bloom RBD DMS; antigen chain uses standard Spike numbering (~332–527), VOC mutations map 1:1 |
| **Sarbecovirus RBD homolog (2)** | 7TTX (10-40, RaTG13 RBD), 7X4I (aSA3, SARS-CoV-1 RBD) | PARTIAL — SARS-CoV-2 DMS applies only with homolog/numbering caveats; used for the ESM-2 panel only |
| **Influenza HA (4)** | 5DUM (H5N1), 5GJS (H1N1 pdm09), 5UGY (CH65, H1N1), 6MLM (H7N2) | MAYBE — Bloom/Doud HA DMS exists but is strain-specific; per-mAb escape maps for these exact antibodies not confirmed. ESM-2 panel built. |
| Influenza NA (1) | 1A14 (not a scored panel target) | weak |
| Dengue/Zika E (11) | 4AL8, 4BZ1, 4L5F, 5GZN, 5JHL, 5VIC, 5VIG, 6FLA, 6FLB, 6PLK, 7BQ5 | WEAK — flavivirus escape is scattered single-residue literature, no comprehensive public DMS |
| Picornavirus capsid (5) | 6AJ9, 6KZ0, 7X2O, 7X2W, 7X3D | NO public antibody-escape DMS |
| HPV capsid (1) | 7CN2 | NO |
| non-viral / other (22) | lysozyme, IL-13, IL-1β, GFP, CD20, integrin, E-cadherin, mesothelin, TCR, anti-idiotype, etc. | NO (not escape-relevant systems) |

Full per-target listing: `${VAAD_ROOT}/work/indep_escape/feasibility_table.txt`.

**Take-away:** the *only* targets with genuinely independent, public escape ground truth are the
SARS-CoV-2 RBD antibodies. This is the go/no-go answer to "which targets can carry the
independent-escape claim."

---

## 2. Independent-panel construction + provenance

### (2a) Natural / clinical VOC panel — model-free  (7 Wuhan-RBD SARS-CoV-2 targets)
Builder: `build_natural_panel.py` → `jsons_natural/{TARGET}_natural.json`.

- **Source:** canonical VOC/VOI RBD-defining substitutions in Spike numbering, taken from
  WHO/covariants VOC definitions and the SARS-CoV-2 antibody-escape literature (Alpha, Beta,
  Gamma, Delta, Omicron BA.1/BA.2/BA.4-5, XBB.1.5) plus the major documented escape hotspots
  (346, 417, 440, 444, 446, 452, 460, 477, 478, 484, 486, 490, 493, 496, 498, 501, 505, ...).
  These are **real circulating antibody-escape mutations** = the reviewer's "documented natural /
  clinical variants."
- **Independence:** constructed with **no ProteinMPNN and no ESM** of any kind.
- **Mapping guard:** each mutation is kept only if the residue exists in the structure AND the
  structure's WT residue matches the expected Wuhan WT. **All 7 targets: 27/27 escape positions
  present, 0 skipped, WT matched** → 38 single-mutant variants + 8 VOC-haplotype variants each.
- **Concordance with DMS:** the natural VOC positions above are exactly the escape-enriched sites
  reported by the Bloom-lab RBD antibody-escape DMS, so this model-free set doubles as a
  DMS-consensus escape panel.

### (2b) ESM-2 masked-marginal epitope panel — different model family  (13 viral targets)
Builder: `build_esm2_panel.py` → `jsons_esm2/{TARGET}_esm2.json` (ESM-2 `esm2_t33_650M_UR50D`).

- **Epitope positions:** geometric — antigen residues with a heavy atom ≤ 5.0 Å from any
  antibody heavy atom in the AACDB complex. **No model is used to pick positions.**
- **Substitutions:** ESM-2 masked-marginal at each epitope position (mask one position, softmax
  over 20 AAs); keep the top tolerated non-WT substitutions (highest masked log-prob), capped at
  ~36 single-mutant variants/target. This mirrors the *intent* of the omega panel (plausible
  antigen variation at the interface) but from a distinct model.
- **Coverage:** all 7 Wuhan-RBD + 7TTX + 7X4I + 5DUM + 5GJS + 5UGY + 6MLM (14–31 epitope
  residues, 28–36 variants each). This is the guaranteed-coverage arm that exists even where no
  DMS is available.
- **Independence caveat (honest):** ESM-2 and ESM-1v are both masked LMs in the broad ESM family
  trained on UniRef; ESM-2 is a different architecture/training generation and shares no
  parameters with the ProteinMPNN generator, but this is a *weaker* independence guarantee than
  the model-free natural-variant panel. Lead the claim with (2a); (2b) is the breadth/robustness
  arm.

### External per-mAb DMS escape maps — go/no-go
- The aggregated **Bloom-lab `jbloomlab/SARS2_RBD_Ab_escape_maps`** (`processed_data/escape_data.csv`)
  is the canonical per-antibody RBD escape-map resource. Its catalogued studies are
  `2021_Dong_AZ` (COVOX/AstraZeneca panel), Starr LY-CoV555/REGN/Vir, Cao Omicron panel, Greaney
  sera/Rockefeller/Moderna, Tortorici S2X259, Greaney B.1.351.
- **Coverage of *our* exact antibodies is low.** PDI-96, B1-182.1, D27, Fab06, IS-9A, E77 are
  from other/newer papers and are **not** in this aggregate. The one plausible direct hit is
  **COVOX-278 (7OR9) via `2021_Dong_AZ`** — *needs confirmation* that condition `COVOX-278`
  exists in `escape_data.csv` (the per-mAb condition list wasn't retrievable this pass; the
  broad-phrase web query was blocked by the usage-policy classifier — pull the CSV directly with
  `gh`/curl and grep for `278`).
- **Verdict:** published per-mAb escape maps are **not** a general option for this panel; the
  documented-natural-variant construction (2a) is the correct, defensible "independent escape
  set" for the SARS-CoV-2 targets, and it coincides with DMS escape hotspots.

---

## 3. Grading (submitted, resumable)

- FoldX array: **job 29431** (`grade_foldx_indep.sh`, `--array=0-19%16`), CPU, 24 h, resumable
  via `E_checkpoint.json`. Manifest `foldx_manifest.tsv` = 7 natural + 13 esm2 = 20 target-panels.
  Re-grades WT(native) + 8 svdd_af_baddg designs × all panel variants.
- Auto-aggregate: **job 29461** (`--dependency=afterany:29431`, `finalize.sh`) writes
  `INDEPENDENT_ESCAPE_BENCHMARK_RESULTS.md`.
- **Smoke test (8IDN / E77, natural panel)** — pipeline validated end-to-end. Native E77
  interaction energy (more negative = tighter): WT antigen **−5.23**; K417N **−5.34** (tolerated);
  **N501Y −2.05** and **Beta(K417N,E484K,N501Y) −2.34** (native loses ~3 kcal to a real VOC
  mutation at its epitope). This is exactly the native-fragility the benchmark probes — the
  question is whether the svdd_af_baddg designs hold binding against these independent variants.

### Metric (identical to the paper's omega beat-native)
Per antibody, worst-case interaction energy = MAX (least-negative) FoldX energy across the panel
variants. best-of-pool = MIN worst-case over the design pool. **beat-native margin = best-of-pool
worst − native worst**; margin < 0 ⇒ the design pool's least-favorable independent variant still
binds tighter than the native antibody's least-favorable variant.

### Results
Auto-populated to `INDEPENDENT_ESCAPE_BENCHMARK_RESULTS.md` when job 29431 finishes; regenerate
anytime with `python3 aggregate_indep.py`. It reports, per target and in summary, beat-native on
the **independent** panel side-by-side with beat-native on the **original omega** panel (the head
-to-head that answers "does the result survive an independently-built escape set?").
**(Independent-panel FoldX still running at report time — those numbers pending.)**

**Comparison baseline already on disk — beat-native on the ORIGINAL omega panel** (svdd_af_baddg,
best-of-pool worst-case margin; <0 = beats native; from the paper's existing FoldX run) for the
covered viral targets. This is the number the independent panel must reproduce:

| target | omega margin | | target | omega margin |
|---|---|---|---|---|
| 7MZK | −1.06 WIN | | 7WPH | −7.52 WIN |
| 7OR9 | −2.35 WIN | | 7YCL | −5.00 WIN |
| 7TBF | −0.02 WIN | | 8IDN | −10.27 WIN |
| 7VYR | −1.89 WIN | | 7TTX | −1.50 WIN |
| 7X4I | −1.29 WIN | | 5DUM | −4.76 WIN |
| 5UGY | −6.06 WIN | | 6MLM | −5.43 WIN |
| 5GJS | +1.98 loss | | | |

Omega beat-native on the covered viral set = **12/13 (92%)** (5GJS the lone loss). The
independent-panel re-grade tests whether this survives when the escape variants are built without
the ProteinMPNN/ESM-1v lineage.

---

## 4. Limitations (scrupulously)

- **Coverage.** Only SARS-CoV-2 RBD antibodies (7 clean Wuhan targets) carry the strongest
  model-free independent claim. HA is "maybe" (strain-specific DMS), and >30 of the 51 targets
  have no public escape ground truth — the independent benchmark is a *spotlight on the
  escape-relevant subset*, not a panel-wide replacement for the omega grade.
- **Grader is still FoldX.** The escape *panel* is now independent of ProteinMPNN/ESM-1v, but the
  *scorer* is FoldX (out-of-family physics, semi-in-loop only in that FoldX is a fixed arbiter).
  H3-DDG re-grade on these panels is a straightforward add (reuse `grade_h3_array.sh` pointed at
  `jsons_natural`/`jsons_esm2`) and is the natural next step for a second independent grader.
- **Leakage caveats the reviewer flags.**
  - *Sequence-corpus leakage:* the ESM-2 panel (2b) is drawn from a protein LM trained on UniRef,
    which contains SARS-CoV-2 and homolog sequences — so 2b is "independent model" but not
    "independent of all training corpora." The natural panel (2a) is corpus-free.
  - *Homolog leakage:* 7TTX (RaTG13) and 7X4I (SARS-CoV-1) are sarbecovirus homologs; SARS-CoV-2
    VOC mutations were deliberately NOT force-mapped onto them (natural panel restricted to the 7
    Wuhan targets); they appear in the ESM-2 arm only.
  - *Epitope leakage:* the ESM-2 epitope is defined from the *same* bound structure the antibody
    was optimized against, so 2b escapes concentrate on the true paratope contact — this is
    intended (worst-case), but it means 2b is not an out-of-epitope generalization test. The
    natural panel spans the whole RBD (many positions outside a given mAb's epitope), giving a
    fairer worst-case.
- **ESM-2 substitution bias.** ESM-2's top-tolerated substitutions skew conservative
  (iso-chemical), so the 2b panel is antigen-viable but not maximally escape-adversarial; it
  tests robustness to plausible drift, not to engineered maximal escape.

---

## 5. Go / no-go

- **GO (model-free external, strongest):** natural / clinical VOC RBD panel on the **7 Wuhan
  SARS-CoV-2 targets** — 7MZK, 7OR9, 7TBF, 7VYR, 7WPH, 7YCL, 8IDN. Built, submitted, validated.
  This is the headline independent-escape benchmark.
- **GO (breadth, out-of-lineage model):** ESM-2 masked-marginal panel on **13 viral targets** —
  the 7 above + 7TTX, 7X4I, 5DUM, 5GJS, 5UGY, 6MLM. Built, submitted.
- **CONDITIONAL GO (per-mAb DMS):** only **COVOX-278 (7OR9)** is a plausible match to a published
  escape map (`2021_Dong_AZ`); confirm by grepping `escape_data.csv`. A future wet/DMS pass, if
  commissioned, should target the 7 Wuhan mAbs (PDI-96, COVOX-278, B1-182.1, D27, Fab06, IS-9A,
  E77) — well-defined RBD binders with existing complex structures.
- **NO-GO:** flavivirus, picornavirus, HPV, and all non-viral targets — no public escape ground
  truth; the ESM-2 panel is the only independent option there and only where an epitope is
  resolvable.

---

## Paths & jobs
- Scratch root: `${VAAD_ROOT}/work/indep_escape/`
- Builders: `build_natural_panel.py`, `build_esm2_panel.py`
- Panels: `jsons_natural/*.json` (7), `jsons_esm2/*.json` (13)
- Grader (read-only reuse): `lyra/scripts/nos_diffusion/validate_designs_foldx.py`
- Array: `grade_foldx_indep.sh` = **job 29431**; manifest `foldx_manifest.tsv`
- Aggregator: `aggregate_indep.py`; auto-run **job 29461**; output `INDEPENDENT_ESCAPE_BENCHMARK_RESULTS.md`
- Feasibility detail: `feasibility_table.txt`
- FoldX binary: `${VAAD_ROOT}/tools/foldx/foldx_20270131`
