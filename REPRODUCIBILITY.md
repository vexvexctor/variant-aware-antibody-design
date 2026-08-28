# Reproducibility

What can be rebuilt from this repository, what needs a cluster and staged dependencies, and what
is genuinely missing. Nothing in this file is aspirational — every "yes" was executed while
assembling the repository.

Two entry points:

```bash
bash scripts/reproduce_all.sh --from-precomputed   # statistics + figures from released CSVs
bash scripts/reproduce_all.sh --full               # generation + scoring; multi-GPU, days
```

## Claim → code → data → output

| Paper claim / artifact | Script | Input | Output | Reproducible from release? |
|---|---|---|---|---|
| RobustBench: 7,491 antigen sequences, 1,484 full-AACDB clusters, 171 represented, 172 benchmark labels | `src/benchmark/extract_antigens.py`, `src/benchmark/run_aacdb_mmseqs.sh` | AACDB structures | `data/manifests/antigen_clusters.csv` | **Manifests yes**; re-clustering needs AACDB + mmseqs2 |
| 200-target panel selection | `src/benchmark/select_targets.py`, `scan_candidates.py` | candidate TSV + escape libraries | `data/manifests/robustbench_targets.csv` | Manifest yes; rebuild needs raw libraries |
| 6 design / 24 held-out split, disjointness | `src/benchmark/build_panels.py` | escape libraries | `data/splits/*.csv` | **Yes** — verified by `tests/test_variant_split.py` |
| ≈3.04 M candidate antigen variants | `src/benchmark/generate_variants.py`, `esm1v_filter.py` | AACDB + ProteinMPNN + ESM-1v | full library (not in git) | **No** — see "Missing" below |
| CVaR objective, sign and tail | `src/steering/cvar.py` | — | — | **Yes** — `tests/test_cvar.py`, incl. equivalence with the in-loop copy |
| Variant-aware decoding | `src/steering/variant_aware.py`, `decoding.py` | staged generators + BA-DDG | design JSONs | Needs GPU + staged deps |
| Held-out evaluation (primary endpoint) | `src/evaluators/h3_ddg.py` via `scripts/04_evaluate_heldout.sh` | design JSONs | grade JSONs | Needs GPU + H3-DDG checkpoint |
| Main generator comparison (+6.3 to +22.7 pp) | `src/analysis/aggregate_matrix.py` | grade JSONs | `results/headline/master_matrix_cells.csv`, `MASTER_MATRIX.md` | **Aggregates yes**; regrading needs GPU |
| Beat-native metric + cluster bootstrap | `src/analysis/metrics.py`, `bootstrap.py` | per-target margins | rates + CIs | **Yes** — `tests/test_metrics.py` |
| Temperature-matched control | `src/analysis/temperature_control.py` | graded designs at T=0.2/1.0 | `results/controls/TEMPERATURE_CONTROL.md` | Summary yes; regeneration needs GPU |
| Edit-distance-matched null | `src/analysis/edit_matched_null.py` | design pools | `results/controls/lift_over_random_null.md` | Summary yes |
| Objective ablation (CVaR vs mean vs worst) | `configs/steering.yaml` + `aggregate_matrix.py` | graded arms | `results/controls/EXP4_CVAR_MATRIX.md` | Summary yes |
| Panel-size sensitivity M = 1/3/6/12 | `src/analysis/panel_size.py`, `paired_test.py` | 612 design runs + grades | `results/controls/panel_size_summary.csv`, `panel_size_per_unit.csv` | **Yes from per-unit CSV** |
| Budget-matched generate-and-rank | `src/analysis/budget_matched_search.py` | shared candidate pool | `results/controls/EXP1_BUDGET.md` | Summary yes |
| RFantibody (de-novo backbone, scored on own structure) | `src/generators/rfantibody_prep.py`, `rfantibody_collect.py` | RFantibody install | `results/headline/generator_by_protein_family.csv` (RFantibody row) | Summary yes |
| Pythia-PPI out-of-family check | `src/evaluators/pythia_ppi.py`, `src/analysis/pythia_beatnative.py` | designs + Pythia | `results/generalization/PYTHIA_OUT_OF_FAMILY.txt` | Summary yes |
| Post-cutoff / low-similarity antigens | — (manifest only) | `data/manifests/postcutoff_targets.csv` | `results/generalization/NOVEL_OOD_RESULTS.md`, `EXP10_POSTCUTOFF.md` | Summary yes |
| Natural SARS-CoV-2 RBD variants | `src/analysis/pythia_beatnative.py` + FoldX | `data/manifests/natural_variants.csv` | `results/generalization/INDEPENDENT_ESCAPE_BENCHMARK_RESULTS.md` | Summary yes |
| **Evaluator agreement (Spearman + Cohen's κ + base rates)** | `src/analysis/evaluator_agreement.py` | `data/processed/shared_1785_design_pool.csv` | `results/evaluator_agreement/{spearman,cohen_kappa,base_rates,percent_agreement,pairwise_long}.csv` | **Yes — fully, in seconds** |
| WT affinity vs variant retention | `src/analysis/affinity_retention.py`, `affinity_retention_ancova.py` | graded designs | `results/headline/EXP5_WT_RETENTION.md`, `wt_retention_quadrants.csv` | **Aggregates yes** |
| Generator × protein family | `src/analysis/plots/plot_generator_by_family.py` | `results/headline/generator_by_protein_family.csv` | `figures/fig_generator_by_family.png` | **Yes** |
| Trajectory analysis (descriptive) | — (raw traces released) | `results/trajectories/*.json` | 121 guided/unguided traces | **Partial** — see below |
| Residue composition, tyrosine/aromatic enrichment | `src/analysis/residue_features.py` | graded sequences | `results/biology/` | Summary yes |
| Epitope-conditioned charge | `src/analysis/epitope_charge.py` | graded sequences | `results/biology/charge_complementarity.txt` | Summary yes |
| Compute cost | `src/analysis/compute_cost.py` | per-run `.cost.json` counters | `results/compute/overhead_by_arm.csv`, `OVERHEAD.md` | Summary yes |

## Missing from the current archive

Stated plainly rather than reconstructed.

1. **The manuscript itself.** No `.tex`, `.bib` or `main.pdf` exists on the machine this
   repository was assembled from. `paper/` therefore contains only a README. The manuscript's
   Figure 1 (method schematic) and several appendix figures (trajectory, biology, target-level
   native figure) are consequently absent; the pseudocode in `README.md` documents the method in
   their place.
2. **The full ≈3.04 M variant candidate library** (1.57 GB raw / 85 MB gzipped). Not in git.
   `data/README.md` records the archive checksum and the code path that regenerates it.
   Regeneration is not bit-identical because proposal sampling is unseeded.
3. **Raw design and grading output** (~77 GB of per-design JSONs and grading shards). Only tidy
   aggregates are released. Absolute rates can be recomputed from `results/headline/`, but
   per-design re-analysis beyond the released columns is not possible from this repository alone.
4. **Trajectory capture is incomplete** — 121 of a planned 120 runs across 12 targets × 5
   generators × guided/unguided, i.e. roughly half the intended grid. The trajectory analysis is
   **descriptive only** and the captured subset should not be read as representative of all runs.
5. **Third-party checkpoints and FoldX.** None are redistributable; `scripts/00_download_dependencies.sh`
   documents sources and the expected layout.
6. **`src/generators/abmpnn.py` does not exist as a standalone file.** AbMPNN runs through the
   live-MPNN code path inside `src/steering/variant_aware.py` (`LIVE_MPNN_BASES`), pointing at an
   AbMPNN checkpoint. Recorded in `configs/generators.yaml` rather than invented as a wrapper.

## Determinism and seeds

| Component | Seeded? | Notes |
|---|---|---|
| Held-out panel selection | **Yes** — `--holdout-seed 0` | The 6/24 split is deterministic given the library |
| Cluster bootstrap | **Yes** — `DEFAULT_SEED = 0`, 10,000 draws | `src/analysis/bootstrap.py`; published intervals reproduce exactly |
| Evaluator agreement | Deterministic (no resampling) | Exact matrices reproduce |
| Generator sampling | **No** | Unseeded; run-to-run sd ≈ 0.5–0.7 H3 units. Single-run cells cannot resolve most effects — use ≥3 replicates and compare rep-matched arms |
| Variant proposal / beam search | **No** | Why library regeneration is not bit-identical |

The unseeded sampling is a property of the original runs, not a choice made here. It is why
every headline comparison is replicate-matched and cluster-bootstrapped.

## Verification performed while assembling this repository

- All 38 Python modules under `src/` byte-compile.
- All 27 unit tests pass, including CVaR equivalence with the in-loop implementation in
  `variant_aware.py`, split disjointness across all 200 targets, and the cluster-bootstrap
  resampling unit.
- `scripts/reproduce_all.sh --from-precomputed` runs clean end to end.
- Every figure in `figures/` regenerates from a released CSV; no result value is typed into a
  plotting script.
- Zero occurrences of the original absolute paths or username anywhere in `src/` or `scripts/`.
- Evaluator-agreement matrices regenerated from the shared pool match the manuscript values
  exactly (Spearman .661/.283/.044/.121/.469/.207/.246/.475/.345/.294; κ .429/.170/.067/.004/
  .325/.081/≈0/.317/.023/−.012; base rates 39.7/37.0/33.4/16.5/1.3).
- Panel-size figures regenerate to the manuscript values (70.6/66.7/66.7/64.7%).
