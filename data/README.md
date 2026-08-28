# Data

The benchmark **definition** lives here. Raw structures and the full candidate variant library
do not — see "Not in this repository" below.

| path | rows | what |
|---|---:|---|
| `manifests/robustbench_targets.csv` | 200 | one row per complex: chains, CDR-H3 numbering, native CDR-H3, cluster labels, protein family, library size |
| `manifests/variant_panels.csv` | 6,200 | every (target, variant): `split` ∈ {steer, eval, both}, omega, mutation string |
| `manifests/antigen_clusters.csv` | 200 | **both** cluster labellings — see the 171/172 note below |
| `manifests/postcutoff_targets.csv` | 10 | post-2024-06, low-identity antigens |
| `manifests/natural_variants.csv` | 322 | documented natural SARS-CoV-2 RBD variants, 7 targets |
| `splits/design_variants.csv` | 1,200 | the 6 design-time variants per target |
| `splits/heldout_variants.csv` | 4,800 | the 24 held-out variants per target |
| `sequences/antigens.fasta` | 200 | antigen sequences (extraction rule in `configs/benchmark.yaml`) |
| `sequences/vh.fasta` | 200 | antibody VH sequences |
| `processed/shared_1785_design_pool.csv` | 8,925 | 1,785 designs × 5 evaluators, for the agreement analysis |

Per target: 1 WT (`split=both`) + 6 design + 24 held-out = 31 antigen states. Disjointness is
enforced by construction and checked in `tests/test_variant_split.py`.

`mutation_string` is `CHAIN:WT_RESNUM_MUT` joined by `;` (e.g. `A:E366D;A:H367W`), in the PDB
residue numbering of the source structure.

## The 171 vs 172 cluster counts — both are correct

They answer different questions and neither should be substituted for the other.

* **172** — the fixed benchmark cluster labels. Obtained by clustering the **1,126 structurally
  eligible candidate** antigens and counting distinct labels among the 200 selected targets.
  These are the labels used as the resampling unit for every cluster-bootstrap interval in the
  paper. Column: `benchmark_cluster`.
* **171** — the number of clusters in the **full AACDB** clustering (7,491 sequences → 1,484
  clusters) that contain at least one benchmark target. This is the coverage statistic:
  171/1,484 = 11.5% of AACDB antigen diversity. Column: `full_aacdb_cluster`.

The difference is a property of greedy set-cover (`--cluster-mode 0`): cluster assignment
depends on the whole input set. One target, `7JN5_HLF`, had its own cluster in the candidate-set
clustering while ten other SARS-CoV-2 RBD targets sat under a representative that is not itself
among the 200; clustered in isolation the eleven merge. Both labellings ship so either statistic
can be recomputed.

Clustering parameters are identical in both cases:
`mmseqs easy-cluster --min-seq-id 0.7 -c 0.8 --cov-mode 0 --cluster-mode 0`.

## Provenance and limits

Structures come from **AACDB**. Antigen variants are **model-generated, not measured**:
ProteinMPNN proposes substitutions at epitope positions (4.5 Å heavy-atom contact; 8.0 Å shell),
beam search combines them to depth 3, ESM-1v filters for plausibility.

- The held-out set is **86% triple mutants**; only 47/200 targets contain any single mutant.
- Held-out performance correlates with wild-type binding at **r = +0.72**, so the benchmark does
  not cleanly separate variant-specific robustness from plain affinity.
- Agreement with measured deep-mutational-scanning escape is moderate: Spearman ≈ 0.33, Jaccard
  ≈ 0.38 on epitope positions, across 28 SARS-CoV-2 antibody–RBD DMS datasets.

## Not in this repository

**Full candidate variant library** — ≈3,040,707 accepted variants, 1.57 GB raw / 85 MB gzipped
(sha256 `d795810fe3da01ef4a49277fc9c76a3685c030fc65d5c258a1c0221fa187e69a`). Too large for git
and not required to use the benchmark. Regenerate with the pipeline in `src/benchmark/`
(`generate_variants.py` → `esm1v_filter.py` → `build_panels.py`); regeneration is **not**
bit-identical because proposal sampling is unseeded.

**AACDB structures** — obtain from AACDB directly; not redistributed. Expected at
`$VAAD_ROOT/datasets/AACDB/complex_structure/`.

Manifests derived here are released under CC BY 4.0; cite AACDB alongside this work.
