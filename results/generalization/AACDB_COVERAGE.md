# AACDB antigen-diversity coverage of the 200-target benchmark

_Run 2026-08-25._

## Headline

```
mmseqs version:                              cb12a2d75a9808ee61721029d064a7bd80af6fec
Full AACDB sequences clustered:              7491
Full AACDB clusters:                         1484
Benchmark clusters when clustered alone:     171      <-- NOT 172; see below
Full-AACDB clusters represented by our 200:  171
Cluster coverage: 171/1484 = 11.52%
Complex coverage: 200/7491 = 2.67%
```

**Merging check:** the 200 occupy 171 clusters alone and 171 within full AACDB — **zero net
merging**. The anticipated effect (benchmark-only clusters merging once the other ~7,300
antigens are added) did not occur at these thresholds.

## The sanity check returned 171, and it is not a version/header/extraction problem

Per instruction I changed no parameters. The cause is a **definitional** one, and it reproduces
exactly.

`172` was never "cluster the 200 and count clusters." In `select_panel.py`, the antigen
clustering is run over **all 1,126 eligible candidates**, and each panel target inherits the
resulting cluster label; `172` is the number of *distinct labels among the 200 selected targets*.
Reproducing that path:

| computation | result |
|---|---|
| cluster all 1,126 candidates | 388 clusters |
| distinct clusters among the 200 | **172** ✓ |
| exact rep-label agreement with the stored `ag_cluster` column | **200/200** |
| cluster the 200 alone | **171** |

The off-by-one is a single target: **`7JN5_HLF`**. Historically it sat in its own cluster
(representative `7JN5_HLF`) while ten other SARS-CoV-2 RBD targets — 7TTX, 7MZK, 7OR9, 7TBF,
7VYR, 7WPH, 7X4I, 7YCL, 8IDN, 7N4I — sat under representative **`7WVL_HLF`**, a candidate that
is *not* itself in the 200. Remove the sequences that defined that boundary and greedy set-cover
puts all eleven in one cluster. This is expected behaviour for `--cluster-mode 0`: assignments
depend on the whole input set, not just the members.

So both numbers are correct for their own question:
- **172** = distinct antigen clusters among the 200, as labelled by the candidate-set clustering (what the paper's N behind cluster-bootstrap CIs refers to).
- **171** = clusters formed by the 200 in isolation.
- **171** = full-AACDB clusters containing ≥1 benchmark target — the number for the coverage ratio.

## Extraction rule (identical in both FASTAs)

Taken verbatim from `scan_candidates.py`, the historical path:

```
antibody chains = parse_chains(id)      # 'XLY' -> [X,L]; otherwise [X]
antigen chains  = remaining chars of the id suffix
sequence        = parse_antigen_residues(pdb, antigen_chains)   # gap positions kept as 'X'
ag_seq          = sequence[:1200]                                # truncated at 1200
header          = the complex id, e.g. 1IC7_HLY
eligibility     = drop antigens shorter than 30 residues
```

The historical extraction read the antigen block out of the precomputed
`complex_features/*.pt`, which exist for only 4,238 of 7,695 complexes. Reading from the PDB
with `parse_antigen_residues()` — the same parser the features were built from — was validated
against the recorded `ag_seq` for all 1,126 candidates:

> **PDB-extraction == historical feature-extraction: 1126/1126 (byte-identical).**

Gap placeholders must be **kept**: stripping `X` drops the match to 568/1126.

`benchmark_200_antigens.fasta` is an exact subset of `AACDB_all_antigens.fasta` (same rule, same
headers, same sequences), so no cross-file discrepancy is possible.

## Counts, reconciled

| quantity | value | note |
|---|---:|---|
| PDB files in the AACDB snapshot | 7,695 | `data_zip/complex_structure/*.pdb` |
| antigen sequences extracted | **7,491** | 204 dropped as antigen < 30 residues |
| rows in `manifest.csv` | 4,584 | a partial manifest; not the snapshot size |
| precomputed `complex_features/*.pt` | 4,238 | why extraction goes via the PDB |
| eligible candidates (`candidates.tsv`) | 1,126 | escape library + features + CDR-H3 checks |
| benchmark panel | 200 | 172 candidate-set clusters / 171 standalone |

7,491 is within 7 of the 7,498 figure — consistent with a slightly different snapshot or
eligibility cut on that side.

## Commands

```bash
mmseqs easy-cluster AACDB_all_antigens.fasta mmseqs_aacdb_full/aacdb70 mmseqs_aacdb_full/tmp \
  --min-seq-id 0.7 -c 0.8 --cov-mode 0 --cluster-mode 0
mmseqs easy-cluster benchmark_200_antigens.fasta mmseqs_bench_only/bench70 mmseqs_bench_only/tmp \
  --min-seq-id 0.7 -c 0.8 --cov-mode 0 --cluster-mode 0
```

The historical `select_panel.py` call omitted `--cluster-mode 0` (relying on the default, which
is 0) and passed `-v 1`; otherwise identical.

## Suggested B.1 wording

> Our 200-target benchmark spans 171 of the 1,484 antigen clusters formed by clustering all
> 7,491 AACDB antigen sequences at 70% identity and 80% bidirectional coverage
> (`--min-seq-id 0.7 -c 0.8 --cov-mode 0 --cluster-mode 0`), i.e. **11.5% of the antigen
> diversity in AACDB**, while comprising 2.7% of its complexes. Because cluster assignment under
> greedy set-cover depends on the full input set, the panel resolves into 172 clusters when
> labelled against the 1,126 eligible design candidates — the figure used as the effective N for
> cluster-bootstrap confidence intervals — and 171 when clustered in isolation.

## Files

| file | what |
|---|---|
| `AACDB_all_antigens.fasta` | 7,491 antigen sequences, historical rule |
| `benchmark_200_antigens.fasta` | the 200, exact subset |
| `candidates_1126_antigens.fasta` | the 1,126 eligible candidates (reproduces 172) |
| `extract_antigens.py` | extraction + `validate` mode (1126/1126 check) |
| `run_aacdb_mmseqs.sh` | both clusterings + coverage arithmetic |
| `mmseqs_aacdb_full/`, `mmseqs_bench_only/`, `mmseqs_cand/` | mmseqs outputs incl. `*_cluster.tsv` |
