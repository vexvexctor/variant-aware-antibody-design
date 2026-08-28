# EXP2 — clustered antigen-holdout re-analysis (headline generalization)

The headline cell **svdd/antifold/baddg** has no in-house component trained on the leaked
escape-variant split (AntiFold / BA-DDG / H3-DDG frozen external; SVDD training-free; FoldX
physics). So the leakage does not gate it through a model we train. This re-slices the
existing per-target margins by the antigen-cluster holdout structure (no retrain, no new
compute). FoldX = out-of-family physics grader (training-free, the honest number).
`beat-nat` = fraction with best-of-pool worst-case held-out margin < 0. See LEAKAGE_ANALYSIS.md.

## svdd/antifold/baddg

**FoldX** — naive per-target: **41/51 (80.4%)**, mean margin -2.914
  cluster-aware (antigen_cluster_70, effective N=31 clusters):
    - cluster-majority beat-native: **23/31 (74.2%)**
    - cluster-mean-margin < 0:      23/31 (74.2%)

**H3** — naive per-target: **30/51 (58.8%)**, mean margin -0.351
  cluster-aware (antigen_cluster_70, effective N=31 clusters):
    - cluster-majority beat-native: **20/31 (64.5%)**
    - cluster-mean-margin < 0:      19/31 (61.3%)

### svdd/antifold/baddg — beat-native stratified by leakage class (FoldX)

leakage_class    |     beat-nat | mean margin
NOVEL_ANTIGEN    |    1/2     50% |      -0.897
SIMILAR_ANTIGEN  |    1/1    100% |      -1.496
SEEN_ANTIGEN     |    9/13    69% |      -2.347
SEEN_PAIR        |   30/35    86% |      -3.281

## svdd/abmpnn/baddg

**FoldX** — naive per-target: **13/18 (72.2%)**, mean margin -3.956
  cluster-aware (antigen_cluster_70, effective N=11 clusters):
    - cluster-majority beat-native: **8/11 (72.7%)**
    - cluster-mean-margin < 0:      8/11 (72.7%)

**H3** — naive per-target: **39/51 (76.5%)**, mean margin -0.872
  cluster-aware (antigen_cluster_70, effective N=31 clusters):
    - cluster-majority beat-native: **24/31 (77.4%)**
    - cluster-mean-margin < 0:      25/31 (80.6%)

### svdd/abmpnn/baddg — beat-native stratified by leakage class (FoldX)

leakage_class    |     beat-nat | mean margin
NOVEL_ANTIGEN    |    1/1    100% |      -2.724
SIMILAR_ANTIGEN  |           -- |          --
SEEN_ANTIGEN     |    4/6     67% |      -1.916
SEEN_PAIR        |    8/11    73% |      -5.181

## Least-leaked / novel-antigen targets (the held-out-antigen evidence)

Sorted by antigen identity to nearest non-self train complex (lowest = most novel).
These are the targets a reviewer will read as genuine held-out antigens.

target      | class           | ag_pid | h3_pid | antifold FoldX |   abmpnn FoldX
1KB5_HLB    | NOVEL_ANTIGEN   |   25.9 |  100.0 |          -1.80 |          -2.72
6AJ9_HLC    | NOVEL_ANTIGEN   |   49.7 |   69.2 |          +0.01 |             --
7TTX_HLA    | SIMILAR_ANTIGEN |   89.2 |  100.0 |          -1.50 |             --
8CYH_HLM    | SEEN_ANTIGEN    |   92.0 |  100.0 |          +1.83 |             --
4BZ1_HLA    | SEEN_ANTIGEN    |   95.0 |   63.6 |          +1.37 |          +0.62
6VRQ_HLA    | SEEN_ANTIGEN    |   95.7 |   58.3 |          -2.77 |             --
4AL8_HLC    | SEEN_PAIR       |   97.0 |  100.0 |          -0.52 |          -0.28
8IDN_HLA    | SEEN_ANTIGEN    |   98.5 |   69.2 |         -10.27 |             --

(FoldX margin < 0 = the design's worst held-out variant still binds tighter than the native
antibody's worst. `ag_pid`/`h3_pid` = % identity to nearest non-self train antigen / CDR-H3.)
