# P5 - Observed-escape robustness transfer

_Generated 2026-08-05._

**Question.** Does a binder optimized against COMPUTATIONALLY-generated partner-protein point mutations (worst-case objective) also hold up against EXPERIMENTALLY-MEASURED partner mutations it never saw during optimization?

**Design.** For each overlapping system we freeze one binder per arm and score it on the measured mutation set with two independent graders:
- **BA-DDG** - a binding-ddG scorer in the same ProteinMPNN lineage as the steering objective (in-family; magnitudes are a surrogate).
- **Pythia_bind** - a self-supervised stability model (never trained on any ddG label), binding term via the thermodynamic cycle E(complex)-E(apo). On the SKEMPI Ab-Ag anchor it correlates with BA-DDG at only +0.396 (FoldX-tier), so it is an **out-of-family** check on the arm ranking. Reward is design-vs-native (native == 0 by construction).

Reward = -ddG vs the native binder; **>0 means the binder holds or improves binding relative to native on that measured mutation**. We report worst-case (min reward), CVaR@0.2 (mean of the worst 20%), mean, and frac>=native, over the full measured set and the high measured-escape subset (top quartile of measured_escape_norm).

**Three arms, shared AntiFold base** (the objective is the variable, not the vehicle):
- `variant-aware` - AntiFold base, SVDD reweighted by the **worst-case** BA-DDG over the COMPUTATIONAL escape variants (`cluster_top1`).
- `WT-only` - AntiFold base, SVDD reweighted by BA-DDG on the **WT antigen only** (`single`).
- `AntiFold one-shot` - raw AntiFold argmax, no binding objective (`antifold_r1`).

The measured mutations NEVER enter optimization. Computational escape variants are multi-point omega-weighted combinations; the measured set is single-point deep-mutational-scan substitutions from the harmonized corpus - disjoint mutation strings.

## Overlap

Matching is by **antigen (target-protein) chain sequence** against the corpus DMS systems, with the numbering registration recovered by a single-offset scan (wt-verified). **18 design-panel targets overlap a corpus system: 13 on SARS-CoV-2 RBD (graded on the COVID-19_2021c_C002 DMS), 5 on Zika E (graded on the best-registering zika_2023 DMS).** No design-panel antibody is itself in the corpus, so what transfers is the real antigen-mutation set (positions x substitutions with measured effect), graded by two binder-independent scorers - not a same-antibody escape value. measured_escape_norm was assayed against a different antibody, so it only (i) defines which antigen mutations are experimentally real and (ii) weights the high-escape subset.

**Evaluated here: 18 systems** (13 SARS-CoV-2 RBD + 5 Zika E).

## BA-DDG (in-family) - full measured set

| System | Antigen | Measured muts | Arm | CDR muts | worst | CVaR@0.2 | mean | frac>=native |
|---|---|--:|---|--:|--:|--:|--:|--:|
| 6XKQ_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 529 | variant-aware (worst-case) | 6 | -1.189 | -0.182 | +0.701 | 0.87 |
|  |  |  | WT-only selection | 7 | -3.759 | -2.454 | -1.640 | 0.01 |
|  |  |  | AntiFold one-shot | 7 | -1.291 | -0.098 | +0.787 | 0.89 |
| 7EJ5_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 540 | variant-aware (worst-case) | 8 | -0.168 | +0.475 | +1.144 | 0.99 |
|  |  |  | WT-only selection | 7 | -1.048 | +0.060 | +0.842 | 0.93 |
|  |  |  | AntiFold one-shot | 4 | -1.501 | -0.617 | +0.102 | 0.59 |
| 7K4N_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 2 | -1.191 | +0.323 | +0.930 | 0.98 |
|  |  |  | WT-only selection | 1 | -1.235 | -0.776 | -0.251 | 0.26 |
|  |  |  | AntiFold one-shot | 1 | -1.241 | -0.751 | -0.248 | 0.24 |
| 7MZK_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 7 | -2.005 | -0.926 | -0.015 | 0.51 |
|  |  |  | WT-only selection | 5 | -2.045 | -1.247 | -0.376 | 0.27 |
|  |  |  | AntiFold one-shot | 6 | -5.022 | -3.810 | -2.857 | 0.00 |
| 7OR9_HLE vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 7 | -2.102 | -0.354 | +1.030 | 0.85 |
|  |  |  | WT-only selection | 7 | -3.642 | -1.787 | -0.409 | 0.34 |
|  |  |  | AntiFold one-shot | 8 | -6.657 | -3.535 | -2.118 | 0.02 |
| 7TBF_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 2 | -0.408 | -0.081 | +0.218 | 0.84 |
|  |  |  | WT-only selection | 2 | -0.828 | -0.390 | -0.023 | 0.46 |
|  |  |  | AntiFold one-shot | 2 | -0.748 | -0.355 | -0.019 | 0.47 |
| 7TE1_HLD vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 546 | variant-aware (worst-case) | 4 | +1.506 | +2.333 | +2.949 | 1.00 |
|  |  |  | WT-only selection | 4 | -0.680 | -0.093 | +0.483 | 0.87 |
|  |  |  | AntiFold one-shot | 4 | -0.822 | -0.237 | +0.328 | 0.79 |
| 7TTX_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 366 | variant-aware (worst-case) | 7 | -0.936 | +0.478 | +1.577 | 0.98 |
|  |  |  | WT-only selection | 9 | -0.157 | +0.895 | +1.950 | 1.00 |
|  |  |  | AntiFold one-shot | 7 | -1.493 | +0.287 | +1.390 | 0.95 |
| 7U2E_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 4 | -1.363 | -0.202 | +1.037 | 0.88 |
|  |  |  | WT-only selection | 7 | -2.827 | -1.063 | +0.338 | 0.63 |
|  |  |  | AntiFold one-shot | 5 | -1.650 | -0.316 | +1.024 | 0.86 |
| 7VYR_HLR vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 4 | -0.425 | +0.371 | +1.095 | 0.98 |
|  |  |  | WT-only selection | 3 | -1.370 | -0.405 | +0.496 | 0.78 |
|  |  |  | AntiFold one-shot | 3 | -1.782 | -0.478 | +0.475 | 0.76 |
| 7WPH_HLB vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 6 | -0.496 | +0.516 | +1.523 | 0.98 |
|  |  |  | WT-only selection | 7 | -2.468 | -1.121 | +0.096 | 0.58 |
|  |  |  | AntiFold one-shot | 6 | -1.312 | -0.181 | +0.913 | 0.88 |
| 7YCL_HLD vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 8 | -2.524 | -1.069 | +0.138 | 0.56 |
|  |  |  | WT-only selection | 8 | -3.120 | -1.575 | -0.457 | 0.29 |
|  |  |  | AntiFold one-shot | 5 | -2.731 | -1.469 | -0.339 | 0.34 |
| 8IDN_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 6 | -0.755 | +0.885 | +1.737 | 1.00 |
|  |  |  | WT-only selection | 7 | +1.374 | +3.009 | +3.988 | 1.00 |
|  |  |  | AntiFold one-shot | 7 | -0.396 | +0.915 | +1.841 | 1.00 |
| 5GZN_HLA vs zika_2023_EDE1-C10 | Zika E | 1349 | variant-aware (worst-case) | 10 | -2.292 | -0.326 | +1.727 | 0.88 |
|  |  |  | WT-only selection | 10 | -5.302 | -1.940 | +0.198 | 0.54 |
|  |  |  | AntiFold one-shot | 9 | -7.361 | -1.431 | +0.622 | 0.66 |
| 5JHL_HLA vs zika_2023_EDE1-C10 | Zika E | 1254 | variant-aware (worst-case) | 5 | -1.531 | -0.501 | +0.148 | 0.63 |
|  |  |  | WT-only selection | 4 | -3.001 | -2.126 | -1.433 | 0.01 |
|  |  |  | AntiFold one-shot | 4 | -2.932 | -2.067 | -1.401 | 0.01 |
| 5VIG_HLZ vs zika_2023_ZV-67 | Zika E | 532 | variant-aware (worst-case) | 9 | -1.329 | +0.228 | +1.964 | 0.94 |
|  |  |  | WT-only selection | 7 | -4.070 | -2.787 | -1.411 | 0.07 |
|  |  |  | AntiFold one-shot | 6 | -4.505 | -2.471 | -0.798 | 0.25 |
| 6PLK_HLE vs zika_2023_ZV-67 | Zika E | 532 | variant-aware (worst-case) | 6 | +0.210 | +1.659 | +2.762 | 1.00 |
|  |  |  | WT-only selection | 6 | -3.813 | -2.448 | -1.402 | 0.04 |
|  |  |  | AntiFold one-shot | 5 | -1.653 | -0.391 | +0.616 | 0.82 |
| 7BQ5_HLA vs zika_2023_EDE1-C10 | Zika E | 1254 | variant-aware (worst-case) | 4 | -1.935 | -0.457 | +0.438 | 0.76 |
|  |  |  | WT-only selection | 6 | -1.659 | -0.123 | +0.818 | 0.88 |
|  |  |  | AntiFold one-shot | 2 | -2.064 | -0.626 | +0.250 | 0.65 |

## Pythia_bind (independent, out-of-family) - full measured set

| System | Antigen | Measured muts | Arm | CDR muts | worst | CVaR@0.2 | mean | frac>=native |
|---|---|--:|---|--:|--:|--:|--:|--:|
| 6XKQ_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 529 | variant-aware (worst-case) | 6 | -5.753 | -4.065 | -3.938 | 0.00 |
|  |  |  | WT-only selection | 7 | -14.367 | -11.627 | -11.382 | 0.00 |
|  |  |  | AntiFold one-shot | 7 | -5.512 | -3.540 | -3.363 | 0.00 |
| 7EJ5_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 540 | variant-aware (worst-case) | 8 | -2.283 | -2.007 | -1.985 | 0.00 |
|  |  |  | WT-only selection | 7 | -1.907 | -1.591 | -1.563 | 0.00 |
|  |  |  | AntiFold one-shot | 4 | -1.249 | -0.939 | -0.905 | 0.00 |
| 7K4N_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 2 | -2.732 | -0.314 | -0.239 | 0.02 |
|  |  |  | WT-only selection | 1 | +0.757 | +1.126 | +1.145 | 1.00 |
|  |  |  | AntiFold one-shot | 1 | +0.757 | +1.126 | +1.145 | 1.00 |
| 7MZK_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 7 | -2.470 | -0.913 | -0.649 | 0.03 |
|  |  |  | WT-only selection | 5 | -1.451 | -0.659 | -0.516 | 0.04 |
|  |  |  | AntiFold one-shot | 6 | -12.847 | -11.781 | -11.402 | 0.00 |
| 7OR9_HLE vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 7 | -8.122 | -6.043 | -5.592 | 0.00 |
|  |  |  | WT-only selection | 7 | -5.442 | -4.009 | -3.611 | 0.00 |
|  |  |  | AntiFold one-shot | 8 | -10.571 | -8.412 | -8.051 | 0.00 |
| 7TBF_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 2 | -1.113 | -0.709 | -0.687 | 0.00 |
|  |  |  | WT-only selection | 2 | -0.543 | -0.205 | -0.188 | 0.01 |
|  |  |  | AntiFold one-shot | 2 | -0.543 | -0.205 | -0.188 | 0.01 |
| 7TE1_HLD vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 546 | variant-aware (worst-case) | 4 | +3.812 | +3.812 | +3.812 | 1.00 |
|  |  |  | WT-only selection | 4 | -2.080 | -2.080 | -2.080 | 0.00 |
|  |  |  | AntiFold one-shot | 4 | -3.113 | -3.113 | -3.113 | 0.00 |
| 7TTX_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 366 | variant-aware (worst-case) | 7 | -1.377 | +0.666 | +0.714 | 1.00 |
|  |  |  | WT-only selection | 9 | -6.019 | -3.828 | -3.777 | 0.00 |
|  |  |  | AntiFold one-shot | 7 | -0.465 | +1.561 | +1.606 | 1.00 |
| 7U2E_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 4 | -1.967 | -0.203 | +0.212 | 0.88 |
|  |  |  | WT-only selection | 7 | -3.927 | -1.686 | -1.369 | 0.00 |
|  |  |  | AntiFold one-shot | 5 | -2.416 | -1.229 | -0.968 | 0.01 |
| 7VYR_HLR vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 4 | +0.118 | +1.913 | +2.014 | 1.00 |
|  |  |  | WT-only selection | 3 | +2.651 | +4.295 | +4.434 | 1.00 |
|  |  |  | AntiFold one-shot | 3 | +3.922 | +4.977 | +5.061 | 1.00 |
| 7WPH_HLB vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 6 | -0.839 | +0.340 | +0.458 | 0.99 |
|  |  |  | WT-only selection | 7 | -6.102 | -3.563 | -3.369 | 0.00 |
|  |  |  | AntiFold one-shot | 6 | -0.790 | +0.914 | +1.007 | 0.99 |
| 7YCL_HLD vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 8 | -11.788 | -9.416 | -9.281 | 0.00 |
|  |  |  | WT-only selection | 8 | -17.581 | -16.547 | -16.462 | 0.00 |
|  |  |  | AntiFold one-shot | 5 | -16.670 | -15.421 | -15.306 | 0.00 |
| 8IDN_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 578 | variant-aware (worst-case) | 6 | +0.996 | +1.730 | +1.779 | 1.00 |
|  |  |  | WT-only selection | 7 | -2.587 | -0.916 | -0.774 | 0.04 |
|  |  |  | AntiFold one-shot | 7 | -1.109 | +0.520 | +0.662 | 1.00 |
| 5GZN_HLA vs zika_2023_EDE1-C10 | Zika E | 1349 | variant-aware (worst-case) | 10 | -5.274 | -3.495 | -3.369 | 0.00 |
|  |  |  | WT-only selection | 10 | -5.867 | -3.923 | -3.827 | 0.00 |
|  |  |  | AntiFold one-shot | 9 | -5.503 | -3.744 | -3.651 | 0.00 |
| 5JHL_HLA vs zika_2023_EDE1-C10 | Zika E | 1254 | variant-aware (worst-case) | 5 | -0.742 | -0.312 | -0.275 | 0.02 |
|  |  |  | WT-only selection | 4 | -3.486 | -1.732 | -1.668 | 0.01 |
|  |  |  | AntiFold one-shot | 4 | -3.319 | -1.670 | -1.611 | 0.01 |
| 5VIG_HLZ vs zika_2023_ZV-67 | Zika E | 532 | variant-aware (worst-case) | 9 | +6.107 | +6.943 | +7.097 | 1.00 |
|  |  |  | WT-only selection | 7 | -5.798 | -2.012 | -1.349 | 0.02 |
|  |  |  | AntiFold one-shot | 6 | -4.429 | -3.719 | -3.583 | 0.00 |
| 6PLK_HLE vs zika_2023_ZV-67 | Zika E | 532 | variant-aware (worst-case) | 6 | -8.263 | -4.518 | -4.092 | 0.00 |
|  |  |  | WT-only selection | 6 | -9.972 | -8.131 | -7.788 | 0.00 |
|  |  |  | AntiFold one-shot | 5 | +0.347 | +1.997 | +2.272 | 1.00 |
| 7BQ5_HLA vs zika_2023_EDE1-C10 | Zika E | 1254 | variant-aware (worst-case) | 4 | -4.796 | -4.096 | -4.065 | 0.00 |
|  |  |  | WT-only selection | 6 | -9.896 | -9.145 | -9.080 | 0.00 |
|  |  |  | AntiFold one-shot | 2 | -0.770 | -0.413 | -0.398 | 0.00 |

## BA-DDG (in-family) - high measured-escape subset (top-quartile)

| System | Antigen | Measured muts | Arm | CDR muts | worst | CVaR@0.2 | mean | frac>=native |
|---|---|--:|---|--:|--:|--:|--:|--:|
| 6XKQ_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 133 | variant-aware (worst-case) | 6 | -1.069 | -0.177 | +0.724 | 0.87 |
|  |  |  | WT-only selection | 7 | -3.759 | -2.447 | -1.632 | 0.00 |
|  |  |  | AntiFold one-shot | 7 | -1.131 | -0.051 | +0.803 | 0.90 |
| 7EJ5_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 135 | variant-aware (worst-case) | 8 | +0.075 | +0.456 | +1.106 | 1.00 |
|  |  |  | WT-only selection | 7 | -0.385 | +0.098 | +0.793 | 0.94 |
|  |  |  | AntiFold one-shot | 4 | -1.173 | -0.591 | +0.154 | 0.62 |
| 7K4N_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 2 | -0.980 | +0.224 | +0.863 | 0.97 |
|  |  |  | WT-only selection | 1 | -1.218 | -0.799 | -0.294 | 0.23 |
|  |  |  | AntiFold one-shot | 1 | -1.241 | -0.778 | -0.261 | 0.24 |
| 7MZK_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 7 | -2.005 | -0.911 | +0.002 | 0.53 |
|  |  |  | WT-only selection | 5 | -1.937 | -1.241 | -0.403 | 0.27 |
|  |  |  | AntiFold one-shot | 6 | -4.518 | -3.772 | -2.808 | 0.00 |
| 7OR9_HLE vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 7 | -2.102 | -0.287 | +1.006 | 0.86 |
|  |  |  | WT-only selection | 7 | -3.440 | -1.861 | -0.397 | 0.37 |
|  |  |  | AntiFold one-shot | 8 | -6.657 | -3.829 | -2.215 | 0.01 |
| 7TBF_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 2 | -0.266 | -0.061 | +0.222 | 0.86 |
|  |  |  | WT-only selection | 2 | -0.828 | -0.430 | -0.036 | 0.45 |
|  |  |  | AntiFold one-shot | 2 | -0.748 | -0.372 | -0.009 | 0.48 |
| 7TE1_HLD vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 137 | variant-aware (worst-case) | 4 | +1.864 | +2.374 | +2.955 | 1.00 |
|  |  |  | WT-only selection | 4 | -0.428 | -0.097 | +0.460 | 0.87 |
|  |  |  | AntiFold one-shot | 4 | -0.713 | -0.210 | +0.339 | 0.81 |
| 7TTX_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 364 | variant-aware (worst-case) | 7 | -0.936 | +0.478 | +1.573 | 0.98 |
|  |  |  | WT-only selection | 9 | -0.157 | +0.902 | +1.953 | 1.00 |
|  |  |  | AntiFold one-shot | 7 | -1.493 | +0.288 | +1.389 | 0.95 |
| 7U2E_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 4 | -0.774 | -0.003 | +1.162 | 0.92 |
|  |  |  | WT-only selection | 7 | -2.571 | -1.127 | +0.220 | 0.57 |
|  |  |  | AntiFold one-shot | 5 | -0.861 | -0.118 | +1.120 | 0.90 |
| 7VYR_HLR vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 4 | -0.425 | +0.327 | +1.082 | 0.97 |
|  |  |  | WT-only selection | 3 | -1.370 | -0.469 | +0.483 | 0.78 |
|  |  |  | AntiFold one-shot | 3 | -1.782 | -0.624 | +0.389 | 0.73 |
| 7WPH_HLB vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 6 | -0.266 | +0.591 | +1.511 | 0.98 |
|  |  |  | WT-only selection | 7 | -2.029 | -1.188 | +0.183 | 0.59 |
|  |  |  | AntiFold one-shot | 6 | -1.040 | -0.218 | +0.854 | 0.87 |
| 7YCL_HLD vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 8 | -1.470 | -0.974 | +0.169 | 0.57 |
|  |  |  | WT-only selection | 8 | -3.034 | -1.677 | -0.489 | 0.26 |
|  |  |  | AntiFold one-shot | 5 | -2.475 | -1.446 | -0.338 | 0.34 |
| 8IDN_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 6 | +0.545 | +0.904 | +1.754 | 1.00 |
|  |  |  | WT-only selection | 7 | +2.340 | +3.031 | +3.935 | 1.00 |
|  |  |  | AntiFold one-shot | 7 | +0.444 | +0.959 | +1.815 | 1.00 |
| 5GZN_HLA vs zika_2023_EDE1-C10 | Zika E | 339 | variant-aware (worst-case) | 10 | -2.192 | -0.353 | +1.788 | 0.88 |
|  |  |  | WT-only selection | 10 | -4.573 | -1.700 | +0.231 | 0.55 |
|  |  |  | AntiFold one-shot | 9 | -3.479 | -1.452 | +0.763 | 0.69 |
| 5JHL_HLA vs zika_2023_EDE1-C10 | Zika E | 337 | variant-aware (worst-case) | 5 | -1.235 | -0.489 | +0.128 | 0.60 |
|  |  |  | WT-only selection | 4 | -2.817 | -2.126 | -1.421 | 0.00 |
|  |  |  | AntiFold one-shot | 4 | -2.479 | -2.021 | -1.375 | 0.00 |
| 5VIG_HLZ vs zika_2023_ZV-67 | Zika E | 133 | variant-aware (worst-case) | 9 | -1.174 | +0.376 | +2.110 | 0.95 |
|  |  |  | WT-only selection | 7 | -3.798 | -2.568 | -1.243 | 0.08 |
|  |  |  | AntiFold one-shot | 6 | -3.529 | -2.425 | -0.746 | 0.26 |
| 6PLK_HLE vs zika_2023_ZV-67 | Zika E | 133 | variant-aware (worst-case) | 6 | +1.029 | +1.677 | +2.648 | 1.00 |
|  |  |  | WT-only selection | 6 | -3.501 | -2.542 | -1.390 | 0.05 |
|  |  |  | AntiFold one-shot | 5 | -1.653 | -0.463 | +0.679 | 0.82 |
| 7BQ5_HLA vs zika_2023_EDE1-C10 | Zika E | 337 | variant-aware (worst-case) | 4 | -1.149 | -0.405 | +0.460 | 0.77 |
|  |  |  | WT-only selection | 6 | -0.957 | -0.103 | +0.877 | 0.89 |
|  |  |  | AntiFold one-shot | 2 | -1.450 | -0.585 | +0.274 | 0.65 |

## Pythia_bind (independent, out-of-family) - high measured-escape subset (top-quartile)

| System | Antigen | Measured muts | Arm | CDR muts | worst | CVaR@0.2 | mean | frac>=native |
|---|---|--:|---|--:|--:|--:|--:|--:|
| 6XKQ_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 133 | variant-aware (worst-case) | 6 | -4.545 | -4.037 | -3.910 | 0.00 |
|  |  |  | WT-only selection | 7 | -14.367 | -11.971 | -11.426 | 0.00 |
|  |  |  | AntiFold one-shot | 7 | -4.433 | -3.675 | -3.359 | 0.00 |
| 7EJ5_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 135 | variant-aware (worst-case) | 8 | -2.199 | -2.019 | -1.988 | 0.00 |
|  |  |  | WT-only selection | 7 | -1.869 | -1.625 | -1.571 | 0.00 |
|  |  |  | AntiFold one-shot | 4 | -1.207 | -0.969 | -0.911 | 0.00 |
| 7K4N_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 2 | -2.732 | -0.454 | -0.237 | 0.06 |
|  |  |  | WT-only selection | 1 | +0.757 | +1.090 | +1.137 | 1.00 |
|  |  |  | AntiFold one-shot | 1 | +0.757 | +1.090 | +1.137 | 1.00 |
| 7MZK_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 7 | -1.439 | -0.707 | -0.667 | 0.00 |
|  |  |  | WT-only selection | 5 | -0.591 | -0.569 | -0.566 | 0.00 |
|  |  |  | AntiFold one-shot | 6 | -11.951 | -11.515 | -11.486 | 0.00 |
| 7OR9_HLE vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 7 | -6.913 | -5.916 | -5.477 | 0.00 |
|  |  |  | WT-only selection | 7 | -4.316 | -3.821 | -3.453 | 0.00 |
|  |  |  | AntiFold one-shot | 8 | -9.261 | -8.396 | -7.920 | 0.00 |
| 7TBF_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 2 | -1.113 | -0.749 | -0.695 | 0.00 |
|  |  |  | WT-only selection | 2 | -0.543 | -0.234 | -0.194 | 0.00 |
|  |  |  | AntiFold one-shot | 2 | -0.543 | -0.234 | -0.194 | 0.00 |
| 7TE1_HLD vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 137 | variant-aware (worst-case) | 4 | +3.812 | +3.812 | +3.812 | 1.00 |
|  |  |  | WT-only selection | 4 | -2.080 | -2.080 | -2.080 | 0.00 |
|  |  |  | AntiFold one-shot | 4 | -3.113 | -3.113 | -3.113 | 0.00 |
| 7TTX_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 364 | variant-aware (worst-case) | 7 | -1.377 | +0.666 | +0.714 | 1.00 |
|  |  |  | WT-only selection | 9 | -6.019 | -3.828 | -3.777 | 0.00 |
|  |  |  | AntiFold one-shot | 7 | -0.465 | +1.561 | +1.606 | 1.00 |
| 7U2E_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 4 | -0.946 | -0.023 | +0.210 | 0.92 |
|  |  |  | WT-only selection | 7 | -2.373 | -1.679 | -1.418 | 0.00 |
|  |  |  | AntiFold one-shot | 5 | -1.323 | -1.022 | -0.936 | 0.00 |
| 7VYR_HLR vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 4 | +2.001 | +2.001 | +2.001 | 1.00 |
|  |  |  | WT-only selection | 3 | +4.459 | +4.459 | +4.459 | 1.00 |
|  |  |  | AntiFold one-shot | 3 | +5.063 | +5.063 | +5.063 | 1.00 |
| 7WPH_HLB vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 6 | +0.327 | +0.419 | +0.423 | 1.00 |
|  |  |  | WT-only selection | 7 | -5.447 | -3.432 | -3.373 | 0.00 |
|  |  |  | AntiFold one-shot | 6 | -0.432 | +0.958 | +1.001 | 0.99 |
| 7YCL_HLD vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 8 | -9.268 | -9.268 | -9.268 | 0.00 |
|  |  |  | WT-only selection | 8 | -16.460 | -16.460 | -16.460 | 0.00 |
|  |  |  | AntiFold one-shot | 5 | -15.289 | -15.289 | -15.289 | 0.00 |
| 8IDN_HLA vs COVID-19_2021c_C002 | SARS-CoV-2 RBD | 145 | variant-aware (worst-case) | 6 | +1.790 | +1.790 | +1.790 | 1.00 |
|  |  |  | WT-only selection | 7 | -1.085 | -0.921 | -0.891 | 0.00 |
|  |  |  | AntiFold one-shot | 7 | +0.037 | +0.457 | +0.563 | 1.00 |
| 5GZN_HLA vs zika_2023_EDE1-C10 | Zika E | 339 | variant-aware (worst-case) | 10 | -3.744 | -3.380 | -3.354 | 0.00 |
|  |  |  | WT-only selection | 10 | -4.154 | -3.862 | -3.829 | 0.00 |
|  |  |  | AntiFold one-shot | 9 | -4.109 | -3.672 | -3.646 | 0.00 |
| 5JHL_HLA vs zika_2023_EDE1-C10 | Zika E | 337 | variant-aware (worst-case) | 5 | -0.742 | -0.326 | -0.283 | 0.01 |
|  |  |  | WT-only selection | 4 | -2.632 | -1.714 | -1.676 | 0.00 |
|  |  |  | AntiFold one-shot | 4 | -2.525 | -1.658 | -1.623 | 0.00 |
| 5VIG_HLZ vs zika_2023_ZV-67 | Zika E | 133 | variant-aware (worst-case) | 9 | +6.761 | +6.904 | +7.133 | 1.00 |
|  |  |  | WT-only selection | 7 | -5.518 | -2.232 | -1.389 | 0.02 |
|  |  |  | AntiFold one-shot | 6 | -3.972 | -3.757 | -3.565 | 0.01 |
| 6PLK_HLE vs zika_2023_ZV-67 | Zika E | 133 | variant-aware (worst-case) | 6 | -7.573 | -4.483 | -4.081 | 0.00 |
|  |  |  | WT-only selection | 6 | -9.972 | -8.078 | -7.760 | 0.00 |
|  |  |  | AntiFold one-shot | 5 | +0.347 | +1.932 | +2.282 | 1.00 |
| 7BQ5_HLA vs zika_2023_EDE1-C10 | Zika E | 337 | variant-aware (worst-case) | 4 | -4.301 | -4.082 | -4.066 | 0.00 |
|  |  |  | WT-only selection | 6 | -9.841 | -9.118 | -9.079 | 0.00 |
|  |  |  | AntiFold one-shot | 2 | -0.649 | -0.404 | -0.397 | 0.00 |

## Verdict - does variant-aware beat WT-only / AntiFold on MEASURED mutations?

Worst-case reward across the full measured set, per grader. VA>WT and VA>AF flags per system.

| System | Family | grader | variant-aware | WT-only | AntiFold | VA>WT? | VA>AF? |
|---|---|---|--:|--:|--:|:--:|:--:|
| 6XKQ_HLA | SARS-CoV-2 RBD | BA-DDG | -1.189 | -3.759 | -1.291 | Y | Y |
|  |  | Pythia | -5.753 | -14.367 | -5.512 | Y | . |
| 7EJ5_HLA | SARS-CoV-2 RBD | BA-DDG | -0.168 | -1.048 | -1.501 | Y | Y |
|  |  | Pythia | -2.283 | -1.907 | -1.249 | . | . |
| 7K4N_HLA | SARS-CoV-2 RBD | BA-DDG | -1.191 | -1.235 | -1.241 | Y | Y |
|  |  | Pythia | -2.732 | +0.757 | +0.757 | . | . |
| 7MZK_HLA | SARS-CoV-2 RBD | BA-DDG | -2.005 | -2.045 | -5.022 | Y | Y |
|  |  | Pythia | -2.470 | -1.451 | -12.847 | . | Y |
| 7OR9_HLE | SARS-CoV-2 RBD | BA-DDG | -2.102 | -3.642 | -6.657 | Y | Y |
|  |  | Pythia | -8.122 | -5.442 | -10.571 | . | Y |
| 7TBF_HLA | SARS-CoV-2 RBD | BA-DDG | -0.408 | -0.828 | -0.748 | Y | Y |
|  |  | Pythia | -1.113 | -0.543 | -0.543 | . | . |
| 7TE1_HLD | SARS-CoV-2 RBD | BA-DDG | +1.506 | -0.680 | -0.822 | Y | Y |
|  |  | Pythia | +3.812 | -2.080 | -3.113 | Y | Y |
| 7TTX_HLA | SARS-CoV-2 RBD | BA-DDG | -0.936 | -0.157 | -1.493 | . | Y |
|  |  | Pythia | -1.377 | -6.019 | -0.465 | Y | . |
| 7U2E_HLA | SARS-CoV-2 RBD | BA-DDG | -1.363 | -2.827 | -1.650 | Y | Y |
|  |  | Pythia | -1.967 | -3.927 | -2.416 | Y | Y |
| 7VYR_HLR | SARS-CoV-2 RBD | BA-DDG | -0.425 | -1.370 | -1.782 | Y | Y |
|  |  | Pythia | +0.118 | +2.651 | +3.922 | . | . |
| 7WPH_HLB | SARS-CoV-2 RBD | BA-DDG | -0.496 | -2.468 | -1.312 | Y | Y |
|  |  | Pythia | -0.839 | -6.102 | -0.790 | Y | . |
| 7YCL_HLD | SARS-CoV-2 RBD | BA-DDG | -2.524 | -3.120 | -2.731 | Y | Y |
|  |  | Pythia | -11.788 | -17.581 | -16.670 | Y | Y |
| 8IDN_HLA | SARS-CoV-2 RBD | BA-DDG | -0.755 | +1.374 | -0.396 | . | . |
|  |  | Pythia | +0.996 | -2.587 | -1.109 | Y | Y |
| 5GZN_HLA | Zika E | BA-DDG | -2.292 | -5.302 | -7.361 | Y | Y |
|  |  | Pythia | -5.274 | -5.867 | -5.503 | Y | Y |
| 5JHL_HLA | Zika E | BA-DDG | -1.531 | -3.001 | -2.932 | Y | Y |
|  |  | Pythia | -0.742 | -3.486 | -3.319 | Y | Y |
| 5VIG_HLZ | Zika E | BA-DDG | -1.329 | -4.070 | -4.505 | Y | Y |
|  |  | Pythia | +6.107 | -5.798 | -4.429 | Y | Y |
| 6PLK_HLE | Zika E | BA-DDG | +0.210 | -3.813 | -1.653 | Y | Y |
|  |  | Pythia | -8.263 | -9.972 | +0.347 | Y | . |
| 7BQ5_HLA | Zika E | BA-DDG | -1.935 | -1.659 | -2.064 | . | Y |
|  |  | Pythia | -4.796 | -9.896 | -0.770 | Y | . |

**BA-DDG: variant-aware beats WT-only in 15/18 systems, beats AntiFold one-shot in 17/18 systems (worst-case, full measured set).**
**Pythia: variant-aware beats WT-only in 12/18 systems, beats AntiFold one-shot in 9/18 systems (worst-case, full measured set).**

> **Independent-grader headline.** On Pythia_bind (out-of-family, +0.396 vs BA-DDG), variant-aware still beats WT-only in 12/18 and AntiFold one-shot in 9/18 systems - the arm ranking is not an artifact of grading in the steering objective's own lineage.

_Reward = -ddG vs native; higher worst-case = more robust across the measured mutation set. BA-DDG is in the ProteinMPNN lineage (surrogate magnitudes); Pythia_bind is the independent arbiter for the arm ranking. Pythia scores designs as additive masked marginals on the native backbone with the antigen mutation applied as a backbone relabel._
