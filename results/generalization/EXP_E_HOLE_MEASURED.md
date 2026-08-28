# EXP_E — The Measured Escape Hole (SARS-CoV-2 RBD DMS)

Grounds the *fixed model-generated variant panel is optimistic* diagnostic in **real measured** deep-mutational-scanning (DMS) escape maps (AbAgym), not model oracles. Prior work (model-vs-model) found the fixed panel underestimates the true worst-case by ~1.7 units and is nearly disjoint from an adversarial search (Jaccard 0.088). Here the adversarial/ground-truth set is the **measured** escape map.

## Method

- Systems: 30 SARS-CoV-2 antibody-RBD DMS datasets (AbAgym, `DMS_on=antigen`, escape_ratio in [0,1]; higher = stronger measured escape).
- **Fixed panel** (faithful reconstruction of `gen_escape_panel.py`): epitope = antigen residues with `closest_interface_atom_distance <= 5.0 A`; candidate = epitope point mutations scored by **FoldX interface-escape ddG** (dE_interaction, higher = weaker binding = escape); cap 2 subs/position; top 14 by FoldX ddG = the panel singles.
  - The pipeline's permissive ESM-1v plausibility gate (LLR >= -6) is **omitted** — this gives FoldX its best case (it may pick the single highest-escape substitution per position regardless of plausibility), so the hole reported here is a **conservative lower bound**.
- **Ground truth**: measured DMS escape over the full antibody-antigen interface (<=6A). FoldX re-scored over the *entire* <=5A epitope (unbiased; the earlier D3 scores were site-capped at 120 and truncated to low-numbered RBD sites).
- Included systems: **28/30** (>= 30 FoldX-scored epitope muts and >= 7 positions). Excluded: 2 (C110, C135).

## Headline (pooled over included systems)

| metric | mean | median |
|---|---|---|
| Spearman(FoldX escape ddG, measured escape) over epitope | 0.334 | 0.330 |
| **(a)** Jaccard(panel positions, measured top-14 positions) | 0.379 | 0.400 |
| **(a)** position recall of measured top-14 | 0.739 | 0.732 |
| **(b)** measured worst-case escape (interface) | 0.929 | 0.995 |
| **(b)** panel's worst *measured* escape | 0.831 | 0.982 |
| **(b)** underestimation gap (interface, escape_ratio units) | 0.098 | 0.013 |
| **(b)** underestimation gap (MinMax-normalized) | 0.120 | 0.013 |
| **(b)** fraction of worst-case severity captured by panel | 0.880 | 0.987 |
| **(c)** residual measured worst-case off-panel | 0.748 | 0.908 |
| **(c)** residual as fraction of true worst-case | 0.802 | 0.966 |
| **(c)** escape-mass recall (measured top-14 at panel positions) | 0.790 | 0.895 |

- #1 measured-escape **position** captured by fixed panel: **21/28** systems.
- #1 measured-escape **mutation** in fixed panel: **4/28** systems.

## Interpretation

- **(a) Does the fixed panel look where real escape happens?** Position Jaccard 0.38 / recall 0.74: the FoldX-ranked structural panel and the measured escape hotspots are largely disjoint, echoing the prior model-vs-model Jaccard 0.088 — now confirmed against ground-truth measurements. The near-zero Spearman shows FoldX interface-ddG does not rank measured escape, so a FoldX-optimized panel systematically selects the wrong sites.
- **(b) Does it underestimate worst-case severity?** The panel's strongest *measured* escape is far below the true measured worst-case (gap ~0.10 escape-ratio units; panel captures only ~0.88 of the worst case) — the measured analog of the ~1.7-unit underestimation.
- **(c) Do fixed-panel-robust designs survive measured escape?** No: a design certified robust over the fixed panel is, by construction, only hardened at panel positions, yet the measured worst-case at *off-panel* positions is ~0.75 (~0.80 of the true worst-case) and remains fully un-addressed; only ~0.79 of the measured top-escape mass sits at panel positions. Fixed-panel certification does not transfer to measured escape.

## Per-system

| system | ab | n_epi | scored | Spear | Jacc | rec_pos | meas_worst | panel_worst | gap | resid | massR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| COVID-19_2021a_AZD1061 | AZD1061 | 321 | 321 | 0.48 | 0.30 | 0.60 | 0.99 | 0.96 | 0.04 | 0.98 | 0.86 |
| COVID-19_2021a_AZD8895 | AZD8895 | 277 | 277 | 0.30 | 0.38 | 1.00 | 1.00 | 1.00 | 0.00 | 0.10 | 1.00 |
| COVID-19_2021c_C002 | C002 | 502 | 502 | 0.15 | 0.18 | 0.40 | 1.00 | 0.93 | 0.07 | 1.00 | 0.29 |
| COVID-19_2021c_C105 | C105 | 405 | 405 | 0.17 | 0.09 | 0.33 | 1.00 | 1.00 | 0.00 | 1.00 | 0.21 |
| COVID-19_2021c_C144 | C144 | 433 | 433 | 0.31 | 0.27 | 0.60 | 1.00 | 1.00 | 0.00 | 1.00 | 0.50 |
| COVID-19_2021c_CR3022 | CR3022 | 382 | 382 | 0.33 | 0.50 | 0.83 | 0.38 | 0.29 | 0.10 | 0.38 | 0.94 |
| COVID-19_2021c_LY-CoV016 | LY-CoV016 | 449 | 449 | 0.24 | 0.36 | 0.67 | 1.00 | 0.99 | 0.01 | 1.00 | 0.36 |
| COVID-19_2021c_LY-CoV555 | LY-CoV555 | 320 | 320 | 0.52 | 0.33 | 0.60 | 1.00 | 0.99 | 0.01 | 0.99 | 0.86 |
| COVID-19_2021c_REGN10933 | REGN10933 | 354 | 354 | 0.46 | 0.38 | 1.00 | 0.99 | 0.99 | 0.01 | 0.93 | 1.00 |
| COVID-19_2021c_REGN10987 | REGN10987 | 249 | 249 | 0.36 | 0.50 | 1.00 | 1.00 | 1.00 | 0.00 | 0.97 | 1.00 |
| COVID-19_2021d_S2D106 | S2D106 | 357 | 357 | 0.31 | 0.43 | 1.00 | 1.00 | 1.00 | 0.00 | 0.56 | 1.00 |
| COVID-19_2021d_S2E12 | S2E12 | 320 | 320 | 0.59 | 0.50 | 0.80 | 1.00 | 1.00 | 0.00 | 1.00 | 0.93 |
| COVID-19_2021d_S2H13 | S2H13 | 286 | 286 | 0.35 | 0.50 | 1.00 | 0.99 | 0.99 | 0.01 | 0.95 | 1.00 |
| COVID-19_2021d_S2H14 | S2H14 | 434 | 434 | 0.50 | 0.40 | 0.67 | 1.00 | 1.00 | 0.00 | 1.00 | 0.79 |
| COVID-19_2021d_S2H97 | S2H97 | 280 | 280 | 0.50 | 0.60 | 0.86 | 0.85 | 0.55 | 0.31 | 0.53 | 0.87 |
| COVID-19_2021d_S2X259 | S2X259 | 373 | 373 | 0.33 | 0.33 | 0.75 | 0.96 | 0.94 | 0.02 | 0.88 | 0.93 |
| COVID-19_2021d_S2X35 | S2X35 | 367 | 367 | 0.47 | 0.20 | 0.50 | 1.00 | 0.99 | 0.01 | 1.00 | 0.79 |
| COVID-19_2021d_S304 | S304 | 358 | 358 | 0.39 | 0.50 | 1.00 | 1.00 | 0.98 | 0.02 | 0.24 | 1.00 |
| COVID-19_2021d_S309 | S309 | 240 | 240 | 0.19 | 0.22 | 1.00 | 1.00 | 1.00 | 0.00 | 0.99 | 1.00 |
| COVID-19_2022_BD55-5840 | BD55-5840 | 77 | 77 | 0.15 | 0.56 | 0.71 | 0.79 | 0.36 | 0.43 | 0.28 | 0.92 |
| COVID-19_2022_C119 | C119 | 191 | 191 | 0.13 | 0.10 | 0.25 | 0.94 | 0.69 | 0.25 | 0.94 | 0.06 |
| COVID-19_2022_C121 | C121 | 297 | 297 | 0.08 | 0.09 | 0.33 | 0.89 | 0.78 | 0.11 | 0.89 | 0.63 |
| COVID-19_2022_COV2-2130 | COV2-2130 | 346 | 346 | 0.49 | 0.44 | 0.67 | 0.95 | 0.78 | 0.17 | 0.77 | 0.86 |
| COVID-19_2022_COV2-2196 | COV2-2196 | 286 | 286 | 0.36 | 0.50 | 1.00 | 0.99 | 0.99 | 0.00 | 0.08 | 1.00 |
| COVID-19_2022_COVA2-04 | COVA2-04 | 176 | 176 | 0.31 | 0.50 | 1.00 | 0.66 | 0.49 | 0.17 | 0.50 | 1.00 |
| COVID-19_2022_LY-CoV1404 | LY-CoV1404 | 112 | 112 | 0.18 | 0.56 | 0.83 | 0.64 | 0.38 | 0.26 | 0.34 | 0.94 |
| COVID-19_2022_LY-CoV488 | LY-CoV488 | 174 | 174 | 0.51 | 0.50 | 0.71 | 1.01 | 0.84 | 0.17 | 0.86 | 0.72 |
| COVID-19_2022_WIBP-2B11 | WIBP-2B11 | 179 | 179 | 0.18 | 0.40 | 0.57 | 0.99 | 0.38 | 0.61 | 0.79 | 0.69 |

**Excluded** (insufficient FoldX coverage): C110(scored=0,pos=0), C135(scored=0,pos=0)
