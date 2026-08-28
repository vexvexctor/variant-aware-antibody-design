# Novel out-of-distribution antibody-design eval

Targets reported: 10/10 (position-disjoint held-out escape split; best-of-6-pool worst-case over HELD-OUT eval variants).

Beat-native: design worst-case-over-eval strictly better (tighter) than the native antibody's. H3-DDG native baseline = 0 (win = margin < 0); FoldX = absolute interaction E vs native WT worst (win = margin < 0).

| target | antigen | organism | ag-id%vs-train | n_steer | n_eval | headline H3 | headline FoldX | AntiFold H3 | AntiFold FoldX |
|---|---|---|---|---|---|---|---|---|---|
| 9MA7_HLA | Envelope glycoprotein GP35 | human gammaherpesvirus | 0.0 | 6 | 12 | -0.15 ✅ | -0.05 ✅ | +0.35 ❌ | -0.12 ✅ |
| 9BQW_HLA | Decorin-binding protein A | Borreliella burgdorfer | 0.0 | 6 | 11 | +0.06 ❌ | -0.08 ✅ | +0.53 ❌ | -0.66 ✅ |
| 9JT1_HLA | Large envelope protein | HBV genotype D3 | 2.6 | 6 | 11 | +0.27 ❌ | -0.38 ✅ | +0.55 ❌ | +0.98 ❌ |
| 9OG2_HLA | A28L | Monkeypox virus | 6.4 | 6 | 12 | -0.70 ✅ | -3.03 ✅ | -1.67 ✅ | -1.66 ✅ |
| 9YIO_HLA | PvRipr EGF7-8 | Plasmodium vivax | 15.1 | 6 | 8 | -0.42 ✅ | -3.08 ✅ | -0.24 ✅ | -3.09 ✅ |
| 9LY5_HLA | CD276 antigen | Homo sapiens | 16.3 | 6 | 8 | -0.82 ✅ | -4.73 ✅ | +0.22 ❌ | -3.21 ✅ |
| 8RWB_HLA | UL16-binding protein 6 | Homo sapiens | 17.4 | 6 | 14 | -0.35 ✅ | -3.21 ✅ | -1.03 ✅ | -2.67 ✅ |
| 9KKJ_HLA | Nectin-4 | Homo sapiens | 29.8 | 6 | 14 | -0.26 ✅ | -3.11 ✅ | +1.84 ❌ | -0.33 ✅ |
| 9VL2_HLA | Myeloid cell surface antig | Homo sapiens | 40.0 | 6 | 11 | +0.10 ❌ | -3.59 ✅ | +0.39 ❌ | -0.54 ✅ |
| 9PI9_HLA | Tumor-associated calcium s | Homo sapiens | 48.6 | 6 | 8 | -0.42 ✅ | -0.55 ✅ | -0.12 ✅ | -0.90 ✅ |

## Summary

| arm | oracle | beat-native rate |
|---|---|---|
| headline | H3 | 7/10 (70%) |
| headline | FOLDX | 10/10 (100%) |
| antifold | H3 | 4/10 (40%) |
| antifold | FOLDX | 9/10 (90%) |

### Paired headline-vs-AntiFold (both scored on the oracle)

| oracle | headline-only wins | AntiFold-only wins | same (both win/both lose) |
|---|---|---|---|
| H3 | 3 | 0 | 7 |
| FOLDX | 1 | 0 | 9 |
