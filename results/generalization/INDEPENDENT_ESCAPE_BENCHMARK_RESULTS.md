# Independent-escape benchmark — beat-native on an independently-constructed escape set

Headline cell **svdd_af_baddg** designs (already generated) re-graded by the out-of-family
FoldX physics arbiter on escape variants built WITHOUT the ProteinMPNN/ESM-1v lineage of
the optimization (omega) panel. Metric = best-of-pool WORST-CASE interaction energy across
the panel's variants minus the native antibody's worst-case; margin<0 => design pool's
least-favorable variant binds tighter than native's least-favorable variant (more negative
energy = tighter). This is the SAME metric used for the omega-panel beat-native in the paper.

## Per-target

### Independent panel: NATURAL VOC (documented natural/clinical variants; model-free)
- **7MZK_HLA** (46 indep variants): nat_worst=  -2.72 design_bp_worst=  -1.12 margin=  +1.60 [loss] (scored 46v/8d)
    omega panel (24 eval variants):    nat_worst=  -5.72 design_bp_worst=  -6.78 margin=  -1.06 [WIN] (scored 24v/8d)
- **7OR9_HLE** (46 indep variants): nat_worst=  -0.30 design_bp_worst=  -2.39 margin=  -2.09 [WIN] (scored 46v/8d)
    omega panel (24 eval variants):    nat_worst=  -5.10 design_bp_worst=  -7.45 margin=  -2.35 [WIN] (scored 24v/8d)
- **7TBF_HLA** (46 indep variants): nat_worst=  -1.26 design_bp_worst=  -2.23 margin=  -0.97 [WIN] (scored 46v/8d)
    omega panel (24 eval variants):    nat_worst=  -0.74 design_bp_worst=  -0.77 margin=  -0.02 [WIN] (scored 24v/8d)
- **7VYR_HLR** (46 indep variants): nat_worst=  -7.59 design_bp_worst=  -7.59 margin=  +0.01 [loss] (scored 46v/8d)
    omega panel (24 eval variants):    nat_worst=  -6.95 design_bp_worst=  -8.84 margin=  -1.89 [WIN] (scored 24v/8d)
- **7WPH_HLB** (46 indep variants): nat_worst=  -3.25 design_bp_worst=  -5.61 margin=  -2.37 [WIN] (scored 46v/8d)
    omega panel (24 eval variants):    nat_worst=  15.61 design_bp_worst=   8.09 margin=  -7.52 [WIN] (scored 24v/8d)
- **7YCL_HLD** (46 indep variants): nat_worst=  -7.14 design_bp_worst=  -9.09 margin=  -1.95 [WIN] (scored 46v/8d)
    omega panel (24 eval variants):    nat_worst=  -0.23 design_bp_worst=  -5.22 margin=  -5.00 [WIN] (scored 24v/8d)
- **8IDN_HLA** (46 indep variants): nat_worst=  -0.66 design_bp_worst=  -1.96 margin=  -1.30 [WIN] (scored 46v/8d)
    omega panel (24 eval variants):    nat_worst=   7.40 design_bp_worst=  -2.88 margin= -10.27 [WIN] (scored 24v/8d)

### Independent panel: ESM-2 masked-marginal epitope (out-of-lineage model)
- **5DUM_HLA** (36 indep variants): nat_worst=  -3.48 design_bp_worst=  -7.85 margin=  -4.37 [WIN] (scored 36v/8d)
    omega panel (24 eval variants):    nat_worst=  -2.00 design_bp_worst=  -6.76 margin=  -4.76 [WIN] (scored 24v/6d)
- **5GJS_HLB** (36 indep variants): nat_worst=  -6.31 design_bp_worst=  -5.43 margin=  +0.88 [loss] (scored 36v/8d)
    omega panel (24 eval variants):    nat_worst=  -3.29 design_bp_worst=  -1.31 margin=  +1.98 [loss] (scored 24v/8d)
- **5UGY_HLA** (36 indep variants): nat_worst=  10.67 design_bp_worst=   7.15 margin=  -3.52 [WIN] (scored 36v/8d)
    omega panel (24 eval variants):    nat_worst=  18.62 design_bp_worst=  12.56 margin=  -6.06 [WIN] (scored 24v/8d)
- **6MLM_KLH** (36 indep variants): nat_worst=   4.08 design_bp_worst=   0.98 margin=  -3.10 [WIN] (scored 36v/8d)
    omega panel (24 eval variants):    nat_worst=  30.21 design_bp_worst=  24.78 margin=  -5.43 [WIN] (scored 24v/8d)
- **7MZK_HLA** (28 indep variants): nat_worst=  -6.99 design_bp_worst=  -6.94 margin=  +0.05 [loss] (scored 28v/8d)
    omega panel (24 eval variants):    nat_worst=  -5.72 design_bp_worst=  -6.78 margin=  -1.06 [WIN] (scored 24v/8d)
- **7OR9_HLE** (36 indep variants): nat_worst=  -2.61 design_bp_worst=  -4.73 margin=  -2.11 [WIN] (scored 36v/8d)
    omega panel (24 eval variants):    nat_worst=  -5.10 design_bp_worst=  -7.45 margin=  -2.35 [WIN] (scored 24v/8d)
- **7TBF_HLA** (34 indep variants): nat_worst=  -0.77 design_bp_worst=  -0.86 margin=  -0.08 [WIN] (scored 34v/8d)
    omega panel (24 eval variants):    nat_worst=  -0.74 design_bp_worst=  -0.77 margin=  -0.02 [WIN] (scored 24v/8d)
- **7TTX_HLA** (36 indep variants): nat_worst=  -7.91 design_bp_worst=  -9.49 margin=  -1.58 [WIN] (scored 36v/8d)
    omega panel (24 eval variants):    nat_worst=  -8.49 design_bp_worst=  -9.99 margin=  -1.50 [WIN] (scored 24v/8d)
- **7VYR_HLR** (36 indep variants): nat_worst=  -5.66 design_bp_worst=  -6.83 margin=  -1.18 [WIN] (scored 36v/8d)
    omega panel (24 eval variants):    nat_worst=  -6.95 design_bp_worst=  -8.84 margin=  -1.89 [WIN] (scored 24v/8d)
- **7WPH_HLB** (32 indep variants): nat_worst=  -0.94 design_bp_worst=  -2.78 margin=  -1.84 [WIN] (scored 32v/8d)
    omega panel (24 eval variants):    nat_worst=  15.61 design_bp_worst=   8.09 margin=  -7.52 [WIN] (scored 24v/8d)
- **7X4I_HA** (36 indep variants): nat_worst=  -2.23 design_bp_worst=  -1.57 margin=  +0.66 [loss] (scored 36v/8d)
    omega panel (24 eval variants):    nat_worst=  -2.83 design_bp_worst=  -4.12 margin=  -1.29 [WIN] (scored 24v/8d)
- **7YCL_HLD** (36 indep variants): nat_worst=  -8.04 design_bp_worst=  -9.67 margin=  -1.63 [WIN] (scored 36v/8d)
    omega panel (24 eval variants):    nat_worst=  -0.23 design_bp_worst=  -5.22 margin=  -5.00 [WIN] (scored 24v/8d)
- **8IDN_HLA** (36 indep variants): nat_worst=   0.84 design_bp_worst=  -2.71 margin=  -3.55 [WIN] (scored 36v/8d)
    omega panel (24 eval variants):    nat_worst=   7.40 design_bp_worst=  -2.88 margin= -10.27 [WIN] (scored 24v/8d)

## Summary (beat-native rate)
- **natural independent panel**: beat-native 5/7 = 71%  (mean margin -1.01)
    same 7 targets, ORIGINAL omega panel: beat-native 7/7 = 100% (mean margin -4.02); independent on those: 5/7 = 71%
- **esm2 independent panel**: beat-native 10/13 = 77%  (mean margin -1.65)
    same 13 targets, ORIGINAL omega panel: beat-native 12/13 = 92% (mean margin -3.47); independent on those: 10/13 = 77%

