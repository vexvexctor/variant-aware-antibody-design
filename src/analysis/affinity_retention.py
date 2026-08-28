#!/usr/bin/env python3
"""Definitive model for the escape-penalty claim.

Three earlier cuts disagreed on the held-out evaluator (H3), for understandable reasons:
  - per-design Mann-Whitney inside strata : unpaired (weak) AND pseudo-replicated (bad p)
  - paired-by-target inside strata        : properly paired, but residual within-stratum
                                            WT-binding imbalance still favours variant-aware
  - pooled cubic ANCOVA                   : confound removed, but pairing thrown away

This does both at once:

    penalty ~ target fixed effects + poly(s_wt, 3) + beta * 1[variant-aware]

Target fixed effects absorb every between-target difference (backbone, epitope, panel
difficulty) = the pairing. The cubic in s(WT) absorbs the affinity confound continuously = the
matching. beta is what remains: the difference in how much binding is lost under antigen change,
between two designs of the SAME target that bind WT EQUALLY WELL.

CIs bootstrap antigen clusters. A placebo column refits with the arm labels shuffled WITHIN
target, which should give beta ~ 0 if the machinery is sound.
"""
import numpy as np, pandas as pd
from pathlib import Path
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

OUT = Path(f'{VAAD_ROOT}/scratch/paper_checks')
RD = Path(f'{VAAD_ROOT}/results/report_data')
RNG = np.random.default_rng(0)
VA, UN = 'f5_worst_wt1', 'f5_unsteer'
DEG = 3


def design_matrix(g):
    tgt = pd.get_dummies(g.target, drop_first=True).to_numpy(float)
    poly = np.column_stack([g.z.values ** k for k in range(1, DEG + 1)])
    return np.column_stack([np.ones(len(g)), tgt, poly, g.va.values])


def beta_of(g):
    X = design_matrix(g)
    try:
        coef, *_ = np.linalg.lstsq(X, g.penalty.values, rcond=None)
    except np.linalg.LinAlgError:
        return np.nan
    return coef[-1]


def main():
    d = pd.read_csv(OUT / 'escape_penalty_per_design.csv')
    clu = pd.read_csv(RD / 'clustered_holdout_split.csv')[['target', 'antigen_cluster_70']]
    cmap = dict(zip(clu.target, clu.antigen_cluster_70))

    rows = []
    for ev, sub in d.groupby('evaluator'):
        s = sub[sub.arm.isin([VA, UN])].copy()
        s['cluster'] = s.target.map(cmap)
        s = s.dropna(subset=['cluster', 's_wt', 'penalty'])
        # keep only targets where BOTH arms are present, else the target FE eats the contrast
        ok = s.groupby('target').arm.nunique()
        s = s[s.target.isin(ok[ok == 2].index)]
        if len(s) < 60 or s.arm.nunique() < 2:
            print(f'[{ev}] skipped (n={len(s)})')
            continue
        s['z'] = (s.s_wt - s.s_wt.mean()) / (s.s_wt.std() or 1.0)
        s['va'] = (s.arm == VA).astype(float)

        beta = beta_of(s)
        uc = s.cluster.unique()
        boots, placebo = [], []
        for _ in range(4000):
            pick = RNG.choice(uc, len(uc), replace=True)
            g = pd.concat([s[s.cluster == c] for c in pick])
            ok2 = g.groupby('target').va.nunique()
            g = g[g.target.isin(ok2[ok2 == 2].index)]
            if len(g) < 40:
                continue
            b = beta_of(g)
            if np.isfinite(b):
                boots.append(b)
        for _ in range(500):
            g = s.copy()
            g['va'] = g.groupby('target')['va'].transform(lambda v: RNG.permutation(v.values))
            b = beta_of(g)
            if np.isfinite(b):
                placebo.append(b)
        lo, hi = (np.percentile(boots, [2.5, 97.5]) if len(boots) > 100 else (np.nan, np.nan))
        plo, phi = (np.percentile(placebo, [2.5, 97.5]) if len(placebo) > 100 else (np.nan, np.nan))
        rows.append(dict(evaluator=ev, n_designs=len(s), n_targets=s.target.nunique(),
                         n_clusters=len(uc), beta=beta, ci_lo=lo, ci_hi=hi,
                         excludes_zero=bool(np.isfinite(lo) and ((lo < 0) == (hi < 0))),
                         placebo_lo=plo, placebo_hi=phi))

    r = pd.DataFrame(rows)
    r.to_csv(OUT / 'escape_penalty_final.csv', index=False)
    pd.set_option('display.width', 220)
    print('===== ESCAPE PENALTY: target fixed effects + continuous WT-binding adjustment =====')
    print('beta = extra binding lost under antigen change by variant-aware vs unsteered,')
    print('       for designs of the same target that bind WT equally well.')
    print('       NEGATIVE = variant-aware is genuinely more retentive.\n')
    print(r.round(4).to_string(index=False))


if __name__ == '__main__':
    main()
