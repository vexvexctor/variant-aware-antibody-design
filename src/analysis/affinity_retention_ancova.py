#!/usr/bin/env python3
"""Continuous adjustment for the WT-binding confound (the version to put in the paper).

Binning into strata leaves residual imbalance: inside every stratum the variant-aware arm still
binds WT slightly better than the unsteered arm (3-15% of bin width, always the same sign). On
H3 the confound slope is POSITIVE, so that residue flatters variant-aware. Rather than argue the
residue is small, adjust for s(WT) continuously:

    penalty ~ 1 + f(s_wt) + beta * 1[arm == variant-aware]

beta is the escape penalty difference at matched WT binding. f is a cubic polynomial in s_wt so
the adjustment does not assume the confound is linear. CIs come from bootstrapping ANTIGEN
CLUSTERS (designs within a target, and targets within a cluster, are not independent).

Reported for every evaluator. The one that matters for the paper's claim is the evaluator NOT
used during optimization.
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


def fit_beta(swt, pen, is_va):
    """OLS of penalty on [poly(s_wt), arm]; return the arm coefficient."""
    X = np.column_stack([swt ** k for k in range(DEG + 1)] + [is_va])
    try:
        coef, *_ = np.linalg.lstsq(X, pen, rcond=None)
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
        if s.arm.nunique() < 2 or len(s) < 50:
            continue
        s['cluster'] = s.target.map(cmap)
        s = s.dropna(subset=['cluster', 's_wt', 'penalty'])
        # standardise s_wt so the cubic is numerically sane
        mu, sd = s.s_wt.mean(), s.s_wt.std() or 1.0
        s['z'] = (s.s_wt - mu) / sd
        s['va'] = (s.arm == VA).astype(float)

        beta = fit_beta(s.z.values, s.penalty.values, s.va.values)
        raw = s[s.va == 1].penalty.median() - s[s.va == 0].penalty.median()

        uc = s.cluster.unique()
        boots = []
        for _ in range(4000):
            pick = RNG.choice(uc, len(uc), replace=True)
            g = pd.concat([s[s.cluster == c] for c in pick])
            if g.va.nunique() < 2:
                continue
            b = fit_beta(g.z.values, g.penalty.values, g.va.values)
            if np.isfinite(b):
                boots.append(b)
        lo, hi = (np.percentile(boots, [2.5, 97.5]) if len(boots) > 100 else (np.nan, np.nan))
        rows.append(dict(evaluator=ev, n_designs=len(s), n_clusters=len(uc),
                         raw_delta=raw, adjusted_delta=beta, ci_lo=lo, ci_hi=hi,
                         excludes_zero=bool(np.isfinite(lo) and ((lo < 0) == (hi < 0)))))

    r = pd.DataFrame(rows)
    r.to_csv(OUT / 'escape_penalty_ancova.csv', index=False)
    pd.set_option('display.width', 200)
    print('===== ESCAPE PENALTY, ADJUSTED FOR WT BINDING (cubic), antigen-cluster bootstrap =====')
    print('negative = variant-aware loses LESS at matched WT binding\n')
    print(r.round(4).to_string(index=False))
    print('\nraw_delta      = unadjusted median difference (confounded)')
    print('adjusted_delta = arm coefficient after removing the s(WT) trend')


if __name__ == '__main__':
    main()
