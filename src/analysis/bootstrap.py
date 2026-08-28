"""Cluster bootstrap.

Resampling unit is the ANTIGEN CLUSTER, never the PDB entry: repeated antigens
(lysozyme, SARS-CoV-2 RBD, ...) would otherwise inflate significance. Cluster labels
come from data/manifests/antigen_clusters.csv (`benchmark_cluster`).

Seeded by default so published intervals are reproducible.
"""
from __future__ import annotations
import random
from typing import Callable, Mapping, Sequence, Tuple

DEFAULT_SEED = 0
DEFAULT_DRAWS = 10000


def cluster_bootstrap(values: Sequence[float], clusters: Sequence[str],
                      statistic: Callable[[Sequence[float]], float] = None,
                      draws: int = DEFAULT_DRAWS, seed: int = DEFAULT_SEED,
                      ci: float = 0.95) -> Tuple[float, float, float]:
    """-> (point estimate, lo, hi). Resamples clusters with replacement."""
    if statistic is None:
        statistic = lambda v: sum(v) / len(v)
    by: dict = {}
    for v, c in zip(values, clusters):
        by.setdefault(c, []).append(float(v))
    keys = list(by)
    if not keys:
        return float("nan"), float("nan"), float("nan")
    rng, out = random.Random(seed), []
    for _ in range(draws):
        sample = [x for k in (rng.choice(keys) for _ in keys) for x in by[k]]
        out.append(statistic(sample))
    out.sort()
    lo_i = int((1 - ci) / 2 * draws)
    hi_i = int((1 - (1 - ci) / 2) * draws) - 1
    return statistic([x for v in by.values() for x in v]), out[lo_i], out[hi_i]


def paired_delta_bootstrap(a: Mapping[str, float], b: Mapping[str, float],
                           cluster_of: Mapping[str, str], draws: int = DEFAULT_DRAWS,
                           seed: int = DEFAULT_SEED) -> Tuple[float, float, float]:
    """Paired b-minus-a beat-native difference in percentage points, cluster-bootstrapped."""
    shared = sorted(set(a) & set(b))
    by: dict = {}
    for t in shared:
        by.setdefault(cluster_of.get(t, t), []).append((a[t], b[t]))
    keys = list(by)
    rate = lambda pairs, i: 100 * sum(1 for p in pairs if p[i] < 0) / len(pairs)
    rng, out = random.Random(seed), []
    for _ in range(draws):
        s = [p for k in (rng.choice(keys) for _ in keys) for p in by[k]]
        out.append(rate(s, 1) - rate(s, 0))
    out.sort()
    allp = [p for v in by.values() for p in v]
    return rate(allp, 1) - rate(allp, 0), out[int(.025 * draws)], out[int(.975 * draws) - 1]
