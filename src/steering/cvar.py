"""CVaR aggregation for the variant-aware objective.

This is the reference implementation of the aggregation used during steering. The
in-loop copy lives in ``variant_aware.py`` (function ``_agg_reward``); the logic here is
the same, lifted into a dependency-free pure function so it can be unit-tested.

SIGN CONVENTION -- the part that is easy to get wrong
-----------------------------------------------------
Two orientations are in play and they are negatives of each other.

*Margin space* (how results are reported). For antigen variant j:

    delta_j(x) = score(x, antigen_j) - score(native, antigen_j)

with scores oriented so that LOWER is better. So ``delta_j(x) < 0`` means the designed
antibody is predicted to bind variant j better than the native antibody does, and
LARGER delta is WORSE. The bad tail is therefore the UPPER tail:

    m_VA(x) = CVaR_alpha({delta_j(x)})        # mean of the worst (largest) alpha-fraction

*Reward space* (what the decoder maximises). ``reward_j = -delta_j``, so HIGHER reward is
better and the bad tail is the LOWER tail:

    R(x) = CVaR_alpha over the lowest alpha-fraction of rewards = -m_VA(x)

``aggregate()`` below operates in REWARD space, matching the steering loop. Use
``cvar_margin()`` for the margin-space quantity reported in the paper.

alpha = 1.0 recovers the mean; alpha -> 0 recovers the strict worst case.
The tail size is ``k = max(1, ceil(alpha * n))``.
"""
from __future__ import annotations

import math
from typing import Iterable, Sequence


def tail_size(n: int, alpha: float) -> int:
    """Number of variants in the bad tail. Always at least one."""
    if n <= 0:
        raise ValueError("need at least one variant")
    return max(1, math.ceil(alpha * n))


def cvar_reward(rewards: Sequence[float], alpha: float = 0.2) -> float:
    """Mean of the LOWEST ``k`` rewards (the worst tail in reward space)."""
    v = sorted(float(r) for r in rewards)
    k = tail_size(len(v), alpha)
    return sum(v[:k]) / k


def cvar_margin(margins: Sequence[float], alpha: float = 0.2) -> float:
    """m_VA: mean of the HIGHEST ``k`` margins (the worst tail in margin space)."""
    v = sorted(float(m) for m in margins)
    k = tail_size(len(v), alpha)
    return sum(v[len(v) - k:]) / k


def aggregate(rewards: Sequence[float], mode: str = "cvar", alpha: float = 0.2) -> float:
    """Reward-space aggregation. ``mode`` in {mean, worst, cvar}.

    Mirrors ``variant_aware.py::_agg_reward`` for these three modes. The ``fused``
    (escape-likelihood-weighted) mode is only implemented in the steering loop.
    """
    v = [float(r) for r in rewards]
    if not v:
        raise ValueError("need at least one variant")
    if mode == "mean":
        return sum(v) / len(v)
    if mode == "worst":
        return min(v)
    if mode == "cvar":
        return cvar_reward(v, alpha)
    raise ValueError(f"unknown aggregation mode {mode!r} (mean|worst|cvar)")


def objective(margins: Iterable[float], alpha: float = 0.2) -> float:
    """R(x) = -m_VA(x): the quantity the steering loop maximises, from margins."""
    m = list(margins)
    return -cvar_margin(m, alpha)
