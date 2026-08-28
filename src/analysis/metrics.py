"""Primary held-out endpoint: worst held-out margin and the beat-native decision.

Margins are design-minus-native and oriented so LOWER IS BETTER, so:

  * a design's score on the held-out panel is its WORST (maximum) margin;
  * it beats native only when that worst margin is < 0.

Pooling across replicates is monotone in pool size -- always compare rep-matched arms.
"""
from __future__ import annotations
from typing import Iterable, Mapping, Sequence


def worst_heldout_margin(margins: Sequence[float]) -> float:
    """Worst = largest margin over the held-out variants."""
    if not len(margins):
        raise ValueError("no held-out margins")
    return max(float(m) for m in margins)


def beats_native(margins: Sequence[float]) -> bool:
    return worst_heldout_margin(margins) < 0.0


def best_of_pool(per_candidate_margins: Iterable[Sequence[float]]) -> float:
    """Ship the candidate with the least-bad worst held-out margin."""
    worsts = [worst_heldout_margin(m) for m in per_candidate_margins]
    if not worsts:
        raise ValueError("empty candidate pool")
    return min(worsts)


def beat_native_rate(per_target_worst: Mapping[str, float]) -> float:
    """Fraction of targets whose shipped design beats native."""
    if not per_target_worst:
        return float("nan")
    return sum(1 for v in per_target_worst.values() if v < 0) / len(per_target_worst)


def paired_targets(a: Mapping[str, float], b: Mapping[str, float]) -> list:
    """Targets present in BOTH conditions -- all paired stats intersect first."""
    return sorted(set(a) & set(b))
