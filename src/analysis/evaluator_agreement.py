#!/usr/bin/env python3
"""Evaluator agreement on a shared frozen design pool.

Every evaluator scores the SAME 1,785 (target, seed, selector, design) units, so all
pairs are exactly matched. Emits three tidy matrices:

    spearman.csv      rank agreement of continuous design-minus-native margins
    cohen_kappa.csv   chance-corrected agreement on the binary beat-native decision
    base_rates.csv    each evaluator's positive-call rate on the shared pool

Raw percent agreement is reported in the long-form CSV but is NOT the headline statistic:
two evaluators that both almost always say "not beat-native" agree often by construction.
FoldX and MM-GBSA reach 82.4% raw agreement at kappa = -0.012.

    python -m analysis.evaluator_agreement \
        --pool data/processed/shared_1785_design_pool.csv \
        --outdir results/evaluator_agreement
"""
from __future__ import annotations
import argparse, collections, csv, itertools, math, os

DEFAULT_ORDER = ["BA-DDG", "H3", "Pythia", "FoldX", "MM-GBSA"]


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        rk = [0.0] * len(v); i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                rk[order[k]] = avg
            i = j + 1
        return rk
    rx, ry = rank(xs), rank(ys)
    n = len(xs); mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def load(pool_csv):
    vals = collections.defaultdict(dict)
    for r in csv.DictReader(open(pool_csv)):
        unit = (r["target"], r["seed"], r["selector"], r["design"])
        try:
            vals[r["evaluator"]][unit] = float(r["w"])
        except ValueError:
            pass
    return vals


def pairwise(vals, order):
    rows = []
    for a, b in itertools.combinations(order, 2):
        shared = sorted(set(vals[a]) & set(vals[b]))
        if not shared:
            continue
        xs = [vals[a][u] for u in shared]; ys = [vals[b][u] for u in shared]
        da = [x < 0 for x in xs]; db = [y < 0 for y in ys]
        n = len(shared)
        po = sum(1 for p, q in zip(da, db) if p == q) / n
        pa, pb = sum(da) / n, sum(db) / n
        pe = pa * pb + (1 - pa) * (1 - pb)
        rows.append(dict(evaluator_a=a, evaluator_b=b, n=n,
                         spearman=round(spearman(xs, ys), 3),
                         pct_agree=round(100 * po, 1),
                         cohen_kappa=round((po - pe) / (1 - pe), 3) if pe < 1 else float("nan")))
    return rows


def write_matrix(rows, order, key, path):
    M = {(r["evaluator_a"], r["evaluator_b"]): r[key] for r in rows}
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["evaluator"] + order)
        for a in order:
            w.writerow([a] + ["" if a == b else (M.get((a, b)) if (a, b) in M else M.get((b, a)))
                              for b in order])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", default="data/processed/shared_1785_design_pool.csv")
    ap.add_argument("--outdir", default="results/evaluator_agreement")
    ap.add_argument("--evaluators", nargs="*", default=DEFAULT_ORDER)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    vals = load(a.pool)
    order = [e for e in a.evaluators if e in vals]
    rows = pairwise(vals, order)

    write_matrix(rows, order, "spearman", os.path.join(a.outdir, "spearman.csv"))
    write_matrix(rows, order, "cohen_kappa", os.path.join(a.outdir, "cohen_kappa.csv"))
    write_matrix(rows, order, "pct_agree", os.path.join(a.outdir, "percent_agreement.csv"))
    with open(os.path.join(a.outdir, "pairwise_long.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(os.path.join(a.outdir, "base_rates.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["evaluator", "n_designs", "beat_native_rate_pct"])
        for o in order:
            v = list(vals[o].values())
            w.writerow([o, len(v), round(100 * sum(1 for x in v if x < 0) / len(v), 1)])
    n = len(next(iter(vals.values())))
    print(f"{len(order)} evaluators x {n} shared designs -> {a.outdir}")
    for o in order:
        v = list(vals[o].values())
        print(f"  {o:>9}: beat-native {100*sum(1 for x in v if x<0)/len(v):5.1f}%")


if __name__ == "__main__":
    main()
