#!/usr/bin/env python3
"""Generator-zoo master matrix: NATIVE vs GENERATOR-ALONE vs GENERATOR + OUR OBJECTIVE.

Three quantities, all on the SAME held-out escape panel per target (split == "eval", never
steered against), all with the native affinity-matured antibody as the reference:

  beat-native      best-of-pool worst-case margin vs native < 0      (headline)
  paired delta     (generator+ours) - (generator alone), per target  (isolates the objective
                   from the prior: same base distribution, same pool size, same panel)
  cost             wall / cuda seconds and reward-eval counts per design, from .cost.json

Aggregation is CLUSTER-AWARE: the 200 targets collapse to ~172 mmseqs-70 antigen clusters, and
a bootstrap that resamples TARGETS would overstate confidence wherever an antigen repeats. The
bootstrap here resamples antigen CLUSTERS and takes a random target within each, so the
effective N is the number of independent antigens, not the number of PDB entries.

Robust to partial data — safe to run while grading is still filling in.
"""
from __future__ import annotations

import csv
import json
import os
import random
import statistics as st
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

G = f"{VAAD_ROOT}/scratch/genzoo"
DESIGNS = f"{G}/designs"
H3DIR = f"{G}/h3_out"
OUT_MD = f"{G}/MASTER_MATRIX.md"
OUT_CSV = f"{G}/master_matrix_cells.csv"
REPS = [1, 2, 3]
N_BOOT = 2000

# MATCHED pairs: baseline and steered arm must use the SAME generator in the SAME inference
# mode, so the only difference is the objective. AntiFold's native mode IS independent sampling
# from a fixed logit table, so its frozen pair is matched; the MPNN family and ESM-IF1 are
# autoregressive, so their matched pair is (native AR baseline) vs (live re-conditioned SVDD).
HEADLINE_PAIRS = [
    ("AntiFold",            "base_antifold",         "ours_antifold"),
    ("AbMPNN",              "base_abmpnn_ar",        "ours_abmpnn_live"),
    ("ProteinMPNN",         "base_proteinmpnn_ar",   "ours_proteinmpnn_live"),
    ("ESM-IF1",             "base_esmif_ar",         "ours_esmif_live"),
    ("our masked diffusion", "base_diffusion",       "ours_diffusion"),
]
# ABLATION pairs: the same models with a FROZEN marginal table on both sides. Isolates what
# per-step re-conditioning buys, and is NOT the headline for autoregressive models.
ABLATION_PAIRS = [
    ("AbMPNN (frozen)",      "base_abmpnn_npz",  "ours_abmpnn_npz"),
    ("ProteinMPNN (frozen)", "base_proteinmpnn", "ours_proteinmpnn"),
    ("ESM-IF1 (frozen)",     "base_esmif",       "ours_esmif"),
]
# Generators that design a NEW CDR-H3 CONFORMATION (and, for RFantibody, re-dock the Fv).
# Threading their sequence onto the native backbone discards what they actually produced, so
# their row here is NOT a fair measure of them -- observed 0/23 for RFantibody. Their real
# number comes from grade_own_structure.py, which scores each design on its own complex.
DENOVO = {"base_rfantibody", "base_diffab"}
MODELS = ["antifold", "abmpnn_npz", "proteinmpnn", "esmif", "diffusion"]
DISP_MODEL = {"antifold": "AntiFold", "abmpnn_npz": "AbMPNN", "proteinmpnn": "ProteinMPNN",
              "esmif": "ESM-IF1", "diffusion": "our masked diffusion"}
# every real arm gets an edit-distance-matched random twin: `rand_<arm>` mutates the SAME number
# of CDR positions per design, so it is that arm's own null.
_REAL = sorted({c for _n, b, o in HEADLINE_PAIRS + ABLATION_PAIRS for c in (b, o)} |
               {"grad_diffusion", "base_iglm", "base_diffab", "base_rfantibody",
                "base_diffusion_s", "ours_esmif_live", "obj_wt", "obj_mpnn", "obj_mint",
                "base_abmpnn_ar_t1", "base_proteinmpnn_ar_t1", "base_antifold_t1",
                "base_esmif_ar_t1"})
CELLS = sorted(set(_REAL) | {f"rand_{c}" for c in _REAL})


def load_panel():
    rows = {}
    with open(f"{G}/cand/panel_expanded.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["target"]] = r
    return rows


def eval_variants(target):
    """Held-out eval variant names; identical across cells because every arm mirrors the same
    template panel, so read whichever design JSON is present."""
    for cell in CELLS:
        for rep in REPS:
            p = f"{DESIGNS}/{target}_designs_{cell}_r{rep}.json"
            if os.path.exists(p):
                d = json.load(open(p))
                ev = [v["name"] for v in d.get("antigen_variants", []) if v.get("split") == "eval"]
                if ev:
                    return ev
    return []


def pool_worsts(target, cell, ev):
    """{rep: [per-design worst-case H3-DDG over the held-out eval variants]}.

    Kept PER REP rather than flattened. best-of-pool is monotone in pool size, so an arm with 3
    graded reps (~18 designs) beats an arm with 1 (~6) for reasons that have nothing to do with
    the objective -- exactly the artefact that made a steered arm look 9 points WORSE than its
    baseline while grading was still in flight. Paired comparisons intersect the reps available
    to BOTH arms. H3-DDG is signed vs native (native == 0); more negative = tighter.
    """
    out = {}
    for rep in REPS:
        f = f"{H3DIR}/{target}__{cell}_r{rep}.json"
        if not os.path.exists(f):
            continue
        try:
            grade = json.load(open(f)).get("h3ddg_grade", {})
        except Exception:
            continue
        w = []
        for _design, dvals in grade.items():
            vals = [dvals[n] for n in ev if n in dvals]
            # FULL coverage required: a worst-case taken over a SUBSET of the held-out variants
            # is optimistically biased (fewer chances to look bad), and coverage differs slightly
            # between arms, so partial cells would tilt the comparison.
            if len(vals) == len(ev) and vals:
                w.append(max(vals))            # worst-case = least favourable held-out variant
        if w:
            out[rep] = w
    return out


def cost_of(target, cell):
    """Mean per-run cost across reps for this (target, cell)."""
    walls, cudas, rev, dsn = [], [], [], []
    for rep in REPS:
        p = f"{DESIGNS}/{target}_designs_{cell}_r{rep}.cost.json"
        if not os.path.exists(p):
            continue
        try:
            c = json.load(open(p))
        except Exception:
            continue
        walls.append(c.get("wall_s", 0.0))
        cudas.append(c.get("cuda_s", 0.0))
        rev.append(c.get("counters", {}).get("reward_evals", 0))
        dsn.append(c.get("counters", {}).get("designs_out", 0))
    if not walls:
        return None
    return {"wall_s": st.mean(walls), "cuda_s": st.mean(cudas),
            "reward_evals": st.mean(rev), "designs": st.mean(dsn)}


def boot_rate(by_cluster, n_boot=N_BOOT, seed=0):
    """Cluster-aware rate + 95% CI, resampling ANTIGEN CLUSTERS with replacement.

    The point estimate is the mean over clusters of each cluster's own mean, so every
    independent antigen counts once regardless of how many PDB entries it contributed. It is
    DETERMINISTIC -- an earlier version drew one random target per cluster, which made the point
    estimate stochastic and inconsistent with the rates reported beside it (two arms could show
    identical rates yet a non-zero paired delta).
    """
    rng = random.Random(seed)
    clusters = list(by_cluster)
    if not clusters:
        return None, None, None
    cmeans = {c: st.mean(by_cluster[c]) for c in clusters}
    point = st.mean(cmeans.values())
    reps = []
    for _ in range(n_boot):
        pick = [cmeans[rng.choice(clusters)] for _ in clusters]
        reps.append(st.mean(pick))
    reps.sort()
    return (100 * point, 100 * reps[int(0.025 * len(reps))], 100 * reps[int(0.975 * len(reps))])


def cluster_rate(panel, targets_vals):
    """Deterministic cluster-aware rate from {target: 0/1}, matching boot_rate's estimator."""
    by = {}
    for t, v in targets_vals.items():
        by.setdefault(panel[t]["ag_cluster"], []).append(v)
    return 100 * st.mean([st.mean(v) for v in by.values()]) if by else float("nan")


def main():
    panel = load_panel()
    targets = list(panel)

    # per (target, cell): best-of-pool margin, one-shot margin, n graded designs
    cells: dict[tuple[str, str], dict] = {}
    evs = {}
    for t in targets:
        ev = eval_variants(t)
        evs[t] = ev
        if not ev:
            continue
        for cell in CELLS:
            byrep = pool_worsts(t, cell, ev)
            if not byrep:
                continue
            worsts = [x for w in byrep.values() for x in w]
            cells[(t, cell)] = {
                "by_rep": byrep,
                "reps": sorted(byrep),
                "best_of_pool": min(worsts),          # best design across all graded reps
                "one_shot": st.mean(worsts),          # expected single draw
                "n_designs": len(worsts),
                "cost": cost_of(t, cell),
            }

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["target", "ag_cluster", "vh_cluster", "in_old_51", "cell",
                    "best_of_pool_margin", "one_shot_margin", "n_designs",
                    "wall_s", "cuda_s", "reward_evals"])
        for (t, cell), v in sorted(cells.items()):
            c = v["cost"] or {}
            w.writerow([t, panel[t]["ag_cluster"], panel[t]["vh_cluster"], panel[t]["in_old_51"],
                        cell, f"{v['best_of_pool']:.4f}", f"{v['one_shot']:.4f}", v["n_designs"],
                        f"{c.get('wall_s', ''):.1f}" if c else "",
                        f"{c.get('cuda_s', ''):.1f}" if c else "",
                        f"{c.get('reward_evals', ''):.0f}" if c else ""])

    L = []
    L.append("# Generator zoo — master matrix (H3-DDG grader, held-out escape panel)\n")
    L.append("Reference = the NATIVE affinity-matured antibody (H3-DDG margin 0 by construction). "
             "`beat-native` = best-of-pool worst-case margin < 0 on variants never steered "
             "against. CI = 95% bootstrap over ANTIGEN CLUSTERS (mmseqs-70), so N is independent "
             "antigens, not PDB entries.\n")
    n_graded = len({t for (t, _c) in cells})
    L.append(f"Panel: {len(targets)} targets, "
             f"{len({panel[t]['ag_cluster'] for t in targets})} antigen clusters. "
             f"Graded so far: {n_graded} targets.\n")

    L.append("\n## Beat-native by arm\n")
    L.append("Rows flagged \u26a0\ufe0f are de-novo generators: they design a new CDR-H3 backbone "
             "(RFantibody also re-docks the Fv), so grading their SEQUENCE on the native complex "
             "understates them by construction and their row should not be read as their "
             "performance. Their fair number comes from `grade_own_structure.py`.\n")
    L.append("| arm | beat-native (best-of-pool) | 95% CI | n targets | mean margin | "
             "one-shot beat-native | wall s/run | reward evals/run |")
    L.append("|---|---|---|---|---|---|---|---|")
    summary = {}
    for cell in CELLS:
        pts = {t: cells[(t, cell)] for t in targets if (t, cell) in cells}
        if not pts:
            continue
        by_cluster: dict[str, list] = {}
        for t, v in pts.items():
            by_cluster.setdefault(panel[t]["ag_cluster"], []).append(1.0 if v["best_of_pool"] < 0 else 0.0)
        rate, lo, hi = boot_rate(by_cluster)
        one = 100 * st.mean([1.0 if v["one_shot"] < 0 else 0.0 for v in pts.values()])
        mm = st.mean([v["best_of_pool"] for v in pts.values()])
        costs = [v["cost"] for v in pts.values() if v["cost"]]
        wall = st.mean([c["wall_s"] for c in costs]) if costs else float("nan")
        rev = st.mean([c["reward_evals"] for c in costs]) if costs else float("nan")
        summary[cell] = dict(rate=rate, n=len(pts), mm=mm, wall=wall, rev=rev)
        flag = " ⚠️ de-novo backbone: see own-structure protocol" if cell in DENOVO else ""
        L.append(f"| `{cell}`{flag} | {rate:.1f}% | [{lo:.0f}, {hi:.0f}] | {len(pts)} | {mm:+.3f} | "
                 f"{one:.1f}% | {wall:.0f} | {rev:.0f} |")

    L.append("\n## Paired: generator + our objective vs generator alone\n")
    L.append("Same base distribution, same pool size, same panel — the only difference is the "
             "objective. Positive pp = our objective helps that generator.\n")
    L.append("| generator | base alone | + our objective | Δ pp | Δ CI | n paired | "
             "Δ mean margin | extra wall s | extra reward evals |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    def emit_pairs(pairs):
      for name, b, o in pairs:
          paired, bp, op = [], {}, {}
          for t in targets:
              if (t, b) not in cells or (t, o) not in cells:
                  continue
              shared = sorted(set(cells[(t, b)]["reps"]) & set(cells[(t, o)]["reps"]))
              if not shared:                       # never compare unequal pool sizes
                  continue
              bp[t] = min(x for r in shared for x in cells[(t, b)]["by_rep"][r])
              op[t] = min(x for r in shared for x in cells[(t, o)]["by_rep"][r])
              paired.append(t)
          if not paired:
              continue
          by_cluster: dict[str, list] = {}
          for t in paired:
              d = (1.0 if op[t] < 0 else 0.0) - (1.0 if bp[t] < 0 else 0.0)
              by_cluster.setdefault(panel[t]["ag_cluster"], []).append(d)
          delta, lo, hi = boot_rate(by_cluster)
          rb = cluster_rate(panel, {t: (1.0 if bp[t] < 0 else 0.0) for t in paired})
          ro = cluster_rate(panel, {t: (1.0 if op[t] < 0 else 0.0) for t in paired})
          dm = st.mean([op[t] - bp[t] for t in paired])
          cb = [cells[(t, b)]["cost"] for t in paired if cells[(t, b)]["cost"]]
          co = [cells[(t, o)]["cost"] for t in paired if cells[(t, o)]["cost"]]
          dw = (st.mean([c["wall_s"] for c in co]) - st.mean([c["wall_s"] for c in cb])) if cb and co else float("nan")
          dr = (st.mean([c["reward_evals"] for c in co]) - st.mean([c["reward_evals"] for c in cb])) if cb and co else float("nan")
          sig = " **" if (lo > 0 or hi < 0) else ""
          L.append(f"| {name} | {rb:.1f}% | {ro:.1f}% | {delta:+.1f}{sig} | "
                   f"[{lo:+.0f}, {hi:+.0f}] | {len(paired)} | {dm:+.3f} | {dw:+.0f} | {dr:+.0f} |")


    emit_pairs(HEADLINE_PAIRS)
    L.append("\n### Ablation: same models with a FROZEN prior on both sides\n")
    L.append("Isolates what re-conditioning the generator each denoising step buys. For the "
             "autoregressive models this is NOT the headline comparison.\n")
    L.append("| generator | base alone | + our objective | \u0394 pp | \u0394 CI | n paired | "
             "\u0394 mean margin | extra wall s | extra reward evals |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    emit_pairs(ABLATION_PAIRS)

    L.append("\n## Lift over each arm's OWN edit-distance-matched random control\n")
    L.append("Beat-native alone confounds design quality with how far the design moved from the "
             "native CDR: more edits score worse under the grader regardless of who proposed "
             "them. Each arm is therefore compared to a random control that mutates the SAME "
             "number of positions per design. LIFT = arm rate minus its own null rate, on the "
             "targets where both are graded; it is the part of the score not explained by "
             "edit distance.\n")
    L.append("Restricted to targets where ALL FOUR cells are graded (base, ours, and both "
             "nulls), so the two lifts are computed on the SAME targets and are comparable. "
             "Per-arm lift on differing subsets is not.\n")
    L.append("| generator | arm | beat-native | its random null | LIFT pp | \u0394 LIFT (ours-base) | n |")
    L.append("|---|---|---|---|---|---|---|")
    for name, b, o in HEADLINE_PAIRS:
        quad = [t for t in targets
                if all((t, c) in cells for c in (b, o, f"rand_{b}", f"rand_{o}"))]
        if len(quad) < 10:
            continue
        def rate(c):
            return cluster_rate(panel, {t: (1.0 if cells[(t, c)]["best_of_pool"] < 0 else 0.0)
                                        for t in quad})
        rb, nb = rate(b), rate(f"rand_{b}")
        ro, no = rate(o), rate(f"rand_{o}")
        # cluster bootstrap on the per-target lift difference, so delta-LIFT carries a CI
        bl = {}
        for t in quad:
            w = lambda c: 1.0 if cells[(t, c)]["best_of_pool"] < 0 else 0.0
            bl.setdefault(panel[t]["ag_cluster"], []).append(
                (w(o) - w(f"rand_{o}")) - (w(b) - w(f"rand_{b}")))
        dl, dlo, dhi = boot_rate(bl)
        L.append(f"| {name} | `{b}` | {rb:.1f}% | {nb:.1f}% | **{rb - nb:+.1f}** | | {len(quad)} |")
        sig = " **" if (dlo > 0 or dhi < 0) else ""
        L.append(f"| {name} | `{o}` | {ro:.1f}% | {no:.1f}% | **{ro - no:+.1f}** | "
                 f"**{dl:+.1f}**{sig} [{dlo:+.0f}, {dhi:+.0f}] | {len(quad)} |")

    L.append("\n## Old-51 subset vs the new targets (is the expansion changing the story?)\n")
    L.append("| arm | beat-native on old 51 | beat-native on the 149 new | n old | n new |")
    L.append("|---|---|---|---|---|")
    for cell in CELLS:
        old = [cells[(t, cell)] for t in targets
               if (t, cell) in cells and panel[t]["in_old_51"] == "1"]
        new = [cells[(t, cell)] for t in targets
               if (t, cell) in cells and panel[t]["in_old_51"] == "0"]
        if not old and not new:
            continue
        ro = 100 * st.mean([1.0 if v["best_of_pool"] < 0 else 0.0 for v in old]) if old else float("nan")
        rn = 100 * st.mean([1.0 if v["best_of_pool"] < 0 else 0.0 for v in new]) if new else float("nan")
        L.append(f"| `{cell}` | {ro:.1f}% | {rn:.1f}% | {len(old)} | {len(new)} |")

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\n-> {OUT_MD}\n-> {OUT_CSV}")


if __name__ == "__main__":
    main()
