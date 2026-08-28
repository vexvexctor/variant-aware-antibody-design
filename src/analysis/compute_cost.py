#!/usr/bin/env python3
"""Computational-overhead report for every arm in the generator zoo.

The method's selling point is that it is TRAINING-FREE: no fine-tuning, no reward-model
training, no distillation -- just a frozen prior, a frozen reward, and K reward calls per
denoising step. That claim is only worth anything if the alternatives' costs are measured on
the same machine with the same accounting, which is what profile_util.Profiler does for every
design, baseline and export run.

Two axes are reported because they answer different questions:
  * seconds (wall / CUDA)  -- what it actually costs on this cluster's H100s
  * reward_evals, base_forwards -- hardware-independent counts, so a reader can rescale to
    their own setup and so search baselines (GA, best-of-N) are comparable on the budget axis
    they actually spend.

Amortised one-off costs (logit-table exports) are reported separately from the marginal
per-design cost, because they are paid once per target and then reused by every arm and rep.
"""
from __future__ import annotations

import csv
import json
import os
import statistics as st
from collections import defaultdict
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

G = f"{VAAD_ROOT}/scratch/genzoo"
OUT_MD = f"{G}/OVERHEAD.md"
OUT_CSV = f"{G}/overhead_by_arm.csv"


def load(dirpath, suffix=".cost.json"):
    """Return (filename, cost-dict) pairs. The FILENAME is the authority on which matrix cell a
    run belongs to: two different cells can share a profiler tag (base_diffusion is the same
    mechanism as grad_diffusion, just with guidance-scale 0), and merging them would silently
    average an unguided arm into a guided one."""
    rows = []
    if not os.path.isdir(dirpath):
        return rows
    for fn in os.listdir(dirpath):
        if not fn.endswith(suffix):
            continue
        try:
            rows.append((fn, json.load(open(os.path.join(dirpath, fn)))))
        except Exception:
            continue
    return rows


def cell_of(fn, cost):
    """Matrix cell from the output filename: {TARGET}_designs_{CELL}_r{REP}.cost.json,
    or {TARGET}_{MODEL}_logits.npz.cost.json for the one-off exports."""
    if "_designs_" in fn:
        rest = fn.split("_designs_", 1)[1]
        return rest.rsplit("_r", 1)[0] if "_r" in rest else rest.replace(".cost.json", "")
    if "_logits.npz" in fn:
        # target names themselves contain underscores (PDB_CHAINS), so parse the model from the
        # profiler tag ("export/<model>") rather than trying to split the filename
        tag = cost.get("tag", "")
        return "export_" + (tag.split("/", 1)[1] if "/" in tag else "unknown")
    return cost.get("tag", "unknown")


def agg(pairs):
    rows = [r for _fn, r in pairs]
    def m(key, path=None):
        vals = []
        for r in rows:
            v = r.get("counters", {}).get(key) if path == "counters" else r.get(key)
            if isinstance(v, (int, float)):
                vals.append(v)
        return vals
    designs = m("designs_out", "counters") or [1]
    n_designs = st.mean(designs) if designs else 1
    wall = m("wall_s"); cuda = m("cuda_s")
    rev = m("reward_evals", "counters") or [0]
    bf = m("base_forwards", "counters") or [0]
    vram = [v for v in m("peak_vram_gb") if isinstance(v, (int, float))]
    return {
        "n_runs": len(rows),
        "wall_s": st.mean(wall) if wall else float("nan"),
        "wall_s_per_design": (st.mean(wall) / n_designs) if wall and n_designs else float("nan"),
        "cuda_s": st.mean(cuda) if cuda else float("nan"),
        "reward_evals": st.mean(rev),
        "reward_evals_per_design": st.mean(rev) / n_designs if n_designs else float("nan"),
        "base_forwards": st.mean(bf),
        "peak_vram_gb": max(vram) if vram else float("nan"),
        "designs": n_designs,
    }


def main():
    by_cell = defaultdict(list)
    for fn, r in load(f"{G}/designs"):
        if fn.startswith("SMOKE"):
            continue
        by_cell[cell_of(fn, r)].append((fn, r))
    exports = defaultdict(list)
    for fn, r in load(f"{G}/npz"):
        exports[cell_of(fn, r)].append((fn, r))

    rows = []
    for cell, rs in sorted(by_cell.items()):
        a = agg(rs)
        a["cell"] = cell
        rows.append(a)

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cell", "n_runs", "designs_per_run", "wall_s", "wall_s_per_design", "cuda_s",
                    "reward_evals", "reward_evals_per_design", "base_forwards", "peak_vram_gb"])
        for a in rows:
            w.writerow([a["cell"], a["n_runs"], f"{a['designs']:.1f}", f"{a['wall_s']:.2f}",
                        f"{a['wall_s_per_design']:.2f}", f"{a['cuda_s']:.2f}",
                        f"{a['reward_evals']:.0f}", f"{a['reward_evals_per_design']:.0f}",
                        f"{a['base_forwards']:.0f}", f"{a['peak_vram_gb']:.2f}"])

    L = ["# Computational overhead by arm\n",
         "All runs on one H100 per task, identical accounting via `profile_util.Profiler`. "
         "TRAINING COST IS ZERO FOR EVERY ARM HERE: no arm fine-tunes, trains a reward model, or "
         "distills anything -- the priors and the BA-DDG reward are frozen external checkpoints. "
         "So the whole comparison is inference-time cost.\n",
         "\n## Marginal cost per run (one target, one rep, pool of 6 designs)\n",
         "| arm | n runs | wall s | wall s / design | CUDA s | reward evals | reward evals / design | base forwards | peak VRAM GB |",
         "|---|---|---|---|---|---|---|---|---|"]
    for a in sorted(rows, key=lambda x: x["wall_s"]):
        aw, awd, ac = a["wall_s"], a["wall_s_per_design"], a["cuda_s"]
        L.append(f"| `{a['cell']}` | {a['n_runs']} | {aw:.2f} | {awd:.3f} | "
                 f"{ac:.2f} | {a['reward_evals']:.0f} | {a['reward_evals_per_design']:.0f} | "
                 f"{a['base_forwards']:.0f} | {a['peak_vram_gb']:.2f} |")

    if exports:
        L.append("\n## One-off amortised cost: per-target logit-table export\n")
        L.append("Paid once per target, then reused by every arm and every rep that uses that base. "
                 "For the table-based bases (AntiFold, and the frozen-npz ablation arms) this export "
                 "IS the generator's real inference cost -- the 0.0x s 'baseline' below is only the "
                 "cost of drawing samples from an already-computed table, so read the two together.\n")
        L.append("| export | n targets | wall s | CUDA s |")
        L.append("|---|---|---|---|")
        for cell, rs in sorted(exports.items()):
            a = agg(rs)
            L.append(f"| `{cell}` | {a['n_runs']} | {a['wall_s']:.2f} | {a['cuda_s']:.2f} |")

    # the headline overhead statement: our objective's cost ON TOP of each generator
    L.append("\n## Overhead of adding our objective to a generator\n")
    L.append("| generator | baseline wall s | + our objective wall s | multiplier | added reward evals |")
    L.append("|---|---|---|---|---|")
    pairs = [("antifold", "base_antifold", "ours_antifold"),
             ("ProteinMPNN", "base_proteinmpnn_ar", "ours_proteinmpnn_live"),
             ("AbMPNN", "base_abmpnn_ar", "ours_abmpnn_live"),
             ("ESM-IF1", "base_esmif_ar", "ours_esmif"),
             ("our diffusion", "base_diffusion", "ours_diffusion")]
    idx = {a["cell"]: a for a in rows}
    for name, b, o in pairs:
        if b in idx and o in idx:
            bb, oo = idx[b], idx[o]
            mult = oo["wall_s"] / bb["wall_s"] if bb["wall_s"] else float("nan")
            L.append(f"| {name} | {bb['wall_s']:.2f} | {oo['wall_s']:.2f} | {mult:.0f}x | "
                     f"{oo['reward_evals'] - bb['reward_evals']:.0f} |")

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\n-> {OUT_MD}\n-> {OUT_CSV}")


if __name__ == "__main__":
    main()
