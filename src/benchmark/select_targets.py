#!/usr/bin/env python3
"""Phase 1 step B/C: cluster candidate targets and greedily select the expanded panel.

The old 51-target panel's weakness was not its size but its REDUNDANCY: 51 targets collapsed
to ~31 distinct antigen clusters at mmseqs-70, so the effective N behind every beat-native
number was 31. Selection here therefore maximises DISTINCT ANTIGEN CLUSTERS first and raw
count second — a panel of 200 targets spread over 200 clusters is worth far more than 400
targets over 60.

  1. mmseqs easy-cluster antigen sequences at 70% identity  -> antigen cluster id
  2. mmseqs easy-cluster VH sequences at 90% identity       -> antibody cluster id
  3. greedy round-robin over antigen clusters, best-quality target first within each cluster,
     never two targets sharing BOTH an antigen cluster and an antibody cluster
  4. the existing 51 panel targets are force-included (so the new panel is a strict superset
     and every historical number stays directly comparable)

Output: panel_expanded.txt (one target per line) + panel_expanded.tsv (with cluster columns).
"""
from __future__ import annotations
import argparse, csv, os, subprocess, sys
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

G = f"{VAAD_ROOT}/scratch/genzoo"
MMSEQS = f"{VAAD_ROOT}/tools/mmseqs/bin/mmseqs"
OLD_PANEL = f"{VAAD_ROOT}/scratch/svdd_screen/panel_targets.txt"


def mmseqs_cluster(seqs: dict, tag: str, min_id: float, cov: float = 0.8) -> dict:
    """seqs: {target: sequence} -> {target: cluster_rep}. Uses mmseqs easy-cluster."""
    wd = os.path.join(G, "cand", f"mmseqs_{tag}")
    os.makedirs(wd, exist_ok=True)
    fa = os.path.join(wd, "in.fasta")
    with open(fa, "w") as fh:
        for t, s in seqs.items():
            fh.write(f">{t}\n{s}\n")
    pre = os.path.join(wd, "clu")
    cmd = [MMSEQS, "easy-cluster", fa, pre, os.path.join(wd, "tmp"),
           "--min-seq-id", str(min_id), "-c", str(cov), "--cov-mode", "0", "-v", "1"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    out = {}
    with open(pre + "_cluster.tsv") as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            out[mem] = rep
    return out


def quality(row) -> float:
    """Prefer targets with a rich escape library, a plausible CDR-H3, and a paired H+L Fv."""
    nv = min(int(row["n_variants"]), 3000) / 3000.0
    npos = min(int(row["n_positions"]), 30) / 30.0
    cdr = int(row["cdr_len"])
    cdr_ok = 1.0 if 8 <= cdr <= 22 else (0.4 if 5 <= cdr <= 25 else 0.0)
    paired = 1.0 if int(row["n_ab_chains"]) == 2 else 0.6
    return 2.0 * npos + 1.0 * nv + 1.5 * cdr_ok + 0.5 * paired


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--ag-id", type=float, default=0.7)
    ap.add_argument("--vh-id", type=float, default=0.9)
    args = ap.parse_args()

    rows = {}
    with open(os.path.join(G, "cand", "candidates.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows[row["target"]] = row
    print(f"[select] {len(rows)} eligible candidates", flush=True)

    ag = mmseqs_cluster({t: r["ag_seq"] for t, r in rows.items()}, "antigen", args.ag_id)
    vh = mmseqs_cluster({t: r["vh_seq"] for t, r in rows.items()}, "vh", args.vh_id)
    print(f"[select] antigen clusters @{args.ag_id}: {len(set(ag.values()))}; "
          f"VH clusters @{args.vh_id}: {len(set(vh.values()))}", flush=True)

    old = [t.strip() for t in open(OLD_PANEL) if t.strip()]
    old_ok = [t for t in old if t in rows]
    print(f"[select] existing panel: {len(old)} targets, {len(old_ok)} pass eligibility", flush=True)

    # bucket candidates by antigen cluster, best-first
    buckets: dict[str, list] = {}
    for t, r in rows.items():
        buckets.setdefault(ag.get(t, t), []).append(t)
    for k in buckets:
        buckets[k].sort(key=lambda t: -quality(rows[t]))

    chosen, taken_pairs = [], set()

    def take(t, force=False):
        key = (ag.get(t, t), vh.get(t, t))
        if t in chosen or (key in taken_pairs and not force):
            return False
        chosen.append(t)
        taken_pairs.add(key)
        return True

    # Force-include unconditionally (bypassing the near-duplicate rule) so the new panel is a
    # STRICT superset of the old 51 and every historical number stays directly comparable.
    # A handful of the old 51 are near-duplicates of each other -- that redundancy is exactly
    # what cluster-aware aggregation downstream is for, so we keep them and let the analysis
    # collapse them rather than silently dropping targets from the published panel.
    for t in old_ok:
        take(t, force=True)
    n_forced = len(chosen)

    order = sorted(buckets, key=lambda k: -max(quality(rows[t]) for t in buckets[k]))
    depth = 0
    while len(chosen) < args.n:
        added = 0
        for k in order:
            if len(chosen) >= args.n:
                break
            if depth < len(buckets[k]) and take(buckets[k][depth]):
                added += 1
        depth += 1
        if added == 0:
            break

    outtsv = os.path.join(G, "cand", "panel_expanded.tsv")
    with open(outtsv, "w") as fh:
        fh.write("target\tag_cluster\tvh_cluster\tn_variants\tn_positions\tcdr_len\tcdr_seq\tin_old_51\n")
        for t in chosen:
            r = rows[t]
            fh.write(f"{t}\t{ag.get(t,t)}\t{vh.get(t,t)}\t{r['n_variants']}\t{r['n_positions']}\t"
                     f"{r['cdr_len']}\t{r['cdr_seq']}\t{int(t in old)}\n")
    with open(os.path.join(G, "cand", "panel_expanded.txt"), "w") as fh:
        fh.write("\n".join(chosen) + "\n")

    nag = len({ag.get(t, t) for t in chosen})
    nvh = len({vh.get(t, t) for t in chosen})
    old_ag = len({ag.get(t, t) for t in old_ok})
    print(f"[select] SELECTED {len(chosen)} targets ({n_forced} carried from the old panel)")
    print(f"[select]   distinct antigen clusters: {nag}  (old panel: {old_ag})")
    print(f"[select]   distinct VH clusters     : {nvh}")
    print(f"[select] -> {outtsv}")


if __name__ == "__main__":
    main()
