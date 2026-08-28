#!/usr/bin/env python3
"""Cluster steering (the end goal): design a CDR robust across a target's weighted
antigen-mutant cluster, using direct differentiable MPNN (no distillation).

For one target complex it builds one `feats` per antigen mutant — the antigen block
of `S` (positions where chain_encoding > #antibody chains) replaced by the mutant
antigen sequence from clustered.csv — plus its omega. Then it compares three CDR
designs by scoring each against EVERY cluster member (ground-truth MPNN log-ratio):

  plain   — diffusion.sample()
  single  — sample_nos_mpnn_direct() steered on the WT complex only
  cluster — sample_nos_mpnn_cluster_direct() steered on the weighted mutant set

Robustness = MEAN and WORST-CASE log-ratio across the cluster. The end-goal claim
is that `cluster` improves the worst-case (and/or mean) over `single`.

With `--holdout-frac F` the cluster is split: `cluster` steers on the training
members only, and ALL three designs are scored on the HELD-OUT members — testing
generalization to UNSEEN escape variants (the convincing robustness claim). F=0
(default) scores on the full cluster = the in-cluster behavior above.

Usage
-----
    PYTHONPATH=src python3 scripts/nos_diffusion/cluster_steer.py \
        --target 1A14_HLN --clustered-csv $RESULTS/clustered.csv \
        --nos-checkpoint $RESULTS/finetune/nos_4layer/finetune_best.pt \
        --features-dir $SCRATCH/complex_features \
        --proteinmpnn-dir $LYRA/vendor/ProteinMPNN \
        --proteinmpnn-weights $LYRA/vendor/ProteinMPNN/vanilla_model_weights/v_48_002.pt \
        --cdr-positions 95 96 97 98 99 100 101 102 103 104 105 106 \
        --max-cluster 16 --n-reps 3 --guidance-scale 1.0 [--holdout-frac 0.3 --seed 0]
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import statistics as stats
import sys

import torch

AA = "ACDEFGHIKLMNPQRSTVWYX"
N_AA = 20

_3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
    "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
    "TRP": "W", "TYR": "Y",
}


def parse_antigen_residues(pdb_path, antigen_chains):
    """Ordered [(chain, resnum_1based|None, icode, aa1)] matching parse_PDB's layout:
    residues filled by number from min..max (missing -> 'X' gap), insertion codes
    sorted within a number, residue present on ANY atom line, MSE->MET. resnum is
    None for gap positions. This must mirror parse_PDB_biounits so block indices
    line up with feats['S']."""
    ordered = []
    for chain in sorted(antigen_chains):
        seq: dict[int, dict[str, str]] = {}
        min_resn, max_resn = 10**6, -10**6
        with open(pdb_path, errors="ignore") as fh:
            for line in fh:
                if line[:6] == "HETATM" and line[17:20] == "MSE":
                    line = line.replace("HETATM", "ATOM  ").replace("MSE", "MET")
                if line[:4] != "ATOM" or line[21:22] != chain:
                    continue
                field = line[22:27].strip()
                if not field:
                    continue
                if field[-1].isalpha():
                    resa = field[-1]
                    try:
                        resn = int(field[:-1]) - 1
                    except ValueError:
                        continue
                else:
                    resa = ""
                    try:
                        resn = int(field) - 1
                    except ValueError:
                        continue
                min_resn, max_resn = min(min_resn, resn), max(max_resn, resn)
                seq.setdefault(resn, {}).setdefault(resa, line[17:20].strip())
        if max_resn < min_resn:
            continue
        for resn in range(min_resn, max_resn + 1):
            if resn in seq:
                for k in sorted(seq[resn]):
                    ordered.append((chain, resn + 1, k, _3TO1.get(seq[resn][k], "X")))
            else:
                ordered.append((chain, None, "", "X"))   # gap fill, like parse_PDB
    return ordered


def parse_mutation(tok: str):
    """'N:A369I' -> (chain, wt, resnum, icode, mut), or None."""
    tok = tok.strip()
    if ":" not in tok:
        return None
    chain, rest = tok.split(":", 1)
    m = re.match(r"^([A-Z])(-?\d+)([A-Za-z]?)([A-Z])$", rest.strip())
    if not m:
        return None
    wt, resnum, icode, mut = m.groups()
    return chain.strip(), wt, int(resnum), icode.upper(), mut


def parse_chains(complex_id: str):
    """Same convention as featurize_complexes.py: antibody vs antigen chains."""
    chains_str = complex_id.split("_")[-1]
    if len(chains_str) >= 3 and chains_str[1] == "L":
        antibody = [chains_str[0], chains_str[1]]
    else:
        antibody = [chains_str[0]]
    return antibody


def to_seq(toks) -> str:
    return "".join(AA[int(t)] if int(t) < N_AA else "X" for t in toks)


def tokenize(seq: str, device):
    return torch.tensor([AA.index(a) if a in AA[:N_AA] else N_AA for a in seq],
                        dtype=torch.long, device=device)


def cdr_h3_positions(base, n_ab_chains, min_len=3, max_len=30):
    """Derive CDR-H3 flat indices for the heavy chain from packed features, by anchoring
    on the conserved FR3 Cys and the FR4 'WGxG' motif: CDR-H3 = (Cys+1 .. W-1).

    Replaces the hardcoded flat slice (95..106) which mis-registered across targets: the
    fixed window (a) included the conserved Cys anchor on some targets (e.g. 1A14, where the
    VH has one extra FR3 residue so flat-95 lands on the Cys) and (b) truncated the loop's
    C-terminal residues on all targets (real CDR-H3 here is 14-15 long, not 12). Anchoring
    on Cys->WGxG yields the full, consistently-registered loop per target.

    The heavy chain is among the first `n_ab_chains` chains; the light-chain FR4 is 'FGxG',
    so the 'WGxG' anchor uniquely selects the heavy chain. Returns flat indices into
    feats['S'][0].
    """
    S = base["S"][0].tolist()
    ce = base["chain_encoding_all"][0].tolist()
    seq = to_seq(S)
    for cval in sorted(set(ce))[:n_ab_chains]:
        idx = [i for i in range(len(S)) if ce[i] == cval]
        cseq = "".join(seq[i] for i in idx)
        wmatches = [m.start() for m in re.finditer(r"WG.G", cseq)]
        cmatches = [m.start() for m in re.finditer("C", cseq)]
        if not wmatches or not cmatches:
            continue
        w = wmatches[-1]                                    # VH FR4 (last WGxG in the V-region)
        c = max((x for x in cmatches if x < w), default=None)
        if c is None or not (min_len <= w - c - 1 <= max_len):
            continue
        return [idx[i] for i in range(c + 1, w)]            # Cys+1 .. W-1
    raise ValueError("cdr_h3_positions: no heavy-chain Cys->WGxG anchor found")


def gt_logratio(mpnn_model, feats, cdr_positions, toks, native, device):
    """Ground-truth MPNN complex log-ratio of a CDR vs native, for given feats."""
    S = feats["S"].clone()
    S[0, cdr_positions] = toks
    chain_M = torch.zeros_like(feats["chain_M"])
    chain_M[0, cdr_positions] = 1
    randn = torch.randn(1, feats["X"].shape[1], device=device)
    with torch.no_grad():
        logp = mpnn_model.conditional_probs(
            feats["X"], S, feats["mask"], chain_M, feats["residue_idx"],
            feats["chain_encoding_all"], randn, backbone_only=False)[0]
    return float(sum(logp[p, int(toks[k])] - logp[p, int(native[k])]
                     for k, p in enumerate(cdr_positions)))


def _refold_name(mutation_string):
    """refold-panel filename stem for a variant: same sanitizer as refold stage3's `safe()`."""
    return re.sub(r"[^A-Za-z0-9]", "_", (mutation_string or "").strip())[:60]


def _featurize_refold(pdb_path, ab_letters, visible_letters, device):
    """Featurize a re-folded variant PDB -> feats dict whose X is the RE-FOLDED backbone.
    Mirrors the refold panel's stage3 loader (parse_PDB -> tied_featurize; designed=antibody,
    visible=antigen, which must be DISJOINT). Returns None on any parse failure. The vendor
    ProteinMPNN dir must already be on sys.path (both build_cluster callers add it)."""
    from protein_mpnn_utils import parse_PDB, tied_featurize
    pdb_list = parse_PDB(pdb_path)
    if not pdb_list:
        return None
    nm = pdb_list[0]["name"]
    r = tied_featurize(pdb_list, "cpu", {nm: (list(ab_letters), list(visible_letters))},
                       None, None, None, None, None)
    X, S, mask, chain_M, residue_idx, chain_enc = r[0], r[1], r[2], r[4], r[12], r[5]
    return {"X": X.to(device), "S": S.to(device), "mask": mask.to(device),
            "chain_M": chain_M.to(device), "residue_idx": residue_idx.to(device),
            "chain_encoding_all": chain_enc.to(device)}


def build_cluster(base_feats, target, rows, max_cluster, pdb_path, device, pick="top", seed=0,
                  refold_dir=None, cdr_positions=None):
    """Build [(feats_i, omega_i)] by stamping each mutant's point mutations into the
    antigen block of S, mapping PDB residue numbers -> block index via the PDB.

    `pick` selects which `max_cluster` mutants to take from `rows`:
      - "top"    : highest omega (WT-proximal) — the default steering panel
      - "bottom" : lowest omega (most divergent) — held-out divergence eval set
      - "random" : random sample at `seed` — i.i.d. held-out floor

    `refold_dir` (needs `cdr_positions`): per-variant re-folded PDBs from the refold panel
    (`{refold_dir}/{_refold_name(mutation_string)}.pdb`). When a variant's refold exists and
    is index-aligned (same length AND native CDR at `cdr_positions`, mirroring stage3's guards),
    that member's backbone `X` is REPLACED by its own re-folded structure — fixing the
    frozen-backbone degeneracy where every escape variant shared the native backbone (only S
    differed). S/topology stay frozen-aligned; only X carries the per-variant geometry. Misses
    fall back to the frozen backbone and are logged.
    """
    antibody = parse_chains(target)
    n_ab = len(antibody)
    chains_str = target.split("_")[-1]
    antigen_chains = [c for c in chains_str if c not in antibody]

    enc = base_feats["chain_encoding_all"][0]
    antigen_idx = (enc > n_ab).nonzero(as_tuple=True)[0]   # flat S indices of the antigen block
    wt_block = to_seq(base_feats["S"][0, antigen_idx])

    ordered = parse_antigen_residues(pdb_path, antigen_chains)
    parsed_seq = "".join(aa for _, _, _, aa in ordered)
    aligned = (len(ordered) == len(antigen_idx) and parsed_seq == wt_block)
    res_to_idx = {(c, rn, ic): i for i, (c, rn, ic, _) in enumerate(ordered) if rn is not None}

    _omega = lambda r: float(r.get("omega_normalized", 0) or 0)  # noqa: E731
    if pick == "top":
        rows = sorted(rows, key=_omega, reverse=True)[:max_cluster]
    elif pick == "bottom":
        rows = sorted(rows, key=_omega)[:max_cluster]
    elif pick == "random":
        import random as _random
        rows = list(rows)
        _random.Random(seed).shuffle(rows)
        rows = rows[:max_cluster]
    else:
        raise ValueError(f"unknown pick={pick!r} (top|bottom|random)")
    use_refold = refold_dir is not None and cdr_positions is not None
    if use_refold:
        base_L = base_feats["S"].shape[1]
        native_cdr = to_seq(base_feats["S"][0, cdr_positions])
        n_refold = n_refold_miss = 0

    out, labels = [], []
    n_parse = n_map = n_wt = 0
    for r in rows:
        toks = [parse_mutation(t) for t in r.get("mutation_string", "").split(";") if t.strip()]
        if not toks or any(t is None for t in toks):
            n_parse += 1
            continue
        f = {k: v.clone() for k, v in base_feats.items()}
        ok = True
        for chain, wt, resnum, icode, mut in toks:
            bi = res_to_idx.get((chain, resnum, icode))
            if bi is None:
                ok = False; n_map += 1; break
            if wt_block[bi] != wt:          # WT-residue verification (alignment guard)
                ok = False; n_wt += 1; break
            f["S"][0, antigen_idx[bi]] = AA.index(mut) if mut in AA[:N_AA] else N_AA
        if ok:
            if use_refold:                  # swap in this variant's OWN re-folded backbone (X only)
                rpdb = f"{refold_dir}/{_refold_name(r.get('mutation_string', ''))}.pdb"
                rf = _featurize_refold(rpdb, antibody, antigen_chains, device) \
                    if os.path.exists(rpdb) else None
                if (rf is not None and rf["S"].shape[1] == base_L
                        and to_seq(rf["S"][0, cdr_positions]) == native_cdr):
                    f["X"] = rf["X"]
                    n_refold += 1
                else:
                    n_refold_miss += 1
            out.append((f, float(r.get("omega_normalized", 1.0) or 1.0)))
            labels.append((r.get("mutation_string", "").strip(), toks))   # name + parsed point mutations
    if use_refold:
        print(f"[build_cluster] refold backbones: {n_refold} swapped, {n_refold_miss} "
              f"fell back to frozen (missing/misaligned) from {refold_dir}", flush=True)
    return out, len(antigen_idx), aligned, (n_parse, n_map, n_wt), labels


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, help="complex_id to design for, e.g. 1A14_HLN")
    ap.add_argument("--clustered-csv", required=True)
    ap.add_argument("--aacdb-dir", required=True, help="dir with {complex_id}.pdb (antigen residue numbers)")
    ap.add_argument("--nos-checkpoint", required=True, help="BASE diffusion checkpoint")
    ap.add_argument("--features-dir", required=True)
    ap.add_argument("--proteinmpnn-dir", required=True)
    ap.add_argument("--proteinmpnn-weights", required=True)
    ap.add_argument("--cdr-positions", nargs="+", type=int, required=True)
    ap.add_argument("--max-cluster", type=int, default=16, help="top-omega mutants to steer over")
    ap.add_argument("--n-reps", type=int, default=3)
    ap.add_argument("--guidance-scale", type=float, default=1.0)
    ap.add_argument("--holdout-frac", type=float, default=0.0,
                    help="fraction of cluster held out for scoring; cluster steers on the rest (0=score on full cluster)")
    ap.add_argument("--seed", type=int, default=0, help="rng seed for the holdout split")
    ap.add_argument("--seeds", nargs="+", type=int, default=None,
                    help="sweep these holdout-split seeds (overrides --seed); prints a per-seed table + win tally")
    ap.add_argument("--n-layers", type=int, default=4)
    ap.add_argument("--hidden-dim", type=int, default=128)
    args = ap.parse_args()

    sys.path.insert(0, args.proteinmpnn_dir)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    from lyra_mutants.Model.diffusion import CDRMaskedDiffusion
    from lyra_mutants.Model.mpnn_adapter import ProteinMPNNAdapter

    n_cdr = len(args.cdr_positions)
    diffusion = CDRMaskedDiffusion(
        n_cdr=n_cdr, hidden_dim=args.hidden_dim, n_layers=args.n_layers,
        binder_dim=128, scorer_dim=128).to(device)
    ckpt = torch.load(args.nos_checkpoint, map_location=device, weights_only=True)
    diffusion.load_state_dict(ckpt["diffusion_state_dict"], strict=False)
    diffusion.eval()

    adapter = ProteinMPNNAdapter(
        args.proteinmpnn_dir, args.proteinmpnn_weights, args.features_dir, args.cdr_positions).to(device)
    for p in adapter.model.parameters():
        p.requires_grad_(False)

    base = adapter._features(args.target)
    struct = adapter.encoder_hidden_states(args.target)
    if base is None or struct is None:
        raise SystemExit(f"No features for {args.target} in {args.features_dir}")
    struct = struct.to(device).detach()
    native = base["S"][0, args.cdr_positions].long().clamp(max=N_AA - 1)

    rows = [r for r in csv.DictReader(open(args.clustered_csv, newline=""))
            if r["complex_id"] == args.target]
    pdb_path = f"{args.aacdb_dir}/{args.target}.pdb"
    cluster, ag_len, aligned, (n_parse, n_map, n_wt), _labels = build_cluster(
        base, args.target, rows, args.max_cluster, pdb_path, device)
    print(f"target {args.target}: {len(rows)} mutants in csv, antigen length {ag_len}, "
          f"PDB↔feats aligned={aligned}")
    print(f"cluster size {len(cluster)}  (skipped: parse {n_parse}, resnum-not-found {n_map}, wt-mismatch {n_wt})")
    if not aligned:
        print("WARNING: parsed PDB antigen seq != feats antigen block — resnum mapping may be off.")
    if not cluster:
        raise SystemExit("empty cluster — check aacdb-dir / mutation_string / alignment")

    # ---- evaluation (optionally swept over holdout-split seeds) ----
    import random
    do_holdout = args.holdout_frac > 0
    seeds = (args.seeds if args.seeds else [args.seed]) if do_holdout else [None]

    def split_cluster(seed):
        """-> (steer_members, score_members). No holdout: both = full cluster."""
        if not do_holdout:
            return cluster, cluster
        idx = list(range(len(cluster)))
        random.Random(seed).shuffle(idx)
        n_test = max(1, round(len(cluster) * args.holdout_frac))
        if len(cluster) - n_test < 1:
            raise SystemExit(f"holdout-frac {args.holdout_frac} leaves no training members (cluster size {len(cluster)})")
        return [cluster[i] for i in idx[n_test:]], [cluster[i] for i in idx[:n_test]]

    def eval_split(steer_members, score_members):
        """Mean over reps of (mean, worst-case) log-ratio across score_members, per design."""
        def score_across(toks):
            s = [gt_logratio(adapter.model, f, args.cdr_positions, toks, native, device) for f, _ in score_members]
            return stats.mean(s), min(s)
        res = {"plain": [], "single": [], "cluster": []}
        for rep in range(args.n_reps):
            plain = diffusion.sample(struct)[:n_cdr]
            single = diffusion.sample_nos_mpnn_direct(struct, adapter.model, base, args.cdr_positions, args.guidance_scale)[:n_cdr]
            clust = diffusion.sample_nos_mpnn_cluster_direct(struct, adapter.model, steer_members, args.cdr_positions, args.guidance_scale)[:n_cdr]
            for name, toks in (("plain", plain), ("single", single), ("cluster", clust)):
                res[name].append(score_across(toks))
            print(f"    rep {rep+1}/{args.n_reps} done", flush=True)
        return {name: (stats.mean([m for m, _ in res[name]]), stats.mean([w for _, w in res[name]]))
                for name in res}

    per_seed, win_worst, win_mean = [], 0, 0
    for seed in seeds:
        steer_members, score_members = split_cluster(seed)
        tag = (f"seed {seed}: steer {len(steer_members)} / score {len(score_members)} HELD-OUT"
               if do_holdout else f"full cluster ({len(cluster)} variants)")
        print(f"\n[{tag}]")
        agg = eval_split(steer_members, score_members)
        per_seed.append((seed, len(score_members), agg))
        win_worst += agg["cluster"][1] > agg["single"][1]
        win_mean += agg["cluster"][0] > agg["single"][0]

    scope = "HELD-OUT antigen variants" if do_holdout else "antigen variants"
    print(f"\n=== cluster steering: WORST-CASE log-ratio across {scope} ===")
    print(f"{'seed':>4} {'n':>3} | {'plain':>8} {'single':>8} {'cluster':>8} | {'Δworst':>7} {'Δmean':>7}  win")
    for seed, n, agg in per_seed:
        dw = agg["cluster"][1] - agg["single"][1]
        dm = agg["cluster"][0] - agg["single"][0]
        print(f"{str(seed):>4} {n:>3} | {agg['plain'][1]:>8.2f} {agg['single'][1]:>8.2f} {agg['cluster'][1]:>8.2f} | "
              f"{dw:>+7.2f} {dm:>+7.2f}   {'Y' if dw > 0 else '.'}")
    if do_holdout and len(seeds) > 1:
        print(f"\ncluster beats single — worst-case in {win_worst}/{len(seeds)} seeds, mean in {win_mean}/{len(seeds)} seeds")
    print("\n(end-goal claim: cluster improves WORST-CASE — and ideally mean — over single)")


if __name__ == "__main__":
    main()
