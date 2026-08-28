#!/usr/bin/env python3
"""Budget-matched black-box optimization baselines for the SVDD/AntiFold/BA-DDG headline arm
(reviewer concern #5).

Fix the SAME base (AntiFold antigen-conditioned per-CDR-H3-position marginals, loaded from the
precomputed npz), the SAME reward (BA-DDG binding ddG vs native, CVaR-alpha=0.2 aggregation over
the escape-variant cluster), the SAME held-out omega eval split, and the SAME eval budget
(number of reward-function evaluations), and compare SVDD value-guided decoding against standard
black-box optimizers: best-of-N, random-search, greedy coordinate ascent, beam search, GA.

One invocation = one (target, rep). It:
  1. reproduces SVDD (base=antifold, oracle=baddg_ddg, agg=cvar, k=6, n_designs=6) with an
     instrumented reward that COUNTS reward calls -> B = SVDD's actual reward-eval budget for
     this (target, rep). Writes the SVDD-repro design JSON too (self-consistency vs the screen).
  2. runs every requested optimizer at MATCHED budget B on the SAME base + SAME reward, keeps
     the best `--keep` distinct CDRs by the CVaR reward as that arm's design pool.
  3. writes one design JSON per arm in the identical screen schema so it grades through the SAME
     FoldX / H3-DDG path, and a per-arm compute record (candidates, reward_evals, wall-clock).

Schema + reward + base are reused from the production code (design_for_validation.build_base_dist,
cluster_steer.build_cluster, baddg_oracle.ddg) so the comparison is controlled.
"""
from __future__ import annotations
import argparse, json, math, os, sys, time, csv, importlib.util
from types import SimpleNamespace
import numpy as np
import torch

LYRA = f"{VAAD_ROOT}/projects/mutation_sampling/lyra"
MPNN_DIR = f"{LYRA}/vendor/ProteinMPNN"
sys.path.insert(0, MPNN_DIR)
sys.path.insert(0, f"{LYRA}/src")
sys.path.insert(0, f"{LYRA}/scripts/nos_diffusion")
sys.path.insert(0, f"{LYRA}/scripts/surrogate")

from cluster_steer import (AA, N_AA, build_cluster, cdr_h3_positions, parse_chains,
                           parse_antigen_residues, to_seq)  # noqa: E402
from design_for_validation import build_base_dist, full_residue_map  # noqa: E402
from lyra_mutants.Model.mpnn_adapter import ProteinMPNNAdapter  # noqa: E402
from lyra_mutants.Model.diffusion import sample_svdd  # noqa: E402
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

N_TIMESTEPS = 20   # CDRMaskedDiffusion default; the screen constructs it with this default


def _load_baddg(device):
    p = f"{LYRA}/scripts/surrogate/baddg_oracle.py"
    spec = importlib.util.spec_from_file_location("baddg_oracle", p)
    bo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bo)
    return bo, bo.load_baddg_mpnn(device)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--rep", type=int, required=True)
    ap.add_argument("--npz", required=True, help="AntiFold logits npz (frozen base)")
    ap.add_argument("--clustered-csv", required=True)
    ap.add_argument("--aacdb-dir", default=f"{VAAD_ROOT}/datasets/AACDB/data_zip/complex_structure")
    ap.add_argument("--features-dir", default=f"{VAAD_ROOT}/scratch/complex_features")
    ap.add_argument("--outdir", default=f"{VAAD_ROOT}/scratch/budget_baselines/validation")
    ap.add_argument("--compute-dir", default=f"{VAAD_ROOT}/scratch/budget_baselines/compute")
    ap.add_argument("--max-cluster", type=int, default=6)
    ap.add_argument("--n-designs", type=int, default=6)   # SVDD reproduction n_designs
    ap.add_argument("--svdd-k", type=int, default=6)
    ap.add_argument("--cvar-alpha", type=float, default=0.2)
    ap.add_argument("--objective", default="cvar", choices=["cvar", "mean", "wt"],
                    help="reward aggregation over the escape cluster: cvar=CVaR-alpha tail (default, "
                         "preserves prior behaviour), mean=omega-flat mean over steer variants, "
                         "wt=BA-DDG on the WT antigen only (no variant awareness).")
    ap.add_argument("--eval-cluster", type=int, default=24)
    ap.add_argument("--keep", type=int, default=4, help="designs kept per arm (screen keeps top3+single=4)")
    ap.add_argument("--arms", default="svdd,bestofn,random,greedy,beam,ga")
    ap.add_argument("--bestofn-temp", type=float, default=1.0,
                    help="temperature on the AntiFold base marginals for best-of-N sampling "
                         "(p^(1/T) renorm; T=1.0 = the default arm). Temp sweep uses 0.1/0.5/1.0.")
    ap.add_argument("--budget-B", type=int, default=None,
                    help="skip the SVDD budget-defining run and use this B (temp sweep / re-runs).")
    ap.add_argument("--arm-suffix", default="",
                    help="appended to output arm label so variants (e.g. temp sweep) don't clobber.")
    ap.add_argument("--anytime", action="store_true",
                    help="log best-so-far incumbent at budget fractions -> anytime curve data.")
    ap.add_argument("--anytime-dir",
                    default=f"{VAAD_ROOT}/scratch/budget_baselines/anytime")
    args = ap.parse_args()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.compute_dir, exist_ok=True)
    t = args.target
    seed = 100 * args.rep + 7
    torch.manual_seed(seed)
    np.random.seed(seed)

    # ---------- setup: adapter, base features, native CDR-H3, cluster, held-out eval ----------
    feats0 = torch.load(f"{args.features_dir}/{t}.pt", map_location="cpu", weights_only=False)
    n_ab_chains = len(parse_chains(t))
    cdr_positions = cdr_h3_positions(feats0, n_ab_chains)
    n_cdr = len(cdr_positions)
    adapter = ProteinMPNNAdapter(MPNN_DIR, f"{MPNN_DIR}/vanilla_model_weights/v_48_002.pt",
                                 args.features_dir, cdr_positions).to(device)
    for p in adapter.model.parameters():
        p.requires_grad_(False)
    base = adapter._features(t)
    if base is None:
        raise SystemExit(f"no features for {t}")
    native = base["S"][0, cdr_positions].long().clamp(max=N_AA - 1)
    native_cdr = to_seq(native)
    n_ab = n_ab_chains
    pdb_path = f"{args.aacdb_dir}/{t}.pdb"

    rows = [r for r in csv.DictReader(open(args.clustered_csv, newline="")) if r["complex_id"] == t]
    cluster, ag_len, aligned, _c, labels = build_cluster(
        base, t, rows, args.max_cluster, pdb_path, device, cdr_positions=cdr_positions)
    if not cluster:
        raise SystemExit("empty cluster")

    # held-out omega-divergent eval panel (never used by any optimizer)
    steer_names = {name for name, _ in labels}
    ec, _, _, _, el = build_cluster(base, t, rows, args.eval_cluster + len(steer_names),
                                    pdb_path, device, pick="bottom")
    eval_cluster, eval_labels = [], []
    for (name, toks), member in zip(el, ec):
        if name in steer_names:
            continue
        eval_labels.append((name, toks)); eval_cluster.append(member)
        if len(eval_cluster) >= args.eval_cluster:
            break
    print(f"[{t} r{args.rep}] n_cdr={n_cdr} native={native_cdr} cluster={len(cluster)} "
          f"eval={len(eval_cluster)} aligned={aligned}", flush=True)

    # ---------- CDR residue map + antigen variants (exact screen schema) ----------
    antibody = parse_chains(t)
    chains_str = t.split("_")[-1]
    antigen_chains = [c for c in chains_str if c not in antibody]
    ordered = full_residue_map(pdb_path, antibody + antigen_chains)
    decoded = to_seq(base["S"][0])
    cdr_ok = (len(ordered) == base["S"].shape[1]
              and all(ordered[p][3] == decoded[p] for p in cdr_positions))
    cdr_residues = [{"chain": ordered[p][0], "resnum": ordered[p][1], "icode": ordered[p][2]}
                    for p in cdr_positions]

    def _variant(name, toks, omega, split):
        return {"name": name, "omega": round(omega, 5), "split": split,
                "mutations": [{"chain": c, "wt": wt, "resnum": rn, "icode": ic, "mut": mut}
                              for (c, wt, rn, ic, mut) in toks]}
    antigen_variants = [{"name": "WT", "mutations": [], "split": "both"}]
    for (name, toks), (_, omega) in zip(labels, cluster):
        antigen_variants.append(_variant(name, toks, omega, "steer"))
    for (name, toks), (_, omega) in zip(eval_labels, eval_cluster):
        antigen_variants.append(_variant(name, toks, omega, "eval"))

    # ---------- BA-DDG reward: CVaR-alpha over the escape cluster of (-ddg vs native) ----------
    bo, baddg = _load_baddg(device)

    def _cvar(vals):
        v = torch.tensor(vals, dtype=torch.float)
        k = max(1, math.ceil(args.cvar_alpha * v.numel()))
        return float(torch.sort(v).values[:k].mean())   # mean of the worst tail (higher=better)

    _traj = []          # (reward_value, toks_tuple) in evaluation order, per arm (anytime)
    _traj_on = [False]

    def reward(toks):
        if args.objective == "wt":
            r = -float(bo.ddg(baddg, base, cdr_positions, toks, native, n_ab))
        elif args.objective == "mean":
            svec = [-float(bo.ddg(baddg, f, cdr_positions, toks, native, n_ab)) for f, _ in cluster]
            r = sum(svec) / len(svec)
        else:
            r = _cvar([-float(bo.ddg(baddg, f, cdr_positions, toks, native, n_ab)) for f, _ in cluster])
        if _traj_on[0]:
            _traj.append((r, tuple(int(v) for v in toks)))
        return r

    def heldout_worst(toks):
        """worst-case (max) BA-DDG ddg over the HELD-OUT eval variants (design vs native)."""
        vals = [float(bo.ddg(baddg, f, cdr_positions, toks, native, n_ab)) for f, _ in eval_cluster]
        return max(vals) if vals else None

    def heldout_cvar20(toks):
        """held-out BA-DDG CVaR20 in the SAME sign as the reward (higher=better); for anytime."""
        return _cvar([-float(bo.ddg(baddg, f, cdr_positions, toks, native, n_ab))
                      for f, _ in eval_cluster])

    def score_across(toks):
        s = [-float(bo.ddg(baddg, f, cdr_positions, toks, native, n_ab)) for f, _ in cluster]
        return sum(s) / len(s), min(s)

    def wt_ddg(toks):
        return round(float(bo.ddg(baddg, base, cdr_positions, toks, native, n_ab)), 3)

    def baddg_grade(toks):
        g = {"WT": round(float(bo.ddg(baddg, base, cdr_positions, toks, native, n_ab)), 3)}
        for (name, _tk), (f, _o) in list(zip(labels, cluster)) + list(zip(eval_labels, eval_cluster)):
            g[name] = round(float(bo.ddg(baddg, f, cdr_positions, toks, native, n_ab)), 3)
        return g

    # frozen AntiFold base per-position distribution p(x_i | context) over the 20 AAs
    # base_logits/base_server/live_mpnn_weights were added to build_base_dist by the genzoo
    # campaign; supply them explicitly so this namespace stays compatible with the shared fn.
    bd_args = SimpleNamespace(base="antifold", antifold_logits=args.npz,
                              base_logits=args.npz, base_server=None,
                              live_mpnn_weights=None, svdd_base_temp=1.0)
    base_dist_fn = build_base_dist(bd_args, None, None, adapter, base, native, cdr_positions, device)
    base_probs = base_dist_fn(None, None).clone()               # [n_cdr, 20]
    af_argmax = base_probs.argmax(dim=-1).to(torch.long)        # AntiFold argmax CDR
    # temperature on the base marginals for best-of-N sampling: p^(1/T) renormalized.
    # T=1.0 -> identity (the default best-of-N arm). Lower T -> sharper (nearer argmax).
    if abs(args.bestofn_temp - 1.0) > 1e-9:
        tp = base_probs.clamp_min(1e-12).pow(1.0 / args.bestofn_temp)
        sample_probs = (tp / tp.sum(dim=-1, keepdim=True))
    else:
        sample_probs = base_probs

    def _tok(x):
        return torch.as_tensor(x, device=device, dtype=torch.long)

    def _key(toks):
        return tuple(int(v) for v in toks)

    # ---------- arm runners: each returns (list_of_(tokens,reward), n_reward_evals, wall) ----------
    gen = torch.Generator(device="cpu"); gen.manual_seed(seed)

    def run_svdd(nd):
        """Reproduce the screen's SVDD-PM run; count reward calls -> budget B. Keeps nd finals."""
        calls = [0]
        def rwd(toks):
            calls[0] += 1
            return reward(toks)
        t0 = time.time()
        finals = []
        for _i in range(nd):
            # generator=None: reproduce the screen's UNSEEDED SVDD run; the global torch.manual_seed
            # (set in main) still makes it deterministic per (target,rep) without a device mismatch.
            toks = sample_svdd(base_dist_fn, rwd, n_cdr, N_TIMESTEPS, device,
                               svdd_k=args.svdd_k, svdd_temp=1.0, generator=None)
            finals.append((toks.detach().cpu(), reward(toks)))   # ranking score not counted in B
        return finals, calls[0], time.time() - t0

    def _sample_af(g):
        """One CDR drawn i.i.d. from the AntiFold base per-position marginals (temp-adjusted
        for best-of-N via sample_probs; identical to base_probs when bestofn_temp==1.0)."""
        out = torch.empty(n_cdr, dtype=torch.long)
        for i in range(n_cdr):
            out[i] = torch.multinomial(sample_probs[i].cpu(), 1, generator=g)
        return out.to(device)

    def run_bestofn(B):
        t0 = time.time(); scored = []
        for _ in range(B):
            tk = _sample_af(gen); scored.append((tk.cpu(), reward(tk)))
        return scored, B, time.time() - t0

    def run_random(B):
        t0 = time.time(); scored = []
        for _ in range(B):
            tk = torch.randint(0, N_AA, (n_cdr,), generator=gen).to(device)
            scored.append((tk.cpu(), reward(tk)))
        return scored, B, time.time() - t0

    def run_greedy(B):
        """Greedy steepest-ascent coordinate search with restarts until budget B is spent.
        Restart 0 = AntiFold argmax; further restarts = AntiFold samples. Collect every
        distinct full CDR scored so the pool = best-of-everything-evaluated at budget B."""
        t0 = time.time(); scored = {}; used = [0]
        def ev(tk):
            k = _key(tk)
            if k not in scored:
                scored[k] = (tk.cpu(), reward(tk)); used[0] += 1
            return scored[k][1]
        restart = 0
        while used[0] < B:
            cur = af_argmax.clone() if restart == 0 else _sample_af(gen)
            cur_r = ev(cur)
            improved = True
            while improved and used[0] < B:
                improved = False; best_move = None; best_r = cur_r
                for i in range(n_cdr):
                    if used[0] >= B:
                        break
                    orig = int(cur[i])
                    for aa in range(N_AA):
                        if aa == orig or used[0] >= B:
                            continue
                        cur[i] = aa; r = ev(cur)
                        if r > best_r + 1e-9:
                            best_r = r; best_move = (i, aa)
                    cur[i] = orig
                if best_move is not None:
                    cur[best_move[0]] = best_move[1]; cur_r = best_r; improved = True
            restart += 1
        return list(scored.values()), used[0], time.time() - t0

    def run_beam(B):
        """Beam search over positions (fixed order). Un-decided positions filled by AntiFold
        argmax for full-CDR evaluation. Width chosen so ~ n_cdr*W*N_AA == B."""
        t0 = time.time(); scored = {}; used = [0]
        def ev(tk):
            k = _key(tk)
            if k not in scored:
                scored[k] = (tk.cpu(), reward(tk)); used[0] += 1
            return scored[k][1]
        W = max(2, round(B / (n_cdr * N_AA)))
        beam = [af_argmax.clone()]
        for pos in range(n_cdr):
            cand = []
            for st in beam:
                for aa in range(N_AA):
                    nx = st.clone(); nx[pos] = aa
                    cand.append((nx, ev(nx)))
                if used[0] >= B:
                    break
            cand.sort(key=lambda c: -c[1])
            beam = [c[0] for c in cand[:W]]
            if used[0] >= B:
                break
        return list(scored.values()), used[0], time.time() - t0

    def run_ga(B):
        """Genetic algorithm at matched budget. Init = AntiFold samples; uniform crossover;
        per-position mutation resamples from the AntiFold marginal (stays near the base manifold)."""
        t0 = time.time(); scored = {}; used = [0]
        def ev(tk):
            k = _key(tk)
            if k not in scored:
                scored[k] = (tk.cpu(), reward(tk)); used[0] += 1
            return scored[k][1]
        P = 40; elite = 8; p_mut = 1.5 / n_cdr
        pop = [(_sample_af(gen)) for _ in range(P)]
        pr = [ev(x) for x in pop]
        while used[0] < B:
            order = sorted(range(len(pop)), key=lambda j: -pr[j])
            elites = [pop[j].clone() for j in order[:elite]]
            elite_r = [pr[j] for j in order[:elite]]
            newpop, newr = list(elites), list(elite_r)
            while len(newpop) < P and used[0] < B:
                a = elites[int(torch.randint(0, elite, (1,), generator=gen))]
                b = elites[int(torch.randint(0, elite, (1,), generator=gen))]
                child = a.clone()
                mask = torch.rand(n_cdr, generator=gen) < 0.5
                child[mask.to(device)] = b[mask.to(device)]
                for i in range(n_cdr):
                    if float(torch.rand(1, generator=gen)) < p_mut:
                        child[i] = torch.multinomial(base_probs[i].cpu(), 1, generator=gen).to(device)
                newpop.append(child); newr.append(ev(child))
            pop, pr = newpop, newr
        return list(scored.values()), used[0], time.time() - t0

    RUN = {"bestofn": run_bestofn, "random": run_random, "greedy": run_greedy,
           "beam": run_beam, "ga": run_ga}
    arms = args.arms.split(",")

    FRACTIONS = [0.1, 0.25, 0.5, 0.75, 1.0]

    def _anytime(arm, traj, B):
        recs, designs = [], []
        if not traj:
            return recs, designs
        for fr in FRACTIONS:
            k = max(1, math.ceil(fr * B))
            window = traj[:k]
            bi = max(range(len(window)), key=lambda j: window[j][0])
            bval, btoks = window[bi]
            tk = _tok(list(btoks))
            hw = heldout_worst(tk); hc = heldout_cvar20(tk)
            recs.append({"frac": fr, "n_evals": len(window), "best_steer_reward": round(bval, 4),
                         "heldout_worst_baddg": round(hw, 4) if hw is not None else None,
                         "heldout_cvar20_baddg": round(hc, 4)})
            designs.append({"name": f"{arm}_f{int(fr * 100)}", "cdr_seq": to_seq(tk),
                            "heldout_worst_baddg": round(hw, 4) if hw is not None else None})
        return recs, designs

    suffix = args.arm_suffix
    anytime_all = {}
    anytime_designs = {}

    # ---------- SVDD first (defines budget B) unless a fixed budget is supplied ----------
    B = args.budget_B; compute = {}
    if "svdd" in arms:
        _traj.clear(); _traj_on[0] = args.anytime
        finals, Bsv, wall = run_svdd(args.n_designs)
        _traj_on[0] = False
        if B is None:
            B = Bsv
        _write_arm("svdd", finals, args, t, native_cdr, antibody, antigen_chains, cdr_positions,
                   cdr_residues, antigen_variants, aligned and cdr_ok, score_across, wt_ddg,
                   baddg_grade, keep=args.keep, suffix=suffix)
        compute["svdd"] = {"reward_evals": Bsv, "candidates": args.n_designs, "wall_s": round(wall, 1),
                           "note": f"n_designs={args.n_designs}, svdd_k={args.svdd_k}, "
                                   f"n_reveal_steps~{n_cdr}"}
        if args.anytime:
            anytime_all["svdd"], anytime_designs["svdd"] = _anytime("svdd", list(_traj), Bsv)
        print(f"[{t} r{args.rep}] SVDD budget B={Bsv} reward-evals over {args.n_designs} designs "
              f"({wall:.0f}s)", flush=True)
    if B is None:
        raise SystemExit("need budget B: include 'svdd' in --arms or pass --budget-B")

    for arm in arms:
        if arm == "svdd":
            continue
        _traj.clear(); _traj_on[0] = args.anytime
        scored, nev, wall = RUN[arm](B)
        _traj_on[0] = False
        _write_arm(arm, scored, args, t, native_cdr, antibody, antigen_chains, cdr_positions,
                   cdr_residues, antigen_variants, aligned and cdr_ok, score_across, wt_ddg,
                   baddg_grade, keep=args.keep, suffix=suffix)
        compute[arm] = {"reward_evals": nev, "candidates": len(scored), "wall_s": round(wall, 1)}
        if args.anytime:
            anytime_all[arm], anytime_designs[arm] = _anytime(arm, list(_traj), B)
        print(f"[{t} r{args.rep}] {arm}: {nev} reward-evals, {len(scored)} distinct cands ({wall:.0f}s)",
              flush=True)

    if args.anytime:
        os.makedirs(args.anytime_dir, exist_ok=True)
        json.dump({"target": t, "rep": args.rep, "budget_B": B, "fractions": FRACTIONS,
                   "anytime": anytime_all},
                  open(f"{args.anytime_dir}/{t}_r{args.rep}_anytime.json", "w"), indent=1)
        # H3-gradeable incumbent design JSONs per arm (designs = incumbents at each fraction)
        for arm, designs in anytime_designs.items():
            out = {"target": t, "pdb": f"{t}.pdb", "antibody_chains": antibody,
                   "antigen_chains": antigen_chains, "cdr_positions_flat": list(cdr_positions),
                   "cdr_residues": cdr_residues, "native_cdr": native_cdr,
                   "antibody_designs": designs, "antigen_variants": antigen_variants,
                   "mechanism": f"{arm}_anytime", "base": "antifold", "oracle": "baddg_ddg",
                   "holdout_by": "omega", "alignment_ok": bool(aligned and cdr_ok)}
            json.dump(out, open(f"{args.anytime_dir}/{t}_designs_bb_{arm}_anytime_r{args.rep}.json",
                                "w"), indent=1)

    # suffix/anytime variant runs write their own compute file so they never clobber the
    # canonical per-(target,rep) compute record used by compute_efficiency.csv.
    ctag = suffix if suffix else ("_anytime" if args.anytime else "")
    cpath = f"{args.compute_dir}/{t}{ctag}_r{args.rep}.json"
    json.dump({"target": t, "rep": args.rep, "budget_B": B, "n_cdr": n_cdr,
               "cluster": len(cluster), "eval": len(eval_cluster), "compute": compute},
              open(cpath, "w"), indent=2)
    print(f"wrote compute {cpath}", flush=True)


def _write_arm(arm, scored, args, t, native_cdr, antibody, antigen_chains, cdr_positions,
               cdr_residues, antigen_variants, alignment_ok, score_across, wt_ddg, baddg_grade,
               keep, suffix=""):
    """Keep the top-`keep` DISTINCT CDRs by reward as this arm's design pool; write screen-schema JSON."""
    # dedup by sequence, best reward first
    best = {}
    for tk, r in scored:
        k = tuple(int(v) for v in tk)
        if k not in best or r > best[k][1]:
            best[k] = (tk, r)
    ranked = sorted(best.values(), key=lambda x: -x[1])[:keep]
    designs = []
    for i, (tk, r) in enumerate(ranked):
        tkd = tk.to("cuda" if torch.cuda.is_available() else "cpu").long()
        m, w = score_across(tkd)
        designs.append({"name": f"{arm}_top{i+1}", "cdr_seq": to_seq(tk),
                        "mpnn_mean": round(m, 3), "mpnn_worst": round(w, 3),
                        "cvar_reward": round(float(r), 4),
                        "wt_ddg": wt_ddg(tkd), "baddg_grade": baddg_grade(tkd)})
    out = {"target": t, "pdb": f"{t}.pdb", "antibody_chains": antibody,
           "antigen_chains": antigen_chains, "cdr_positions_flat": list(cdr_positions),
           "cdr_residues": cdr_residues, "native_cdr": native_cdr,
           "antibody_designs": designs, "antigen_variants": antigen_variants,
           "mechanism": arm, "base": "antifold", "oracle": "baddg_ddg",
           "steer_agg": args.objective, "cvar_alpha": args.cvar_alpha, "wt_weight": 0.0,
           "holdout_by": "omega", "alignment_ok": bool(alignment_ok)}
    p = f"{args.outdir}/{t}_designs_bb_{arm}{suffix}_r{args.rep}.json"
    json.dump(out, open(p, "w"), indent=2)


if __name__ == "__main__":
    main()
