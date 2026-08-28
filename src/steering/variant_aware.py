#!/usr/bin/env python3
"""Stage 1 of physics validation: generate-and-filter cluster-steered CDR designs and
dump everything the (CPU) structure-validation stage needs.

Generates --n-designs cluster-steered CDRs (direct MPNN, the same machinery as
cluster_steer.py), scores each across the antigen-mutant cluster by ground-truth MPNN
log-ratio, keeps the top --top-k by WORST-CASE (robustness), plus one single-target
design for contrast. Writes a JSON describing:
  * the CDR positions as PDB (chain, resnum, icode)  -- so ChimeraX can mutate them,
  * the native CDR + each design's CDR sequence + its MPNN mean/worst,
  * the antigen variants (WT + each cluster member's point mutations + omega).

Stage 2 (validate_designs_prodigy.py) consumes this JSON, mutates+relaxes+PRODIGY-scores
each (antibody x antigen) pair, and compares designs vs WT on potency and robustness.

Usage
-----
    PYTHONPATH=src python3 scripts/nos_diffusion/design_for_validation.py \
        --target 1A14_HLN --clustered-csv $RESULTS/clustered.csv \
        --aacdb-dir ${VAAD_ROOT}/datasets/AACDB/data_zip/complex_structure \
        --nos-checkpoint $RESULTS/finetune/nos_4layer/finetune_best.pt \
        --features-dir $SCRATCH/complex_features \
        --proteinmpnn-dir $LYRA/vendor/ProteinMPNN \
        --proteinmpnn-weights $LYRA/vendor/ProteinMPNN/vanilla_model_weights/v_48_002.pt \
        --cdr-positions 95 96 97 98 99 100 101 102 103 104 105 106 \
        --max-cluster 24 --n-designs 32 --top-k 3 --guidance-scale 1.0 \
        --out $RESULTS/validation/1A14_HLN_designs.json
"""
from __future__ import annotations
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

import argparse
import csv
import json
import os
import sys

import numpy as np
import torch

# same-directory helpers (run as a script => this dir is on sys.path)
from profile_util import Profiler
from cluster_steer import (
    AA, N_AA, build_cluster, cdr_h3_positions, gt_logratio, parse_antigen_residues,
    parse_chains, to_seq,
)


# Bases supplied as a precomputed per-CDR-position logit table (.npz, the contract in
# antifold_logits_export.py). Any inverse-folding model becomes a usable SVDD base by writing
# one exporter to this format -- no import of the model into the steering process, which is
# what lets each model keep its own (mutually incompatible) venv.
NPZ_BASES = {"antifold", "proteinmpnn", "ligandmpnn", "esmif", "abmpnn_npz", "iglm"}

# LIVE bases are re-queried at EVERY denoising step with the partially-revealed CDR as context,
# so the generator is genuinely inside the guided-denoising loop rather than supplying a frozen
# marginal up front. This is the training-free method proper: the prior re-conditions as the
# loop fills in, and SVDD's reward reweights that step-conditional prior. The npz variants above
# are the same models with a FROZEN table -- kept as the ablation that isolates what
# re-conditioning buys.
LIVE_MPNN_BASES = {"proteinmpnn_live": None, "abmpnn_live": f"{VAAD_ROOT}/tools/abmpnn/abmpnn.pt"}

# Bases whose model cannot be imported alongside the steering process (conflicting venvs) but
# which still re-condition every denoising step: they run as a resident subprocess speaking
# JSON lines (logit_server.py) and answer one request per step.
SERVER_BASES = {"esmif_live"}


def build_base_dist(args, diffusion, struct, adapter, base, native, cdr_positions, device):
    """Return base_dist_fn(tokens[n_cdr], t) -> probs[n_cdr,20], the per-position next-token
    distribution SVDD samples candidate reveals from, for the selected --base generator.

      diffusion : the masked-diffusion model's own decode softmax(decode(encode(...,tokens,t))).
                  Sequence/step-aware (re-decoded each step from the partial CDR context).
      antifold  : AntiFold's antigen-conditioned per-CDR-position marginal (fixed table loaded
                  from --antifold-logits, the antifold_logits_export.py npz). Unusable positions
                  (unmapped / wt-mismatch) fall back to a native one-hot so SVDD keeps native
                  there — same convention as the AntiFold baseline generator. Antigen structure
                  enters through AntiFold; SVDD reweights by the (BA-DDG) reward -> stays on
                  AntiFold's SOTA inverse-folding manifold.
      abmpnn    : AbMPNN (antibody-tuned ProteinMPNN) native-context complex-conditioned
                  per-position marginal (fixed table), computed via the same ProteinMPNN adapter
                  pointed at the AbMPNN weights.
    Both fixed-table bases are position-marginals (context-independent), so SVDD's value
    reweighting supplies all the cross-position coupling; the diffusion base re-conditions each
    step. All three keep candidate reveals ON the base's manifold (we only sample from it)."""
    import torch.nn.functional as F
    n_cdr = len(cdr_positions)

    if args.base == "diffusion":
        @torch.no_grad()
        def fn(tokens, t):
            h = diffusion.encode(struct, tokens, t - 1)
            return F.softmax(diffusion.decode(h), dim=-1)
        return fn

    if args.base in NPZ_BASES:
        npz = args.base_logits or (args.antifold_logits if args.base == "antifold" else None)
        if not npz:
            raise SystemExit(f"--base {args.base} requires --base-logits <npz> (export it with "
                             f"the matching exporter; see NPZ_BASES in this file)")
        d = np.load(npz, allow_pickle=True)
        logits = d["logits"].astype(np.float64)                 # [n_cdr,20], NaN where unusable
        usable = d["usable"].astype(bool)
        npz_native = str(d["native_cdr"])
        run_native = "".join(AA[int(t)] for t in native.tolist())
        if logits.shape[0] != n_cdr:
            raise SystemExit(f"--base-logits n_cdr={logits.shape[0]} != run n_cdr={n_cdr}")
        if npz_native != run_native:
            raise SystemExit(f"--base-logits native_cdr={npz_native!r} != run "
                             f"native={run_native!r} (CDR window/registration mismatch)")
        tab = np.zeros((n_cdr, N_AA), dtype=np.float64)
        for i in range(n_cdr):
            if usable[i] and np.isfinite(logits[i]).all():
                r = logits[i] - logits[i].max()
                tab[i] = np.exp(r) / np.exp(r).sum()            # softmax over 20 AAs
            else:
                tab[i, int(native[i])] = 1.0                    # keep native where unusable
        probs = torch.tensor(tab, device=device, dtype=torch.float)
        print(f"SVDD base={args.base}: usable {int(usable.sum())}/{n_cdr} positions "
              f"(from {os.path.basename(npz)})", flush=True)

        def fn(tokens, t):
            return probs
        return fn

    if args.base in SERVER_BASES:
        import json as _json
        import subprocess as _sp
        if not args.base_server:
            raise SystemExit(f"--base {args.base} requires --base-server '<cmd>' "
                             f"(see logit_server.py)")
        proc = _sp.Popen(args.base_server, shell=True, stdin=_sp.PIPE, stdout=_sp.PIPE,
                         text=True, bufsize=1)
        hello = _json.loads(proc.stdout.readline())
        if not hello.get("ready"):
            raise SystemExit(f"base server failed to start: {hello}")
        if int(hello.get("n_cdr", -1)) != n_cdr:
            raise SystemExit(f"base server n_cdr={hello.get('n_cdr')} != run n_cdr={n_cdr}")
        print(f"SVDD base={args.base} (LIVE via server, re-conditioned every denoising step)",
              flush=True)
        import atexit
        atexit.register(lambda: proc.stdin.write('{"cmd":"quit"}\n') if proc.poll() is None else None)

        fill = native.clone().to(device)

        def fn(tokens, t):
            nonlocal fill
            revealed = tokens != N_AA
            ctx = fill.clone()
            ctx[revealed] = tokens[revealed]
            cdr = "".join(AA[int(i)] for i in ctx.tolist())
            proc.stdin.write(_json.dumps({"cdr": cdr}) + "\n")
            proc.stdin.flush()
            resp = _json.loads(proc.stdout.readline())
            if "error" in resp:
                raise RuntimeError(f"base server: {resp['error']}")
            lg = torch.tensor(resp["logits"], device=device, dtype=torch.float)
            tab = torch.softmax(lg, dim=-1)
            nxt = tab.argmax(dim=-1)
            nxt[revealed] = tokens[revealed]
            fill = nxt
            return tab

        return fn

    if args.base in LIVE_MPNN_BASES:
        from lyra_mutants.Model.mpnn_adapter import ProteinMPNNAdapter
        w = args.live_mpnn_weights or LIVE_MPNN_BASES[args.base] or args.proteinmpnn_weights
        live_adapter = ProteinMPNNAdapter(
            args.proteinmpnn_dir, w, args.features_dir, cdr_positions).to(device)
        for p in live_adapter.model.parameters():
            p.requires_grad_(False)
        idx = torch.as_tensor(cdr_positions, device=device, dtype=torch.long)
        # Unrevealed CDR positions still have to carry SOME identity, because ProteinMPNN
        # conditions the decoded position on every other residue in S. We fill them with the
        # running predicted-clean (argmax) sequence -- the same predicted-clean convention
        # sample_svdd itself uses for its value estimate -- seeded at the native CDR, which is
        # exactly the context the frozen-table variant of this model was built with. So the two
        # arms start from identical information and differ ONLY in re-conditioning.
        fill = native.clone().to(device)

        @torch.no_grad()
        def fn(tokens, t):
            nonlocal fill
            revealed = tokens != N_AA
            S = base["S"].clone()
            ctx = fill.clone()
            ctx[revealed] = tokens[revealed]
            S[0, idx] = ctx
            chain_M = torch.zeros_like(base["chain_M"])
            todo = ~revealed
            if not bool(todo.any()):            # defensive: the sampler stops before this
                todo = torch.ones_like(revealed)
            chain_M[0, idx[todo]] = 1
            randn = torch.randn(1, base["X"].shape[1], device=device)
            logp = live_adapter.model.conditional_probs(
                base["X"], S, base["mask"], chain_M, base["residue_idx"],
                base["chain_encoding_all"], randn, backbone_only=False)[0]
            tab = logp[idx, :N_AA].exp()
            tab = tab / tab.sum(-1, keepdim=True).clamp_min(1e-12)
            nxt = tab.argmax(dim=-1)
            nxt[revealed] = tokens[revealed]
            fill = nxt
            return tab

        print(f"SVDD base={args.base} (LIVE, re-conditioned every denoising step) weights={w}",
              flush=True)
        return fn

    if args.base == "abmpnn":
        from lyra_mutants.Model.mpnn_adapter import ProteinMPNNAdapter
        w = args.abmpnn_weights
        ab_adapter = ProteinMPNNAdapter(
            args.proteinmpnn_dir, w, args.features_dir, cdr_positions).to(device)
        for p in ab_adapter.model.parameters():
            p.requires_grad_(False)
        S = base["S"].clone()
        S[0, cdr_positions] = native
        chain_M = torch.zeros_like(base["chain_M"])
        chain_M[0, cdr_positions] = 1
        randn = torch.randn(1, base["X"].shape[1], device=device)
        with torch.no_grad():
            logp = ab_adapter.model.conditional_probs(
                base["X"], S, base["mask"], chain_M, base["residue_idx"],
                base["chain_encoding_all"], randn, backbone_only=False)[0]
        idx = torch.as_tensor(cdr_positions, device=device, dtype=torch.long)
        tab = logp[idx, :N_AA].exp()
        probs = (tab / tab.sum(-1, keepdim=True)).detach()
        print(f"SVDD base=abmpnn: native-context table from {w}", flush=True)

        def fn(tokens, t):
            return probs
        return fn

    raise SystemExit(f"unknown --base {args.base!r}")


def full_residue_map(pdb_path, ordered_chains):
    """Flat residue list across chains in the given order, mirroring parse_PDB layout
    so index i lines up with feats['S'][0, i]. -> [(chain, resnum|None, icode, aa1)]."""
    ordered = []
    for ch in ordered_chains:
        ordered.extend(parse_antigen_residues(pdb_path, [ch]))
    return ordered


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--clustered-csv", required=True)
    ap.add_argument("--aacdb-dir", required=True)
    ap.add_argument("--nos-checkpoint", required=True)
    ap.add_argument("--features-dir", required=True)
    ap.add_argument("--proteinmpnn-dir", required=True)
    ap.add_argument("--proteinmpnn-weights", required=True)
    ap.add_argument("--cdr-positions", nargs="+", type=int, default=None,
                    help="explicit flat CDR indices (legacy fixed window); omit when using "
                         "--auto-cdr-h3")
    ap.add_argument("--auto-cdr-h3", action="store_true",
                    help="derive the CDR-H3 window per target from the conserved Cys->WGxG "
                         "anchors (full, correctly-registered loop) instead of a fixed slice")
    ap.add_argument("--max-cluster", type=int, default=24)
    ap.add_argument("--n-designs", type=int, default=32, help="cluster CDRs to sample before filtering")
    ap.add_argument("--top-k", type=int, default=3, help="best designs (by worst-case) to keep")
    ap.add_argument("--guidance-scale", type=float, default=1.0)
    ap.add_argument("--guidance-mode", choices=["standard", "guo"], default="standard",
                    help="'guo' = forward-prediction (look-ahead) loss guidance, Guo et al. 2024 "
                         "(arXiv:2404.14743): self-limiting step toward a reachable energy target, "
                         "keeps designs on-manifold vs the constant unit-step push of 'standard'.")
    ap.add_argument("--guo-margin", type=float, default=2.0,
                    help="(guidance-mode guo) per-step reachable energy-improvement target (energy "
                         "units); the look-ahead squared-loss target is E0 - guo_margin.")
    ap.add_argument("--oracle", choices=["mpnn", "baddg_ddg", "consensus", "ptm_energy", "mint"],
                    default="mpnn",
                    help="steering energy: ProteinMPNN complex log-ratio (default) or BA-DDG "
                         "differentiable binding ddG (Stage-B B2), or the Stage-B B3 CONSENSUS "
                         "z(baddg_w)+lam*z(mpnn_w)[+mu*z(h3_w)] energy, or PTM_ENERGY = in-loop "
                         "differentiable Boltz-2 interface pTMEnergy (BindEnergyCraft, "
                         "arXiv:2505.21241; matrix arms A3/A4), or MINT = out-of-ProteinMPNN-family "
                         "MINT (ESM-2 650M + cross-chain attention) binding ddG on realized CDRs "
                         "(SVDD reward only — score_ddg site-head; no differentiable soft-CDR path). "
                         "baddg_if (B1) was lost to the revert and not restored.")
    ap.add_argument("--mint-site-head",
                    default=f"{VAAD_ROOT}/scratch/steer_tests/mint_anchor/mint_site_head.pt",
                    help="(oracle mint) trained MINTSiteDDGHead checkpoint (state_dict+ymean+ystd)")
    ap.add_argument("--ptm-panel-k", type=int, default=4,
                    help="(ptm_energy) TRACTABILITY CAP: fold only the top-k (highest-omega) escape "
                         "variants per guided step with Boltz-2 (default 4). Folding the whole "
                         "cluster every step is the cost driver. A3/A4 therefore aggregate over a "
                         "SUBSET of the panel that A1/A2 (BA-DDG) see in full — a comparison caveat.")
    ap.add_argument("--ptm-guided-steps", type=int, default=4,
                    help="(ptm_energy) number of denoising steps to actually steer, converted to an "
                         "effective --guide-last-frac given the CDR length (the masked-diffusion "
                         "schedule unmasks one position/step, so the active window is the first "
                         "n_cdr high-t steps, NOT the last fraction). 0 = use --guide-last-frac "
                         "verbatim. Bounds Boltz calls to guided_steps x 3 x panel_k per design.")
    ap.add_argument("--lam", type=float, default=1.0,
                    help="(consensus) weight on the ProteinMPNN complex-IF log-ratio term")
    ap.add_argument("--mu", type=float, default=1.0,
                    help="(consensus) intended weight on the H3-DDG term; H3-DDG has no "
                         "differentiable soft-CDR forward, so the h3 term is NOT applied during "
                         "steering — mu is recorded in metadata and reported as omitted.")
    ap.add_argument("--guide-last-frac", type=float, default=1.0,
                    help="(baddg_ddg) steer only the last fraction of denoising steps")
    ap.add_argument("--prior-weight", type=float, default=0.0,
                    help="(baddg_ddg) KL-to-prior manifold anchor strength (step 3); 0 = off")
    ap.add_argument("--steer-agg", choices=["mean", "worst", "cvar", "fused"], default="mean",
                    help="(baddg_ddg) how to combine the escape panel during steering: "
                         "omega-weighted mean (current), single worst-case, CVaR worst-tail, or "
                         "'fused' = escape-likelihood-weighted CVaR (tolerability x binding, "
                         "P(escape) weights from --pesc-json; EXP-D generation objective)")
    ap.add_argument("--pesc-json", default=None,
                    help="(steer-agg=fused) JSON {steer_variant_name: p_escape_weight} used as the "
                         "CVaR tail measure; keys are cluster_subsets mutation_strings")
    ap.add_argument("--cvar-alpha", type=float, default=0.2,
                    help="(steer-agg=cvar/fused) worst-tail fraction; α=1->mean, α->0->worst")
    ap.add_argument("--wt-weight", type=float, default=0.0,
                    help="(baddg_ddg) weight on a WT-antigen binding term added to the steering "
                         "energy, so robustness steering keeps WT potency. 0=pure robustness")
    ap.add_argument("--refold-dir", default=None,
                    help="(baddg_ddg|consensus) dir of per-variant re-folded PDBs "
                         "({safe(mutation_string)}.pdb from the refold panel stage2). When set, each "
                         "escape cluster member's backbone X is replaced by its OWN re-folded "
                         "structure — fixing the frozen-backbone degeneracy (frozen-backbone-gradient-"
                         "rootcause) where escape variants differed only in antigen sequence, not "
                         "geometry. Also disables the antibody-alone term cache (per-variant antibody "
                         "coords differ). Variants without a refold fall back to the frozen backbone.")
    ap.add_argument("--antifold-logits", default=None,
                    help="(baddg_ddg|consensus) .npz from antifold_logits_export.py: AntiFold's "
                         "antigen-conditioned per-CDR-position logits, row-aligned to this target's "
                         "CDR window. Loaded as the product-of-experts prior; requires matching "
                         "native_cdr (same target + CDR registration).")
    ap.add_argument("--antifold-prior-weight", type=float, default=0.0,
                    help="(baddg_ddg|consensus) strength of the AntiFold PoE prior term "
                         "+w·(-Σ probs·logP_af) added to the steering energy; 0 = off. Pulls the "
                         "design toward AntiFold's SOTA distribution (the base/floor) while BA-DDG "
                         "supplies the binding/robustness signal AntiFold is blind to.")
    ap.add_argument("--warmstart", action="store_true",
                    help="(baddg_ddg|consensus, needs --antifold-logits) Fix B: freeze a scaffold of "
                         "AntiFold's argmax residues and re-optimize only the OPEN positions with "
                         "BA-DDG. Hard manifold scaffold (real AntiFold residues), not a soft prior.")
    ap.add_argument("--warmstart-open-frac", type=float, default=0.5,
                    help="(warmstart) fraction of CDR positions left OPEN for BA-DDG steering; the "
                         "rest are frozen at AntiFold argmax. Unusable-by-AntiFold positions are always "
                         "open. 0 = pure AntiFold (sanity: reproduces the baseline), 1 = pure STEER.")
    ap.add_argument("--holdout-by", choices=["none", "random", "omega"], default="none",
                    help="held-out generalization test: steer on the top-omega panel, then ALSO "
                         "FoldX-score the design on a held-out eval panel never steered on. "
                         "omega=divergent (low-omega) eval = the decisive distribution-shift test; "
                         "random=i.i.d. eval = a weaker floor.")
    ap.add_argument("--eval-cluster", type=int, default=24,
                    help="(holdout) size of the held-out eval panel")
    ap.add_argument("--holdout-seed", type=int, default=0,
                    help="(holdout-by random) rng seed for the eval sample")
    ap.add_argument("--mechanism", choices=["gradient", "svdd", "smc"], default="gradient",
                    help="steering MECHANISM. 'gradient' (default) = the existing NOS gradient "
                         "guidance: descend a DIFFERENTIABLE energy_fn(soft_probs) on the hidden "
                         "states. 'svdd' = SVDD-PM value-based decoding (arXiv:2408.08252): at each "
                         "unmask step draw --svdd-k candidate reveals from the base model's own "
                         "distribution, score each on its REALIZED discrete predicted-clean CDR with "
                         "the reward, and resample ∝ exp(reward/--svdd-temp). Derivative-free, so it "
                         "works with NON-differentiable rewards and on any --base generator.")
    ap.add_argument("--base", choices=["diffusion", "antifold", "abmpnn",
                                       "proteinmpnn", "ligandmpnn", "esmif", "abmpnn_npz", "iglm",
                                       "proteinmpnn_live", "abmpnn_live", "esmif_live"],
                    default="diffusion",
                    help="(svdd) BASE generator whose per-position distribution SVDD samples "
                         "candidate reveals from. 'diffusion' = the masked-diffusion model (default); "
                         "'antifold' = AntiFold antigen-conditioned marginals (needs --antifold-logits); "
                         "'abmpnn' = AbMPNN native-context marginals (--abmpnn-weights). "
                         "IGNORED by --mechanism gradient (always the diffusion model).")
    ap.add_argument("--svdd-k", type=int, default=12,
                    help="(svdd) number of candidate next-states drawn from the base per step")
    ap.add_argument("--svdd-temp", type=float, default=1.0,
                    help="(svdd) resampling temperature; ->0 greedy on reward, large -> base sampling")
    ap.add_argument("--n-particles", type=int, default=12,
                    help="(smc) number of SMC particles (twisted-SMC / FK-steering; the multi-particle "
                         "generalisation of svdd). Uses the same reward/base_dist wiring as svdd.")
    ap.add_argument("--smc-temp", type=float, default=1.0, help="(smc) twist temperature")
    ap.add_argument("--n-timesteps", type=int, default=0,
                    help="override the diffusion n_timesteps (0 = use the model default, 20). "
                         "Sweepable robustness knob (F); the checkpoint does not constrain it.")
    ap.add_argument("--base-logits", default=None,
                    help="npz logit table for any NPZ_BASES base (proteinmpnn/ligandmpnn/esmif/"
                         "iglm/antifold). Same contract as antifold_logits_export.py.")
    ap.add_argument("--svdd-base-temp", type=float, default=1.0,
                    help="sharpen the base distribution before drawing SVDD candidates "
                         "(1.0 = raw prior; 0.2 matches the baselines' sampling temperature)")
    ap.add_argument("--trajectory-out", default=None,
                    help="write the full per-step inference trajectory (prior proposals, reward "
                         "of every candidate, resample weights, which survived) as JSON")
    ap.add_argument("--base-server", default=None,
                    help="(--base *_live served) shell command starting logit_server.py in the "
                         "model's own venv; spoken to over stdin/stdout, one request per step")
    ap.add_argument("--live-mpnn-weights", default=None,
                    help="(--base *_live) ProteinMPNN-format checkpoint re-queried each "
                         "denoising step; defaults per LIVE_MPNN_BASES")
    ap.add_argument("--abmpnn-weights", default=f"{VAAD_ROOT}/tools/abmpnn/abmpnn.pt",
                    help="(svdd --base abmpnn) AbMPNN ProteinMPNN-format checkpoint")
    ap.add_argument("--n-layers", type=int, default=4)
    ap.add_argument("--hidden-dim", type=int, default=128)
    ap.add_argument("--out", required=True, help="output JSON path")
    args = ap.parse_args()

    # --- compute-cost accounting (Phase 3): every arm is profiled identically so the
    # training-free/sampling-cheap claim can be checked against the baselines' real cost.
    prof = Profiler(tag=f"{args.mechanism}/{args.base}/{args.oracle}/{args.steer_agg}",
                    target=args.target,
                    meta={"n_designs": args.n_designs, "svdd_k": args.svdd_k,
                          "max_cluster": args.max_cluster, "guidance_scale": args.guidance_scale,
                          "cvar_alpha": args.cvar_alpha, "top_k": args.top_k})
    if args.mechanism in ("svdd", "smc") and args.oracle not in ("baddg_ddg", "consensus", "mpnn", "mint"):
        ap.error(f"--mechanism {args.mechanism} currently supports --oracle baddg_ddg|consensus|mpnn|mint "
                 "(the reward is evaluated on realized discrete CDRs)")
    if args.oracle == "mint" and args.mechanism != "svdd":
        ap.error("--oracle mint is realized-sequence only (no differentiable soft-CDR forward); "
                 "use it with --mechanism svdd")
    if args.warmstart:
        if not args.antifold_logits:
            ap.error("--warmstart requires --antifold-logits (the scaffold source)")
        if args.oracle not in ("baddg_ddg", "consensus"):
            ap.error("--warmstart is only implemented for --oracle baddg_ddg|consensus")
        if not (0.0 <= args.warmstart_open_frac <= 1.0):
            ap.error("--warmstart-open-frac must be in [0,1]")
    if args.refold_dir and args.oracle not in ("baddg_ddg", "consensus", "ptm_energy"):
        ap.error("--refold-dir is only implemented for --oracle baddg_ddg|consensus|ptm_energy")

    sys.path.insert(0, args.proteinmpnn_dir)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    from lyra_mutants.Model.diffusion import CDRMaskedDiffusion
    from lyra_mutants.Model.mpnn_adapter import ProteinMPNNAdapter

    # derive the CDR-H3 window per target (correct registration) before sizing the model
    if args.auto_cdr_h3:
        feats0 = torch.load(f"{args.features_dir}/{args.target}.pt",
                            map_location="cpu", weights_only=False)
        args.cdr_positions = cdr_h3_positions(feats0, len(parse_chains(args.target)))
        print(f"auto CDR-H3 [{args.target}]: flat {args.cdr_positions[0]}..{args.cdr_positions[-1]} "
              f"({len(args.cdr_positions)} residues) = "
              f"{to_seq(feats0['S'][0, args.cdr_positions])}", flush=True)
    elif not args.cdr_positions:
        ap.error("provide --cdr-positions or --auto-cdr-h3")

    n_cdr = len(args.cdr_positions)
    _ndt = {"n_timesteps": args.n_timesteps} if args.n_timesteps > 0 else {}
    diffusion = CDRMaskedDiffusion(
        n_cdr=n_cdr, hidden_dim=args.hidden_dim, n_layers=args.n_layers,
        binder_dim=128, scorer_dim=128, **_ndt).to(device)
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
    cluster, ag_len, aligned, (n_parse, n_map, n_wt), labels = build_cluster(
        base, args.target, rows, args.max_cluster, pdb_path, device,
        refold_dir=args.refold_dir, cdr_positions=args.cdr_positions)
    print(f"target {args.target}: cluster size {len(cluster)}, antigen len {ag_len}, aligned={aligned} "
          f"(skipped parse {n_parse}, resnum {n_map}, wt {n_wt})")
    if not cluster:
        raise SystemExit("empty cluster — check aacdb-dir / mutation_string / alignment")

    # --- held-out eval panel (never steered on): divergence (low-omega) or random ---
    eval_cluster, eval_labels = [], []
    if args.holdout_by != "none":
        steer_names = {name for name, _ in labels}
        pick = "bottom" if args.holdout_by == "omega" else "random"
        ec, _, _, _, el = build_cluster(
            base, args.target, rows, args.eval_cluster + len(steer_names), pdb_path, device,
            pick=pick, seed=args.holdout_seed)
        # drop any eval member that is also in the steer panel, then cap to --eval-cluster
        for (name, toks), member in zip(el, ec):
            if name in steer_names:
                continue
            eval_labels.append((name, toks))
            eval_cluster.append(member)
            if len(eval_cluster) >= args.eval_cluster:
                break
        print(f"holdout '{args.holdout_by}': steer on {len(cluster)} / FoldX-eval on "
              f"{len(eval_cluster)} HELD-OUT variants", flush=True)

    # --- map CDR flat indices -> PDB (chain, resnum, icode), with an alignment guard ---
    antibody = parse_chains(args.target)
    chains_str = args.target.split("_")[-1]
    antigen_chains = [c for c in chains_str if c not in antibody]
    ordered = full_residue_map(pdb_path, antibody + antigen_chains)
    decoded = to_seq(base["S"][0])
    cdr_ok = (len(ordered) == base["S"].shape[1]
              and all(ordered[p][3] == decoded[p] for p in args.cdr_positions))
    if not cdr_ok:
        print("WARNING: CDR PDB-residue map does not align with feats['S'] — "
              f"parsed {len(ordered)} residues vs S length {base['S'].shape[1]}; "
              "ChimeraX residue targets may be wrong. Inspect before trusting Stage 2.")
    cdr_residues = [{"chain": ordered[p][0], "resnum": ordered[p][1], "icode": ordered[p][2]}
                    for p in args.cdr_positions]
    native_cdr = to_seq(native)

    # --- oracle: scoring (ranking) + sampling, per --oracle ---
    n_ab = len(antibody)
    # per-variant re-folded backbones => antibody-alone coords differ per member => can't share the cache
    share_ab_logp = args.refold_dir is None
    if args.oracle in ("baddg_ddg", "consensus"):
        import importlib.util as _ilu
        _bo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "surrogate", "baddg_oracle.py")
        _spec = _ilu.spec_from_file_location("baddg_oracle", _bo_path)
        bo = _ilu.module_from_spec(_spec)
        _spec.loader.exec_module(bo)
        baddg = bo.load_baddg_mpnn(device)

        def score_across(toks):
            s = [-bo.ddg(baddg, f, args.cdr_positions, toks, native, n_ab) for f, _ in cluster]
            return sum(s) / len(s), min(s)

        def wt_ddg(toks):
            # design-vs-native binding ddG on the WT antigen (base complex); <0 = tighter than native
            return float(bo.ddg(baddg, base, args.cdr_positions, toks, native, n_ab))

        # AntiFold product-of-experts prior table (optional): [n_cdr,20] log-probs, unusable
        # positions zeroed, aligned to this run's CDR window. Guarded on native_cdr so a stale
        # npz (wrong target/registration) fails loudly instead of steering toward the wrong loop.
        _af_table = None
        _ws_init = _ws_frozen = None                # warm-start scaffold (init_tokens, frozen_mask)
        if (args.antifold_prior_weight or args.warmstart) and args.antifold_logits:
            _d = np.load(args.antifold_logits, allow_pickle=True)
            _af_logits = _d["logits"].astype(np.float64)          # [n_cdr,20], NaN where unusable
            _af_usable = _d["usable"].astype(bool)
            _npz_native = str(_d["native_cdr"])
            _run_native = "".join(AA[int(t)] for t in native.tolist())
            if _af_logits.shape[0] != n_cdr:
                raise SystemExit(f"--base-logits n_cdr={_af_logits.shape[0]} != run n_cdr={n_cdr}")
            if _npz_native != _run_native:
                raise SystemExit(
                    f"--antifold-logits native_cdr={_npz_native!r} != run native={_run_native!r} "
                    "(CDR window/registration mismatch — regenerate the npz for this run)")
            _tab = np.zeros((n_cdr, N_AA), dtype=np.float64)
            for _i in range(n_cdr):
                if _af_usable[_i] and np.isfinite(_af_logits[_i]).all():
                    _r = _af_logits[_i] - _af_logits[_i].max()
                    _tab[_i] = _r - np.log(np.exp(_r).sum())    # log-softmax; unusable rows stay 0
            if args.antifold_prior_weight:
                _af_table = torch.tensor(_tab, device=device, dtype=torch.float)
                print(f"AntiFold prior: weight={args.antifold_prior_weight} "
                      f"usable {int(_af_usable.sum())}/{n_cdr} positions", flush=True)
            if args.warmstart:
                # AntiFold argmax scaffold; freeze the most-confident usable positions, OPEN the rest.
                _argmax = _tab.argmax(axis=1)                       # [n_cdr]
                _init = native.clone()                             # fallback for unusable = native
                for _i in range(n_cdr):
                    if _af_usable[_i]:
                        _init[_i] = int(_argmax[_i])
                _conf = _tab.max(axis=1)                            # log max-prob; higher = more confident
                _usable_idx = [i for i in range(n_cdr) if _af_usable[i]]
                _n_open = max(0, round(args.warmstart_open_frac * n_cdr) - (n_cdr - len(_usable_idx)))
                # open = the least-confident usable positions (+ all unusable, always open)
                _open_usable = sorted(_usable_idx, key=lambda i: _conf[i])[:_n_open]
                _frozen = torch.zeros(n_cdr, dtype=torch.bool)
                for i in _usable_idx:
                    if i not in _open_usable:
                        _frozen[i] = True
                _ws_init = _init.to(device=device, dtype=torch.long)
                _ws_frozen = _frozen.to(device)
                print(f"WARM-START: open_frac={args.warmstart_open_frac} -> frozen "
                      f"{int(_frozen.sum())}/{n_cdr} at AntiFold argmax, open {n_cdr-int(_frozen.sum())} "
                      f"(usable {len(_usable_idx)}/{n_cdr}); scaffold="
                      f"{''.join(AA[int(t)] for t in _ws_init.tolist())}", flush=True)

        if args.oracle == "baddg_ddg":
            print(f"oracle: BA-DDG binding ddG (guide_last_frac={args.guide_last_frac})", flush=True)

            def sample_cluster():
                ef = bo.make_ddg_energy_fn(baddg, cluster, args.cdr_positions, n_ab,
                                           agg=args.steer_agg, cvar_alpha=args.cvar_alpha,
                                           wt_feats=base, wt_weight=args.wt_weight,
                                           antifold_table=_af_table,
                                           antifold_weight=args.antifold_prior_weight,
                                           share_ab_logp=share_ab_logp)
                return diffusion.sample_nos_energy(
                    struct, ef, args.cdr_positions, args.guidance_scale,
                    args.guide_last_frac, args.prior_weight,
                    init_tokens=_ws_init, frozen_mask=_ws_frozen,
                    guidance_mode=args.guidance_mode, guo_margin=args.guo_margin)[:n_cdr]

            def sample_single():
                ef = bo.make_ddg_energy_fn(baddg, [(base, 1.0)], args.cdr_positions, n_ab,
                                           agg=args.steer_agg, cvar_alpha=args.cvar_alpha,
                                           wt_feats=base, wt_weight=args.wt_weight,
                                           antifold_table=_af_table,
                                           antifold_weight=args.antifold_prior_weight)
                return diffusion.sample_nos_energy(
                    struct, ef, args.cdr_positions, args.guidance_scale,
                    args.guide_last_frac, args.prior_weight,
                    init_tokens=_ws_init, frozen_mask=_ws_frozen,
                    guidance_mode=args.guidance_mode, guo_margin=args.guo_margin)[:n_cdr]
        else:  # consensus: BA-DDG binding (+) lam * vendor-ProteinMPNN complex-IF log-ratio
            print(f"oracle: CONSENSUS baddg + lam*mpnn (lam={args.lam}, mu={args.mu}, "
                  f"agg={args.steer_agg}, guide_last_frac={args.guide_last_frac})", flush=True)
            if args.mu:
                print("  NOTE: mu>0 but the H3-DDG term is NOT wired for steering (grader-only, "
                      "no differentiable soft-CDR forward); the h3 term is OMITTED.", flush=True)

            @torch.no_grad()
            def _mpnn_table(f):
                """[n_cdr,20] vendor-ProteinMPNN complex log-prob table in f's antigen context at
                the native CDR — the fixed reward table for the consensus mpnn log-ratio term."""
                S = f["S"].clone()
                S[0, args.cdr_positions] = native
                chain_M = torch.zeros_like(f["chain_M"])
                chain_M[0, args.cdr_positions] = 1
                randn = torch.randn(1, f["X"].shape[1], device=device)
                logp = adapter.model.conditional_probs(
                    f["X"], S, f["mask"], chain_M, f["residue_idx"],
                    f["chain_encoding_all"], randn, backbone_only=False)[0]
                idx = torch.as_tensor(args.cdr_positions, device=device, dtype=torch.long)
                return logp[idx, :N_AA].detach()

            _cluster_tables = [_mpnn_table(f) for f, _ in cluster]
            _wt_table = _mpnn_table(base)

            def sample_cluster():
                ef = bo.make_consensus_energy_fn(
                    baddg, cluster, args.cdr_positions, n_ab, _cluster_tables, native,
                    lam=args.lam, agg=args.steer_agg, cvar_alpha=args.cvar_alpha,
                    wt_feats=base, wt_mpnn_table=_wt_table, wt_weight=args.wt_weight,
                    antifold_table=_af_table, antifold_weight=args.antifold_prior_weight,
                    share_ab_logp=share_ab_logp)
                return diffusion.sample_nos_energy(
                    struct, ef, args.cdr_positions, args.guidance_scale,
                    args.guide_last_frac, args.prior_weight,
                    init_tokens=_ws_init, frozen_mask=_ws_frozen,
                    guidance_mode=args.guidance_mode, guo_margin=args.guo_margin)[:n_cdr]

            def sample_single():
                ef = bo.make_consensus_energy_fn(
                    baddg, [(base, 1.0)], args.cdr_positions, n_ab, [_wt_table], native,
                    lam=args.lam, agg=args.steer_agg, cvar_alpha=args.cvar_alpha,
                    wt_feats=base, wt_mpnn_table=_wt_table, wt_weight=args.wt_weight,
                    antifold_table=_af_table, antifold_weight=args.antifold_prior_weight)
                return diffusion.sample_nos_energy(
                    struct, ef, args.cdr_positions, args.guidance_scale,
                    args.guide_last_frac, args.prior_weight,
                    init_tokens=_ws_init, frozen_mask=_ws_frozen,
                    guidance_mode=args.guidance_mode, guo_margin=args.guo_margin)[:n_cdr]
    elif args.oracle == "ptm_energy":
        # ---- A3/A4: in-loop differentiable Boltz-2 interface pTMEnergy (BindEnergyCraft) ----
        import importlib.util as _ilu
        _pe_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "surrogate", "ptm_energy_oracle.py")
        _spec = _ilu.spec_from_file_location("ptm_energy_oracle", _pe_path)
        pe = _ilu.module_from_spec(_spec)
        _spec.loader.exec_module(pe)
        chain_ids = antibody + antigen_chains          # ordered: antibody chains then antigen

        # the masked-diffusion schedule unmasks 1 position/step, so the active (still-masked) window
        # is the first n_cdr high-t steps; --guide-last-frac<~ (21-n_cdr)/n_timesteps disables ALL
        # guidance. Convert --ptm-guided-steps into the frac that steers exactly that many steps.
        n_ts = diffusion.n_timesteps
        if args.ptm_guided_steps > 0:
            ptm_frac = min(1.0, max(1.0 / n_ts, (n_ts - n_cdr + args.ptm_guided_steps) / n_ts))
            print(f"oracle: PTM_ENERGY (Boltz-2), panel_k={args.ptm_panel_k}, agg={args.steer_agg}, "
                  f"guided_steps~{min(args.ptm_guided_steps, n_cdr)} -> guide_last_frac={ptm_frac:.3f} "
                  f"(n_cdr={n_cdr}, n_timesteps={n_ts})", flush=True)
        else:
            ptm_frac = args.guide_last_frac
            print(f"oracle: PTM_ENERGY (Boltz-2), panel_k={args.ptm_panel_k}, agg={args.steer_agg}, "
                  f"guide_last_frac={ptm_frac} (verbatim)", flush=True)
        if args.antifold_prior_weight or args.warmstart:
            ap.error("--antifold-prior-weight/--warmstart are not implemented for --oracle ptm_energy")

        def wt_ddg(toks):
            return None

        # rank the sampled designs with the cheap ProteinMPNN complex log-ratio (out-of-oracle,
        # so the pTM steerer does not also grade its own designs); pTM only steers.
        def score_across(toks):
            s = [gt_logratio(adapter.model, f, args.cdr_positions, toks, native, device) for f, _ in cluster]
            return sum(s) / len(s), min(s)

        def sample_cluster():
            ef = pe.make_ptm_energy_fn(
                args.target, cluster, args.cdr_positions, n_ab, chain_ids,
                agg=args.steer_agg, cvar_alpha=args.cvar_alpha,
                wt_pmpnn_feats=base, wt_weight=args.wt_weight,
                panel_k=args.ptm_panel_k, device=device, labels=labels)
            return diffusion.sample_nos_energy(
                struct, ef, args.cdr_positions, args.guidance_scale,
                ptm_frac, args.prior_weight,
                guidance_mode=args.guidance_mode, guo_margin=args.guo_margin)[:n_cdr]

        def sample_single():
            ef = pe.make_ptm_energy_fn(
                args.target, [(base, 1.0)], args.cdr_positions, n_ab, chain_ids,
                agg=args.steer_agg, cvar_alpha=args.cvar_alpha,
                wt_pmpnn_feats=base, wt_weight=args.wt_weight,
                panel_k=args.ptm_panel_k, device=device, labels=[("WT", None)])
            return diffusion.sample_nos_energy(
                struct, ef, args.cdr_positions, args.guidance_scale,
                ptm_frac, args.prior_weight,
                guidance_mode=args.guidance_mode, guo_margin=args.guo_margin)[:n_cdr]
    elif args.oracle == "mint":
        # ---- MINT (out-of-ProteinMPNN-family, ESM-2 650M + cross-chain attention) binding ddG on
        #      REALIZED CDRs, as an SVDD reward. Antigen-conditioned + sequence-only: each escape
        #      variant is just a SWAPPED ANTIGEN SEQUENCE (no refold). A whole-CDR redesign is scored
        #      as the sum of per-CDR-position site ddGs (mint_oracle.cdr_redesign_ddg), reusing the
        #      validated site-head. mint is SVDD-only (guarded above); the sample_* defined here are
        #      overwritten by the SVDD block below — they only guard against accidental gradient use.
        if args.antifold_prior_weight or args.warmstart:
            ap.error("--antifold-prior-weight/--warmstart are not implemented for --oracle mint")
        import importlib.util as _ilu
        _mo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "surrogate", "mint_oracle.py")
        _spec = _ilu.spec_from_file_location("mint_oracle", _mo_path)
        MO = _ilu.module_from_spec(_spec)
        _spec.loader.exec_module(MO)
        _mint = MO.load_mint(device=device)
        _mint_head, _mint_ymean, _mint_ystd = MO.load_site_head(args.mint_site_head, device)

        enc0 = base["chain_encoding_all"][0]
        _antigen_idx = (enc0 > n_ab).nonzero(as_tuple=True)[0]
        _antibody_idx = (enc0 <= n_ab).nonzero(as_tuple=True)[0]
        _ab_pos = {int(fi): i for i, fi in enumerate(_antibody_idx.tolist())}
        _site_idx = [_ab_pos[int(p)] for p in args.cdr_positions]   # CDR pos -> index within ab chain
        _ab_native = to_seq(base["S"][0, _antibody_idx])
        print(f"oracle: MINT binding ddG (site-head={os.path.basename(args.mint_site_head)}, "
              f"ab_len={len(_ab_native)}, ag_len={len(_antigen_idx)}, agg={args.steer_agg})",
              flush=True)

        def _mint_prep(feats_list):
            """Precompute (antigen_seq, native-antibody site reps) per complex — reused across every
            candidate (only the mutant forward depends on the candidate CDR)."""
            prep = []
            for f, _ in feats_list:
                ag = to_seq(f["S"][0, _antigen_idx])
                wt_reps = MO.embed_cdr_sites(_mint, _ab_native, ag, _site_idx, device)
                prep.append((ag, wt_reps))
            return prep

        def _mint_ddg_one(ag_seq, wt_reps, toks):
            """ddG of substituting the candidate CDR into the antibody, on antigen `ag_seq`.
            Positive = destabilising (weaker binding than native)."""
            ab_cand = list(_ab_native)
            for k, si in enumerate(_site_idx):
                ab_cand[si] = AA[int(toks[k])]
            changed = torch.tensor([int(toks[k]) != int(native[k]) for k in range(n_cdr)],
                                   device=device)
            mut_reps = MO.embed_cdr_sites(_mint, "".join(ab_cand), ag_seq, _site_idx, device)
            return MO.cdr_redesign_ddg(_mint_head, wt_reps, mut_reps, changed,
                                       _mint_ymean, _mint_ystd)

        _mint_cluster_prep = _mint_prep(cluster)
        _mint_wt_prep = _mint_prep([(base, 1.0)])[0]

        def wt_ddg(toks):
            # design-vs-native MINT binding ddG on the WT antigen; <0 = tighter than native
            return _mint_ddg_one(_mint_wt_prep[0], _mint_wt_prep[1], toks)

        def score_across(toks):
            # rank designs by MINT reward across the panel (higher -ddg = binds tighter than native)
            s = [-_mint_ddg_one(ag, wr, toks) for ag, wr in _mint_cluster_prep]
            return sum(s) / len(s), min(s)

        def sample_cluster():
            raise SystemExit("--oracle mint requires --mechanism svdd")

        def sample_single():
            raise SystemExit("--oracle mint requires --mechanism svdd")
    else:
        # The mpnn reward path goes through the SAME _agg_reward helper as baddg/mint, so
        # worst/cvar aggregation is fully supported here — the old blanket guard was stale and
        # silently killed every obj_mpnn run at argparse. wt_weight genuinely is baddg-only
        # (it needs a WT binding ddG the mpnn log-ratio does not provide).
        if args.wt_weight:
            ap.error("--wt-weight is only implemented for --oracle baddg_ddg")
        if args.antifold_prior_weight:
            ap.error("--antifold-prior-weight is only implemented for --oracle baddg_ddg|consensus")

        def wt_ddg(toks):
            return None

        def score_across(toks):
            s = [gt_logratio(adapter.model, f, args.cdr_positions, toks, native, device) for f, _ in cluster]
            return sum(s) / len(s), min(s)

        def sample_cluster():
            return diffusion.sample_nos_mpnn_cluster_direct(
                struct, adapter.model, cluster, args.cdr_positions, args.guidance_scale)[:n_cdr]

        def sample_single():
            return diffusion.sample_nos_mpnn_direct(
                struct, adapter.model, base, args.cdr_positions, args.guidance_scale)[:n_cdr]

    # --- SVDD mechanism: replace the gradient sample_cluster/sample_single with derivative-free
    #     value-based decoding on the chosen --base generator. The oracle above still supplies the
    #     RANKING (score_across) and the BA-DDG grade; here it (re)appears only as the SVDD REWARD,
    #     evaluated on realized discrete CDRs (works for non-differentiable rewards too). ---
    if args.mechanism in ("svdd", "smc"):
        import math as _math
        from lyra_mutants.Model.diffusion import sample_svdd, sample_smc

        _pesc_map = None
        if args.steer_agg == "fused":
            if not args.pesc_json:
                raise SystemExit("--steer-agg fused requires --pesc-json")
            import json as _json
            _pesc_map = _json.load(open(args.pesc_json))

        def _agg_reward(vals, weights=None):
            v = torch.tensor(vals, dtype=torch.float)     # reward: higher = binds tighter
            if args.steer_agg == "mean":
                return float(v.mean())
            if args.steer_agg == "worst":
                return float(v.min())                      # worst-case robustness (min reward)
            if args.steer_agg == "fused" and weights is not None:
                # escape-likelihood-weighted CVaR at cvar_alpha over the WORST tail (lowest
                # reward = worst binding on that escape variant); weights = P(escape) measure.
                w = torch.tensor(weights, dtype=torch.float)
                order = torch.argsort(v)                   # ascending: worst reward first
                acc, num, den, alpha = 0.0, 0.0, 0.0, args.cvar_alpha
                for i in order.tolist():
                    take = min(float(w[i]), max(0.0, alpha - acc))
                    if take <= 0.0:
                        break
                    num += take * float(v[i]); den += take; acc += float(w[i])
                return num / den if den > 0 else float(v[int(order[0])])
            k = max(1, _math.ceil(args.cvar_alpha * v.numel()))
            return float(torch.sort(v).values[:k].mean())  # CVaR: mean of the worst tail

        if args.oracle in ("baddg_ddg", "consensus"):
            def make_reward(feats_list, names=None):
                w = None
                if args.steer_agg == "fused" and names is not None and _pesc_map is not None:
                    raw = [float(_pesc_map.get(nm, 0.0)) for nm in names]
                    s = sum(raw)
                    w = [x / s for x in raw] if s > 0 else None
                    cov = sum(1 for nm in names if nm in _pesc_map)
                    print(f"SVDD fused reward: pesc coverage {cov}/{len(names)} steer variants "
                          f"(weights={'ok' if w else 'FALLBACK-uniform'})", flush=True)
                def reward(toks):                          # higher = binds tighter than native
                    return _agg_reward([-float(bo.ddg(baddg, f, args.cdr_positions, toks, native, n_ab))
                                        for f, _ in feats_list], w)
                return reward
        elif args.oracle == "mint":
            # MINT reward on realized CDRs: CVaR (via _agg_reward) over the escape panel of
            # -ddg (higher = candidate binds the antigen variant tighter than native). Precompute the
            # per-variant WT reps once per closure; only the candidate's mutant forward is per-call.
            def make_reward(feats_list, names=None):
                prep = _mint_prep(feats_list)
                def reward(toks):                          # higher = binds tighter than native
                    return _agg_reward([-_mint_ddg_one(ag, wr, toks) for ag, wr in prep])
                return reward
        else:  # mpnn
            def make_reward(feats_list, names=None):
                def reward(toks):
                    return _agg_reward([gt_logratio(adapter.model, f, args.cdr_positions,
                                                    toks, native, device)
                                        for f, _ in feats_list])
                return reward

        _reward_cluster = make_reward(cluster, [nm for nm, _ in labels])
        _reward_single = make_reward([(base, 1.0)], None)
        _base_dist_fn = build_base_dist(args, diffusion, struct, adapter, base, native,
                                        args.cdr_positions, device)

        # cost counters: reward_evals is the hardware-independent budget axis that makes this
        # arm comparable to search baselines (GA / best-of-N) that also pay per reward call;
        # oracle_calls counts the per-variant model evaluations inside one reward.
        def _count_reward(fn, n_variants):
            def wrapped(toks):
                prof.bump("reward_evals")
                prof.bump("oracle_calls", n_variants)
                return fn(toks)
            return wrapped

        def _count_base(fn):
            def wrapped(tokens, t):
                prof.bump("base_forwards")
                return fn(tokens, t)
            return wrapped

        _reward_cluster = _count_reward(_reward_cluster, len(cluster))
        _reward_single = _count_reward(_reward_single, 1)
        _base_dist_fn = _count_base(_base_dist_fn)
        _n_ts = diffusion.n_timesteps
        if args.mechanism == "smc":
            print(f"mechanism: twisted-SMC (base={args.base}, particles={args.n_particles}, "
                  f"temp={args.smc_temp}, reward={args.oracle}/{args.steer_agg}, cluster={len(cluster)}, "
                  f"n_timesteps={_n_ts})", flush=True)

            def sample_cluster():
                return sample_smc(_base_dist_fn, _reward_cluster, n_cdr, _n_ts, device,
                                  n_particles=args.n_particles, smc_temp=args.smc_temp)[0]

            def sample_single():
                return sample_smc(_base_dist_fn, _reward_single, n_cdr, _n_ts, device,
                                  n_particles=args.n_particles, smc_temp=args.smc_temp)[0]
        else:
            print(f"mechanism: SVDD-PM (base={args.base}, k={args.svdd_k}, temp={args.svdd_temp}, "
                  f"reward={args.oracle}/{args.steer_agg}, cluster={len(cluster)}, n_timesteps={_n_ts})",
                  flush=True)

            _TRAJ = []

            def sample_cluster():
                tr = [] if args.trajectory_out else None
                out = sample_svdd(_base_dist_fn, _reward_cluster, n_cdr, _n_ts, device,
                                  svdd_k=args.svdd_k, svdd_temp=args.svdd_temp,
                                  base_temp=args.svdd_base_temp, trajectory=tr)
                if tr is not None:
                    _TRAJ.append({"design_index": len(_TRAJ), "steps": tr,
                                  "final_cdr": to_seq(out)})
                return out

            def sample_single():
                return sample_svdd(_base_dist_fn, _reward_single, n_cdr, _n_ts, device,
                                   svdd_k=args.svdd_k, svdd_temp=args.svdd_temp,
                                   base_temp=args.svdd_base_temp)

    # --- generate-and-filter: sample N cluster designs, rank by worst-case ---
    print(f"sampling {args.n_designs} {args.oracle}-steered designs "
          f"[mechanism={args.mechanism}] ...", flush=True)
    def _wtd(toks):
        v = wt_ddg(toks)
        return round(v, 3) if v is not None else None

    cand = []
    for i in range(args.n_designs):
        with prof.section("sample_cluster"):
            toks = sample_cluster()
        mean, worst = score_across(toks)
        cand.append((to_seq(toks), mean, worst, _wtd(toks), toks))
        print(f"  design {i+1}/{args.n_designs}: worst {worst:.2f} mean {mean:.2f}"
              + (f" wt_ddg {cand[-1][3]:+.2f}" if cand[-1][3] is not None else ""), flush=True)
    cand.sort(key=lambda c: -c[2])           # best worst-case first
    top = cand[:args.top_k]

    single = sample_single()
    s_mean, s_worst = score_across(single)

    # Diagnostic BA-DDG grading column (oracle=baddg_ddg only): per-design, per-variant
    # design-vs-native binding ddG across the FULL panel incl. held-out eval. <0 = design binds
    # tighter than native on that variant. This is the STEERER grading its own designs, so it is
    # IN-OBJECTIVE / circular — a separation diagnostic, NEVER a production verdict (cf.
    # consensus_arbiter: exclude the steerer from its own jury). Reuses bo.ddg, already exercised
    # by score_across/wt_ddg above; per-variant feats are already in memory (cluster+eval_cluster).
    def _baddg_grade(toks):
        if args.oracle not in ("baddg_ddg", "consensus"):
            return None
        g = {"WT": round(float(bo.ddg(baddg, base, args.cdr_positions, toks, native, n_ab)), 3)}
        for (name, _t), (f, _o) in list(zip(labels, cluster)) + list(zip(eval_labels, eval_cluster)):
            g[name] = round(float(bo.ddg(baddg, f, args.cdr_positions, toks, native, n_ab)), 3)
        return g

    designs = ([{"name": f"cluster_top{i+1}", "cdr_seq": seq,
                 "mpnn_mean": round(m, 3), "mpnn_worst": round(w, 3), "wt_ddg": wtd,
                 "baddg_grade": _baddg_grade(tk)}
                for i, (seq, m, w, wtd, tk) in enumerate(top)]
               + [{"name": "single", "cdr_seq": to_seq(single),
                   "mpnn_mean": round(s_mean, 3), "mpnn_worst": round(s_worst, 3),
                   "wt_ddg": _wtd(single), "baddg_grade": _baddg_grade(single)}])

    def _variant(name, toks, omega, split):
        return {"name": name, "omega": round(omega, 5), "split": split,
                "mutations": [{"chain": c, "wt": wt, "resnum": rn, "icode": ic, "mut": mut}
                              for (c, wt, rn, ic, mut) in toks]}

    # WT (used by Phase A and as the per-split WT-ag reference) + the steered panel; with
    # holdout, also the held-out eval panel. Each variant tagged with its split so the
    # downstream FoldX worst-case can be sliced steer vs held-out.
    antigen_variants = [{"name": "WT", "mutations": [], "split": "both"}]
    for (name, toks), (_, omega) in zip(labels, cluster):
        antigen_variants.append(_variant(name, toks, omega, "steer"))
    for (name, toks), (_, omega) in zip(eval_labels, eval_cluster):
        antigen_variants.append(_variant(name, toks, omega, "eval"))

    out = {
        "target": args.target,
        "pdb": os.path.basename(pdb_path),
        "antibody_chains": antibody,
        "antigen_chains": antigen_chains,
        "cdr_positions_flat": list(args.cdr_positions),
        "cdr_residues": cdr_residues,
        "native_cdr": native_cdr,
        "antibody_designs": designs,
        "antigen_variants": antigen_variants,
        "guidance_scale": args.guidance_scale,
        "mechanism": args.mechanism,
        "base": args.base if args.mechanism == "svdd" else "diffusion",
        "svdd_k": args.svdd_k if args.mechanism == "svdd" else None,
        "svdd_temp": args.svdd_temp if args.mechanism == "svdd" else None,
        # base proposal temperature (T<1 sharpens the prior BEFORE candidates are drawn).
        # Distinct from svdd_temp, which is the reward-resampling tau. Recorded so a
        # (T, K) sweep cell is identifiable from its export alone.
        "svdd_base_temp": args.svdd_base_temp,
        "oracle": args.oracle,
        "lam": args.lam if args.oracle == "consensus" else None,
        "mu": args.mu if args.oracle == "consensus" else None,
        "ptm_panel_k": args.ptm_panel_k if args.oracle == "ptm_energy" else None,
        "ptm_guided_steps": args.ptm_guided_steps if args.oracle == "ptm_energy" else None,
        # for ptm_energy with --ptm-guided-steps>0 the EFFECTIVE frac is derived from CDR length
        "guide_last_frac": (ptm_frac if args.oracle == "ptm_energy" and args.ptm_guided_steps > 0
                            else args.guide_last_frac),
        "prior_weight": args.prior_weight,
        "steer_agg": args.steer_agg,
        "cvar_alpha": args.cvar_alpha if args.steer_agg == "cvar" else None,
        "wt_weight": args.wt_weight,
        "holdout_by": args.holdout_by,
        "alignment_ok": bool(aligned and cdr_ok),
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)

    print(f"\nnative CDR: {native_cdr}")
    for d in designs:
        print(f"  {d['name']:<14} {d['cdr_seq']}  worst {d['mpnn_worst']:+.2f}  mean {d['mpnn_mean']:+.2f}")
    print(f"\nwrote {len(designs)} designs x {len(antigen_variants)} antigen variants -> {args.out}")
    prof.bump("designs_out", len(designs))
    prof.write(args.out)

    if args.trajectory_out and args.mechanism == "svdd":
        os.makedirs(os.path.dirname(os.path.abspath(args.trajectory_out)), exist_ok=True)
        with open(args.trajectory_out, "w") as fh:
            json.dump({"target": args.target, "base": args.base, "oracle": args.oracle,
                       "svdd_k": args.svdd_k, "svdd_temp": args.svdd_temp,
                       "steer_agg": args.steer_agg, "cvar_alpha": args.cvar_alpha,
                       "native_cdr": to_seq(native), "n_cdr": n_cdr,
                       "guided": args.svdd_k > 1,
                       "trajectories": _TRAJ}, fh)
        print(f"wrote {len(_TRAJ)} inference trajectories -> {args.trajectory_out}", flush=True)


if __name__ == "__main__":
    main()
