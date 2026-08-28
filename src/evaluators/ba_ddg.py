#!/usr/bin/env python3
"""BA-DDG energy oracle for diffusion steering — both rungs of Stage B.

The lyra steering loop (diffusion.py) needs two things from an oracle:
  B1 (baddg_if)  : a per-position [n_cdr, 20] reward table evaluated at the current hard
                   CDR guess  -> drop-in for _mpnn_cdr_logprobs, just BA-DDG weights.
  B2 (baddg_ddg) : a DIFFERENTIABLE scalar binding energy of the SOFT CDR distribution,
                   so the gradient flows through BA-DDG itself (the "real" upgrade).

BA-DDG is a fine-tuned ProteinMPNN + a thermodynamic-cycle ΔΔG head (complex score minus
single-chain score). For steering one CDR we use the Boltzmann-aligned BINDING score:
    binding(soft_cdr, antigen) = E_logP_complex[CDR] - E_logP_antibody_alone[CDR]
where E_logP[CDR] = Σ_p Σ_a soft[p,a]·logP(a | structure, context). Higher = the CDR is
more interface-preferred (the single-chain term removes pure antibody-fold preference, which
is what makes it a *binding* signal rather than plain inverse-folding likelihood).

Both antibody chains parse to chain_nb=0 in BA-DDG, but lyra's feats use chain_encoding_all
where antibody chains are 1..n_ab and antigen is >n_ab; we pass antibody_max_enc=n_ab so the
single-chain (antibody-alone) mask is chain_encoding_all <= n_ab.

Runs under system python3. Imports BA-DDG training modules.
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys

import torch
import torch.nn.functional as F
from easydict import EasyDict
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

BADDG = f"{VAAD_ROOT}/tools/BA-DDG"
_TRAIN = os.path.join(BADDG, "training")


def _load_baddg_isolated():
    """Load BA-DDG's protein_mpnn_utils + ddg_predictor under unique module names so they
    don't collide with vendor/ProteinMPNN's identically-named `protein_mpnn_utils` (which the
    lyra adapter imports first). BA-DDG's protein_mpnn_utils imports common_utils, so we put
    BA-DDG/training on sys.path transiently, and alias `protein_mpnn_utils` to BA-DDG's only
    while ddg_predictor (which does `from protein_mpnn_utils import ...`) executes."""
    added = _TRAIN not in sys.path
    if added:
        sys.path.append(_TRAIN)
    spec = importlib.util.spec_from_file_location(
        "baddg_pmpnn", os.path.join(_TRAIN, "protein_mpnn_utils.py"))
    pmpnn = importlib.util.module_from_spec(spec)
    sys.modules["baddg_pmpnn"] = pmpnn
    spec.loader.exec_module(pmpnn)

    saved = sys.modules.get("protein_mpnn_utils")
    sys.modules["protein_mpnn_utils"] = pmpnn
    try:
        spec2 = importlib.util.spec_from_file_location(
            "baddg_ddg_predictor", os.path.join(_TRAIN, "ddg_predictor.py"))
        ddgp = importlib.util.module_from_spec(spec2)
        sys.modules["baddg_ddg_predictor"] = ddgp
        spec2.loader.exec_module(ddgp)
    finally:
        if saved is not None:
            sys.modules["protein_mpnn_utils"] = saved
        else:
            sys.modules.pop("protein_mpnn_utils", None)
    return pmpnn, ddgp


_pmpnn, _ddgp = _load_baddg_isolated()
DDGPredictor = _ddgp.DDGPredictor
gather_nodes = _pmpnn.gather_nodes
cat_neighbors_nodes = _pmpnn.cat_neighbors_nodes

N_AA = 20
VOCAB = 21
CKPT_DDG = os.path.join(BADDG, "ckpt", "ddg_model.ckpt")
_CFG = EasyDict({
    "ca_only": False, "hidden_dim": 128, "num_layers": 3,
    "backbone_noise": 0.0, "num_edges": 48, "loss_weight_boltzmann": 1.0,
})


def load_baddg_mpnn(device, fold=0):
    """Return BA-DDG's fine-tuned ProteinMPNN (the .mpnn submodule), eval, frozen."""
    model = DDGPredictor(_CFG)
    sd = torch.load(CKPT_DDG, map_location="cpu", weights_only=False)
    model.load_state_dict(sd["models"][fold], strict=False)
    mpnn = model.mpnn.to(device).eval()
    for p in mpnn.parameters():
        p.requires_grad_(False)
    return mpnn


# --------------------------------------------------------------------- B1: reward table
@torch.no_grad()
def cdr_logprobs_table(mpnn, feats, cdr_positions, cdr_tokens):
    """[n_cdr, 20] BA-DDG log-probs at CDR positions, CDR set designable (decoded last,
    conditioned on the full complex incl. antigen S). Drop-in for _mpnn_cdr_logprobs."""
    S = feats["S"].clone()
    S[0, cdr_positions] = cdr_tokens
    chain_M = torch.zeros_like(feats["chain_M"], dtype=torch.float)
    chain_M[0, cdr_positions] = 1.0
    lp = mpnn.deterministic_forward(
        feats["X"], S, feats["mask"], chain_M,
        feats["residue_idx"], feats["chain_encoding_all"])
    return lp[0, cdr_positions, :N_AA]


@torch.no_grad()
def if_logratio(mpnn, feats, cdr_positions, toks, native):
    """BA-DDG inverse-folding CDR log-ratio (design vs native), scalar. B1 ranking."""
    lp = cdr_logprobs_table(mpnn, feats, cdr_positions, toks)
    return float(sum(lp[k, int(toks[k])] - lp[k, int(native[k])]
                     for k in range(len(cdr_positions))))


# --------------------------------------------------------- soft (differentiable) forward
def _decoding_order(chain_M, mask, device, generator=None):
    """argsort so fixed positions (chain_M=0) decode first, designable (=1) decode last."""
    chain_M = chain_M * mask
    rand = torch.abs(torch.randn(chain_M.shape, device=device, generator=generator))
    return torch.argsort((chain_M + 1e-4) * rand)


def soft_log_probs(mpnn, X, soft_onehot, mask, chain_M, residue_idx, chain_encoding_all,
                   decoding_order):
    """deterministic_forward, but the sequence enters as soft_onehot[1,L,VOCAB] @ W_s.weight
    so log-probs are differentiable in soft_onehot. decoding_order is supplied (fixed within
    a denoising step) for stable multi-step gradients. Returns log_probs [1, L, VOCAB]."""
    device = X.device
    E, E_idx = mpnn.features(X, mask, residue_idx, chain_encoding_all)
    h_V = torch.zeros((E.shape[0], E.shape[1], E.shape[-1]), device=device)
    h_E = mpnn.W_e(E)

    mask_attend = gather_nodes(mask.unsqueeze(-1), E_idx).squeeze(-1)
    mask_attend = mask.unsqueeze(-1) * mask_attend
    for layer in mpnn.encoder_layers:
        h_V, h_E = layer(h_V, h_E, E_idx, mask, mask_attend)

    h_S = soft_onehot @ mpnn.W_s.weight            # [1, L, hidden]  (soft embedding)
    h_ES = cat_neighbors_nodes(h_S, h_E, E_idx)
    h_EX_encoder = cat_neighbors_nodes(torch.zeros_like(h_S), h_E, E_idx)
    h_EXV_encoder = cat_neighbors_nodes(h_V, h_EX_encoder, E_idx)

    mask_size = E_idx.shape[1]
    perm = F.one_hot(decoding_order, num_classes=mask_size).float()
    order_mask_backward = torch.einsum(
        "ij, biq, bjp->bqp",
        (1 - torch.triu(torch.ones(mask_size, mask_size, device=device))), perm, perm)
    mask_attend = torch.gather(order_mask_backward, 2, E_idx).unsqueeze(-1)
    mask_1D = mask.view([mask.size(0), mask.size(1), 1, 1])
    mask_bw = mask_1D * mask_attend
    mask_fw = mask_1D * (1.0 - mask_attend)

    h_EXV_encoder_fw = mask_fw * h_EXV_encoder
    for layer in mpnn.decoder_layers:
        h_ESV = cat_neighbors_nodes(h_V, h_ES, E_idx)
        h_ESV = mask_bw * h_ESV + h_EXV_encoder_fw
        h_V = layer(h_V, h_ESV, mask)

    return F.log_softmax(mpnn.W_out(h_V), dim=-1)


def _soft_onehot(feats, cdr_positions, probs):
    """[1, L, VOCAB] hard one-hot of feats['S'] with CDR rows replaced by `probs` (soft)."""
    S = feats["S"].clamp(max=N_AA)
    soft = F.one_hot(S, VOCAB).float()             # [1, L, VOCAB]
    cdr_block = torch.zeros(len(cdr_positions), VOCAB, device=probs.device, dtype=probs.dtype)
    cdr_block[:, :N_AA] = probs
    soft = soft.clone()
    soft[0, cdr_positions] = cdr_block
    return soft


def binding_energy_soft(mpnn, feats, cdr_positions, probs, antibody_max_enc,
                        decoding_order, ab_logp_cache=None):
    """Differentiable BA-DDG binding score of the soft CDR `probs`[n_cdr,20] against the
    antigen in feats. Returns (binding_scalar, ab_logp) where binding = complex_score -
    antibody_alone_score (higher = tighter). ab_logp can be cached/reused across antigen
    variants (antibody-alone is antigen-independent)."""
    device = probs.device
    soft = _soft_onehot(feats, cdr_positions, probs)
    chain_M = torch.zeros_like(feats["chain_M"], dtype=torch.float)
    chain_M[0, cdr_positions] = 1.0

    lp_complex = soft_log_probs(
        mpnn, feats["X"], soft, feats["mask"], chain_M,
        feats["residue_idx"], feats["chain_encoding_all"], decoding_order)

    def cdr_score(lp):
        idx = torch.as_tensor(cdr_positions, device=device, dtype=torch.long)
        return (probs * lp[0, idx, :N_AA]).sum()

    complex_score = cdr_score(lp_complex)

    if ab_logp_cache is None:
        ab_mask = (feats["chain_encoding_all"] <= antibody_max_enc).float() * feats["mask"]
        X_ab = feats["X"] * ab_mask[..., None, None]
        order_ab = _decoding_order(chain_M, ab_mask, device)
        ab_logp = soft_log_probs(
            mpnn, X_ab, soft, ab_mask, chain_M * ab_mask,
            feats["residue_idx"], feats["chain_encoding_all"], order_ab)
    else:
        ab_logp = ab_logp_cache
    ab_score = cdr_score(ab_logp)
    return complex_score - ab_score, ab_logp


def make_ddg_energy_fn(mpnn, cluster_feats, cdr_positions, antibody_max_enc,
                       agg="mean", cvar_alpha=0.2, wt_feats=None, wt_weight=0.0,
                       antifold_table=None, antifold_weight=0.0, share_ab_logp=True):
    """Build energy_fn(probs[n_cdr,20], t) -> scalar cluster energy to MINIMIZE, differentiable
    in probs (gradient flows through BA-DDG). The antibody-alone term is antigen-independent, so
    it's computed once per call and reused across variants; the decoding order is seeded by the
    denoising step t (stable across the 3 inner grad steps).

    `agg` chooses how the per-member energies e_i = -binding_i are combined across the escape
    panel (higher e_i = weaker binding = worse):
      - "mean" : omega-weighted mean  Σ (ω_i/Σω)·e_i   (the original behaviour)
      - "worst": single worst member   max_i e_i        (pure worst-case; jumpy with FoldX noise)
      - "cvar" : mean of the worst ceil(cvar_alpha·n) members (CVaR_α; the worst tail, smoothed).
                 α=1.0 -> unweighted mean; α->0 -> "worst". Default α=0.2.
    "worst"/"cvar" align steering with the worst-case Phase B objective (train==eval); "mean"
    can leave single-mutant escape holes (see PAPER_PIPELINE_REPORT.md §12.4).

    `wt_weight`>0 adds a WT-antigen binding term  + wt_weight·(-binding_WT)  using `wt_feats`
    (the WT complex), so steering keeps one foot on WT-antigen potency. This matters most under
    worst/cvar, where the near-WT panel members drop out of the gradient and WT would otherwise
    be unconstrained (see §12.6). Pure robustness = wt_weight 0; potency-anchored = wt_weight>0.

    `antifold_weight`>0 adds an AntiFold PRIOR (product-of-experts): + antifold_weight·(-Σ probs·logP_af),
    where `antifold_table`[n_cdr,20] holds AntiFold's antigen-conditioned log-probs (unusable CDR
    positions zeroed, so they contribute nothing). Minimising this pulls the design toward AntiFold's
    (SOTA inverse-folding) sequence distribution while BA-DDG supplies the binding/robustness signal
    AntiFold is blind to — so AntiFold becomes the base/floor and steering is the value-add. The term
    is antigen-independent, so it's applied ONCE to the aggregated energy, not per escape member.

    `share_ab_logp` caches the antibody-alone term from the FIRST member and reuses it across all
    members — valid ONLY when every member shares one backbone (frozen backbone: antibody coords are
    identical). Set False when members carry per-variant re-folded backbones (build_cluster
    refold_dir), where each variant's antibody coordinates differ and the cache would be wrong."""
    device = cluster_feats[0][0]["X"].device
    total = sum(o for _, o in cluster_feats) + 1e-8
    cdr_idx = torch.as_tensor(cdr_positions, device=device, dtype=torch.long)
    f0 = cluster_feats[0][0]
    chain_M = torch.zeros_like(f0["chain_M"], dtype=torch.float)
    chain_M[0, cdr_idx] = 1.0
    af_tab = None if antifold_table is None else antifold_table.to(device=device, dtype=torch.float)

    def energy_fn(probs, t):
        gen = torch.Generator(device=device)
        gen.manual_seed(int(t))
        order_c = _decoding_order(chain_M, f0["mask"], device, gen)
        energies, weights, ab_logp = [], [], None
        for f, o in cluster_feats:
            b, ab = binding_energy_soft(
                mpnn, f, cdr_positions, probs, antibody_max_enc, order_c,
                ab_logp_cache=(ab_logp if share_ab_logp else None))
            if share_ab_logp and ab_logp is None:
                ab_logp = ab
            energies.append(-b)             # e_i; higher = weaker binding = worse
            weights.append(o)
        e_stack = torch.stack(energies)     # [n_members], differentiable
        if agg == "mean" or len(energies) == 1:
            w = torch.tensor(weights, device=device, dtype=e_stack.dtype) / total
            e = (w * e_stack).sum()
        elif agg == "worst":
            e = e_stack.max()
        elif agg == "cvar":
            k = max(1, int(math.ceil(cvar_alpha * e_stack.numel())))
            e = torch.topk(e_stack, k).values.mean()   # mean of the worst-k members
        else:
            raise ValueError(f"unknown agg={agg!r} (mean|worst|cvar)")
        if wt_weight > 0 and wt_feats is not None:   # keep WT-antigen binding in the objective
            b_wt, _ = binding_energy_soft(
                mpnn, wt_feats, cdr_positions, probs, antibody_max_enc, order_c,
                ab_logp_cache=ab_logp)
            e = e + wt_weight * (-b_wt)
        if antifold_weight > 0 and af_tab is not None:   # AntiFold product-of-experts prior
            e = e + antifold_weight * (-(probs * af_tab).sum())
        return e

    return energy_fn


# ------------------------------------------------------- B3: consensus steering energy
def make_consensus_energy_fn(mpnn, cluster_feats, cdr_positions, antibody_max_enc,
                             mpnn_tables, native_tokens, lam=1.0, agg="worst",
                             cvar_alpha=0.2, wt_feats=None, wt_mpnn_table=None,
                             wt_weight=0.0, z_baddg=1.0, z_mpnn=1.0,
                             antifold_table=None, antifold_weight=0.0, share_ab_logp=True):
    """Consensus steering energy (Stage-B B3 / EXPB). Per escape-panel member i, combine the
    differentiable BA-DDG binding energy with a ProteinMPNN complex-conditioned inverse-folding
    log-ratio term, each scaled and summed:

        e_i = z_baddg * (-binding_i)  +  lam * z_mpnn * (-mpnn_logratio_i)

    Both oriented so HIGHER = weaker/worse (energy to MINIMISE), then aggregated across the
    panel (worst / cvar / mean) exactly like make_ddg_energy_fn. This is the differentiable
    steering analogue of EXPB's ranking consensus z(baddg_w)+lam*z(mpnn_w)+mu*z(h3_w).

    mpnn_tables[i] : precomputed [n_cdr,20] vendor-ProteinMPNN complex log-prob table for
    cluster member i, in that member's antigen context at the NATIVE CDR (no_grad; supplied by
    the caller, which holds the vendor adapter). The mpnn term is then linear in probs:
        mpnn_logratio = Σ_p (probs - native_onehot)·logp_p
    differentiable in probs, mirroring the direct-MPNN steering energy
    (diffusion.sample_nos_mpnn_cluster_direct) but as a FIXED native-context reward table
    (sample_nos_energy's energy_fn only sees probs, not the running hard guess).

    NOTE — the H3-DDG consensus term (mu) is NOT wired: H3-DDG is a grader-only hypergraph
    model with no differentiable soft-CDR forward in this codebase, so it cannot steer. The
    caller passes/reports mu but no h3 gradient is applied here.
    """
    device = cluster_feats[0][0]["X"].device
    total = sum(o for _, o in cluster_feats) + 1e-8
    cdr_idx = torch.as_tensor(cdr_positions, device=device, dtype=torch.long)
    f0 = cluster_feats[0][0]
    chain_M = torch.zeros_like(f0["chain_M"], dtype=torch.float)
    chain_M[0, cdr_idx] = 1.0
    native_1h = F.one_hot(native_tokens.to(torch.long).clamp(max=N_AA - 1), N_AA).float().to(device)
    af_tab = None if antifold_table is None else antifold_table.to(device=device, dtype=torch.float)

    def _mpnn_term(probs, lp_table):
        # -(design - native) complex-IF log-ratio; higher = design less interface-preferred = worse
        return -((probs - native_1h) * lp_table).sum()

    def energy_fn(probs, t):
        gen = torch.Generator(device=device)
        gen.manual_seed(int(t))
        order_c = _decoding_order(chain_M, f0["mask"], device, gen)
        energies, weights, ab_logp = [], [], None
        for (f, o), lp_tab in zip(cluster_feats, mpnn_tables):
            b, ab = binding_energy_soft(
                mpnn, f, cdr_positions, probs, antibody_max_enc, order_c,
                ab_logp_cache=(ab_logp if share_ab_logp else None))
            if share_ab_logp and ab_logp is None:
                ab_logp = ab
            energies.append(z_baddg * (-b) + lam * z_mpnn * _mpnn_term(probs, lp_tab))
            weights.append(o)
        e_stack = torch.stack(energies)
        if agg == "mean" or len(energies) == 1:
            w = torch.tensor(weights, device=device, dtype=e_stack.dtype) / total
            e = (w * e_stack).sum()
        elif agg == "worst":
            e = e_stack.max()
        elif agg == "cvar":
            k = max(1, int(math.ceil(cvar_alpha * e_stack.numel())))
            e = torch.topk(e_stack, k).values.mean()
        else:
            raise ValueError(f"unknown agg={agg!r} (mean|worst|cvar)")
        if wt_weight > 0 and wt_feats is not None:
            b_wt, _ = binding_energy_soft(
                mpnn, wt_feats, cdr_positions, probs, antibody_max_enc, order_c, ab_logp_cache=ab_logp)
            e_wt = z_baddg * (-b_wt)
            if wt_mpnn_table is not None:
                e_wt = e_wt + lam * z_mpnn * _mpnn_term(probs, wt_mpnn_table)
            e = e + wt_weight * e_wt
        if antifold_weight > 0 and af_tab is not None:   # AntiFold product-of-experts prior
            e = e + antifold_weight * (-(probs * af_tab).sum())
        return e

    return energy_fn


# ------------------------------------------------------------------ B2 ranking: ddG
@torch.no_grad()
def ddg(mpnn, feats, cdr_positions, toks, native, antibody_max_enc):
    """BA-DDG thermodynamic-cycle ΔΔG of (design CDR) vs (native CDR) on this antigen,
    using hard sequences. Lower = design binds tighter than native. B2 ranking."""
    one = lambda t: _hard_binding(mpnn, feats, cdr_positions, t, antibody_max_enc)  # noqa: E731
    return -(one(toks) - one(native))   # higher binding -> more negative ddG


@torch.no_grad()
def _hard_binding(mpnn, feats, cdr_positions, toks, antibody_max_enc):
    """complex_score - antibody_alone_score for a HARD CDR (sum of log-probs at CDR)."""
    device = feats["X"].device
    S = feats["S"].clone()
    S[0, cdr_positions] = toks.to(S.dtype)
    chain_M = torch.zeros_like(feats["chain_M"], dtype=torch.float)
    chain_M[0, cdr_positions] = 1.0
    idx = torch.as_tensor(cdr_positions, device=device, dtype=torch.long)

    lp_c = mpnn.deterministic_forward(
        feats["X"], S, feats["mask"], chain_M,
        feats["residue_idx"], feats["chain_encoding_all"])
    cs = sum(lp_c[0, p, int(toks[k])] for k, p in enumerate(cdr_positions))

    ab_mask = (feats["chain_encoding_all"] <= antibody_max_enc).float() * feats["mask"]
    X_ab = feats["X"] * ab_mask[..., None, None]
    lp_a = mpnn.deterministic_forward(
        X_ab, S, ab_mask, chain_M * ab_mask,
        feats["residue_idx"], feats["chain_encoding_all"])
    as_ = sum(lp_a[0, p, int(toks[k])] for k, p in enumerate(cdr_positions))
    return float(cs - as_)
