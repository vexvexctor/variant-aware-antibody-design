from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import torch

from ..preprocessing.mutants import parse_mutation_token, sequence_position_index
from ..preprocessing.scoring import score_one_structural
from ..preprocessing.types import Mutant, MutantScore, ResidueSegment

# Add proteinnpt/utils to sys.path so we import esm directly as a standalone
# package, without triggering proteinnpt/__init__.py which pulls in tranception
# and requires the transformers library.
_ESM_UTILS = Path(__file__).resolve().parents[3] / "vendor" / "ProteinNPT" / "proteinnpt" / "utils"
if str(_ESM_UTILS) not in sys.path:
    sys.path.insert(0, str(_ESM_UTILS))


def esm1v_batch_scorer(
    wildtype_sequence: str,
    segments: list[ResidueSegment],
    checkpoint_path: Path,
    device: str | None = None,
    lower_bound: float = -1.0,
    upper_bound: float = 1.5,
) -> Callable[[list[Mutant]], list[MutantScore]]:
    """Return a BatchScorer backed by ESM1v wildtype-marginal scoring.

    The wildtype token probabilities are computed once at construction time
    (~L forward passes for a sequence of length L). Each subsequent call
    scores an entire batch via cheap table lookups — no additional forward
    passes per mutant.

    Parameters
    ----------
    wildtype_sequence: Antigen wildtype sequence (single-letter codes).
    segments:         Residue segments used for structural scores
                      (contact preservation, epitope change, etc.).
    checkpoint_path:  Path to the ESM1v .pt checkpoint file.
    device:           "cuda", "mps", or "cpu". Auto-detected when None.
    """
    device = _resolve_device(device)
    token_probs, alphabet = _compute_wildtype_token_probs(
        wildtype_sequence, checkpoint_path, device
    )
    pos_to_idx = sequence_position_index(segments)

    def _score_batch(mutants: list[Mutant]) -> list[MutantScore]:
        results = []
        for mutant in mutants:
            structural = score_one_structural(mutant, segments)
            esm_score = _esm1v_delta_log_likelihood(
                mutant.mutations, wildtype_sequence, pos_to_idx, token_probs, alphabet
            )
            results.append(_apply_esm_score(structural, esm_score, lower_bound, upper_bound))
        return results

    return _score_batch


# ---------------------------------------------------------------------------
# ESM1v internals
# ---------------------------------------------------------------------------

def _resolve_device(device: str | None) -> str:
    if device is not None:
        return device
    if torch.backends.mps.is_available():
        return "mps"
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"


def _compute_wildtype_token_probs(
    sequence: str,
    checkpoint_path: Path,
    device: str,
) -> tuple[torch.Tensor, object]:
    """Mask each position of the wildtype sequence and collect log-softmax outputs.

    Returns token_probs of shape (1, L+2, vocab_size) where L+2 accounts for
    the BOS and EOS tokens that ESM1v prepends/appends.
    """
    import argparse
    from esm.pretrained import load_model_and_alphabet

    torch.serialization.add_safe_globals([argparse.Namespace])
    model, alphabet = load_model_and_alphabet(str(checkpoint_path))
    model.eval().to(device)

    batch_converter = alphabet.get_batch_converter()
    _, _, batch_tokens = batch_converter([("wildtype", sequence)])
    batch_tokens = batch_tokens.to(device)

    seq_len = batch_tokens.size(1)  # includes BOS + EOS
    if seq_len > 1024:
        raise ValueError(
            f"Sequence length {seq_len} exceeds ESM1v's 1024-token limit. "
            "Use a shorter subsequence covering only the mutable region."
        )

    all_token_probs: list[torch.Tensor] = []
    for i in range(seq_len):
        masked = batch_tokens.clone()
        masked[0, i] = alphabet.mask_idx
        with torch.no_grad():
            logits = model(masked)["logits"]
        all_token_probs.append(
            torch.log_softmax(logits[0, i], dim=-1).cpu()
        )

    # Shape: (1, L+2, vocab_size)
    token_probs = torch.stack(all_token_probs, dim=0).unsqueeze(0)
    return token_probs, alphabet


def _esm1v_delta_log_likelihood(
    mutations: tuple[str, ...],
    wildtype_sequence: str,
    pos_to_idx: dict[str, int],
    token_probs: torch.Tensor,
    alphabet: object,
) -> float:
    """Sum delta log-likelihoods across all mutation positions.

    Positive score = mutation is more likely than wildtype under ESM1v
    (i.e. evolutionarily plausible change).
    Negative score = mutation is less likely (potentially disruptive).
    """
    score = 0.0
    for token in mutations:
        position, wt_aa, mut_aa = parse_mutation_token(token)
        seq_idx = pos_to_idx[position]
        # +1 because token_probs includes BOS at position 0
        token_idx = seq_idx + 1
        wt_encoded = alphabet.get_idx(wt_aa)
        mut_encoded = alphabet.get_idx(mut_aa)
        score += (
            token_probs[0, token_idx, mut_encoded]
            - token_probs[0, token_idx, wt_encoded]
        ).item()
    return score


def _apply_esm_score(
    structural: MutantScore,
    esm_score: float,
    lower_bound: float,
    upper_bound: float,
) -> MutantScore:
    n = len(structural.mutant.mutations)
    per_mutation = esm_score / n
    failed: list[str] = []
    if per_mutation < lower_bound:
        failed.append(f"esm1v_per_mutation_below_{lower_bound}")
    if per_mutation > upper_bound:
        failed.append(f"esm1v_per_mutation_above_{upper_bound}")
    return MutantScore(
        mutant=structural.mutant,
        proposal_model_score=structural.proposal_model_score,
        mutation_effect_score=esm_score,
        uncertainty=None,
        potts_score=None,
        esm_score=esm_score,
        predicted_binding_effect=None,
        accepted=not failed,
        failed_reasons=tuple(failed),
    )
