from __future__ import annotations

from .amino_acids import AMINO_ACIDS, chemistry_distance, proposal_log_probability
from .types import AminoAcidProposal, ResidueSegment


def propose_amino_acids(segments: list[ResidueSegment], top_k: int) -> dict[str, list[AminoAcidProposal]]:
    proposals: dict[str, list[AminoAcidProposal]] = {}
    for segment in segments:
        if not segment.mutable:
            continue
        ranked = sorted(
            (aa for aa in AMINO_ACIDS if aa != segment.wildtype),
            key=lambda aa: (
                abs(chemistry_distance(segment.wildtype, aa) - 0.16),
                aa,
            ),
        )
        position_proposals: list[AminoAcidProposal] = []
        for rank, amino_acid in enumerate(ranked[:top_k], start=1):
            logp = proposal_log_probability(rank)
            position_proposals.append(
                AminoAcidProposal(
                    position=segment.position,
                    wildtype=segment.wildtype,
                    amino_acid=amino_acid,
                    probability=pow(2.718281828, logp),
                    log_probability=logp,
                    rank=rank,
                )
            )
        proposals[segment.position] = position_proposals
    return proposals
