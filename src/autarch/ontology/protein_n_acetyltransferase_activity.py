"""protein N-acetyltransferase activity.

Catalysis of transfer of an acetyl group from acetyl-CoA to a protein amino
group.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_ACETYL_COA, CHEBI_COA
from autarch.ontology.protein_n_acyltransferase_activity import (
    ACYLATED_PROTEIN_IDS,
    PROTEIN_AMINO_IDS,
)
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class ProteinNAcetyltransferaseActivity(ReactionClass):
    """protein N-acetyltransferase activity.

    Catalysis of transfer of an acetyl group from acetyl-CoA to a protein amino
    group.
    """

    GO_ID = "GO:0034212"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward

        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse reaction orientation)",
            )
        return forward

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        left_ids = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        if CHEBI_ACETYL_COA not in left_ids or CHEBI_COA not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires acetyl-CoA consumption with CoA release",
            )
        if not (left_ids & PROTEIN_AMINO_IDS):
            return ClassificationResult(
                is_member=False,
                explanation="No supported protein amino acceptor detected",
            )
        if not (right_ids & ACYLATED_PROTEIN_IDS):
            return ClassificationResult(
                is_member=False,
                explanation="No supported N-acetyl protein product detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Protein N-acetyltransferase activity: acetyl-CoA donates an acetyl group to a protein amino group",
        )
