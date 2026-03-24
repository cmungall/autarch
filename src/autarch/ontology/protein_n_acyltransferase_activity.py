"""protein N-acyltransferase activity.

Catalysis of transfer of an acyl group from an acyl donor to a protein amino
group.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_COA
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

PROTEIN_AMINO_IDS = {
    "CHEBI:29969",
    "CHEBI:78597",
    "CHEBI:64718",
    "CHEBI:83690",
}
ACYLATED_PROTEIN_IDS = {
    "CHEBI:61930",
    "CHEBI:78598",
    "CHEBI:83683",
    "CHEBI:64738",
}


class ProteinNAcyltransferaseActivity(ReactionClass):
    """protein N-acyltransferase activity.

    Catalysis of transfer of an acyl group from an acyl donor to a protein amino
    group.
    """

    GO_ID = "GO:0140186"
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
        donor = next((participant for participant in left_participants if participant.is_thioester()), None)
        if donor is None:
            return ClassificationResult(
                is_member=False,
                explanation="Requires an acyl-thioester donor",
            )
        if not any(participant.chebi_id == CHEBI_COA for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires released coenzyme A",
            )

        left_ids = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        if not (left_ids & PROTEIN_AMINO_IDS):
            return ClassificationResult(
                is_member=False,
                explanation="No supported protein amino acceptor detected",
            )
        if not (right_ids & ACYLATED_PROTEIN_IDS):
            return ClassificationResult(
                is_member=False,
                explanation="No supported N-acyl protein product detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Protein N-acyltransferase activity: acyl-thioester donor transfers an acyl group to a protein amino group",
        )
