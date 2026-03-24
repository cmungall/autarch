"""succinyltransferase activity.

Catalysis of transfer of a succinyl group from succinyl-CoA to an acceptor
substrate.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_COA, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

SUCCINYL_COA_IDS = {"CHEBI:15380", "CHEBI:57292"}


class SuccinyltransferaseActivity(ReactionClass):
    """succinyltransferase activity.

    Catalysis of transfer of a succinyl group from succinyl-CoA to an acceptor
    substrate.
    """

    GO_ID = "GO:0016748"
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
        if not any(participant.chebi_id in SUCCINYL_COA_IDS for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No succinyl-CoA donor detected",
            )
        if not any(participant.chebi_id == CHEBI_COA for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires released coenzyme A",
            )

        left_core = [
            participant
            for participant in left_participants
            if participant.chebi_id not in SUCCINYL_COA_IDS | {CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in right_participants
            if participant.chebi_id not in {CHEBI_COA, CHEBI_H_PLUS}
        ]
        if not left_core or not right_core:
            return ClassificationResult(
                is_member=False,
                explanation="Expected an acceptor substrate and a succinylated product",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Succinyltransferase activity: succinyl-CoA donor releases CoA during succinyl transfer",
        )
