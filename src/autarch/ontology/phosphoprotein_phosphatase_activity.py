"""phosphoprotein phosphatase activity.

Catalysis of phosphate removal from a phosphorylated protein substrate by
hydrolysis.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_PHOSPHATE
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

PHOSPHATE_IDS = {CHEBI_PHOSPHATE, "CHEBI:16838", "CHEBI:18367"}
PHOSPHOPROTEIN_IDS = {
    "CHEBI:61978",
    "CHEBI:83421",
    "CHEBI:43474",
    "CHEBI:29965",
    "CHEBI:64837",
    "CHEBI:30013",
    "CHEBI:83586",
}
DEPHOSPHOPROTEIN_IDS = {
    "CHEBI:46858",
    "CHEBI:29999",
    "CHEBI:29979",
    "CHEBI:43474",
    "CHEBI:83226",
    "CHEBI:61977",
}


class PhosphoproteinPhosphataseActivity(ReactionClass):
    """phosphoprotein phosphatase activity.

    Catalysis of phosphate removal from a phosphorylated protein substrate by
    hydrolysis.
    """

    GO_ID = "GO:0004721"
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
        if not any(participant.chebi_id == CHEBI_H2O for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires water as a hydrolytic substrate",
            )
        if not any(participant.chebi_id in PHOSPHATE_IDS for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires phosphate release",
            )

        left_ids = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        if not (left_ids & PHOSPHOPROTEIN_IDS):
            return ClassificationResult(
                is_member=False,
                explanation="No supported phosphoprotein substrate detected",
            )
        if not (right_ids & DEPHOSPHOPROTEIN_IDS):
            return ClassificationResult(
                is_member=False,
                explanation="No supported dephosphorylated protein product detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Phosphoprotein phosphatase activity: hydrolytic removal of phosphate from a protein substrate",
        )
