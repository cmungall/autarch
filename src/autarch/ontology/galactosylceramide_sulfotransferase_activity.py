"""galactosylceramide sulfotransferase activity.

Catalysis of sulfate transfer from PAPS to galactosylceramide, forming a
sulfated galactosylceramide product and PAP.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_GALACTOSYLCERAMIDE = "CHEBI:18390"
CHEBI_SULFO_GALACTOSYLCERAMIDE = "CHEBI:75956"
CHEBI_PAPS = "CHEBI:58339"
CHEBI_PAP = "CHEBI:58343"


class GalactosylceramideSulfotransferaseActivity(ReactionClass):
    """galactosylceramide sulfotransferase activity.

    Catalysis of sulfate transfer from PAPS to galactosylceramide, forming a
    sulfated galactosylceramide product and PAP.
    """

    GO_ID = "GO:0001733"
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
        left: list[Participant],
        right: list[Participant],
    ) -> ClassificationResult:
        left_ids = {participant.chebi_id for participant in left if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right if participant.chebi_id}
        if not {CHEBI_GALACTOSYLCERAMIDE, CHEBI_PAPS} <= left_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires galactosylceramide and PAPS as substrates",
            )
        if not {CHEBI_SULFO_GALACTOSYLCERAMIDE, CHEBI_PAP} <= right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires sulfated galactosylceramide and PAP as products",
            )
        if CHEBI_H_PLUS not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydron release")
        return ClassificationResult(
            is_member=True,
            explanation="galactosylceramide sulfotransferase activity: PAPS-dependent sulfation of galactosylceramide",
        )
