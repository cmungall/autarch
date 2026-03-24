"""3-beta-hydroxysteroid 3-dehydrogenase (NADP+) activity.

Catalysis of oxidation of a 3-beta-hydroxysteroid to the corresponding 3-oxo
steroid with NADP as acceptor.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_NADPH, CHEBI_NADP_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

HYDROXYSTEROID_IDS = {"CHEBI:18378", "CHEBI:1949"}
OXOSTEROID_IDS = {"CHEBI:16495", "CHEBI:136486"}


class ThreeBetaHydroxysteroid3DehydrogenaseNADPActivity(ReactionClass):
    """3-beta-hydroxysteroid 3-dehydrogenase (NADP+) activity.

    Catalysis of oxidation of a 3-beta-hydroxysteroid to the corresponding 3-oxo
    steroid with NADP as acceptor.
    """

    GO_ID = "GO:0000253"
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
        if not (HYDROXYSTEROID_IDS & left_ids):
            return ClassificationResult(is_member=False, explanation="Requires a supported 3-beta-hydroxysteroid substrate")
        if not (OXOSTEROID_IDS & right_ids):
            return ClassificationResult(is_member=False, explanation="Requires a supported 3-oxosteroid product")
        if CHEBI_NADP_PLUS not in left_ids or CHEBI_NADPH not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires NADP(+)/NADPH coupling")
        if CHEBI_H_PLUS not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydron release")
        return ClassificationResult(
            is_member=True,
            explanation="3-beta-hydroxysteroid 3-dehydrogenase (NADP+) activity: oxidation of a 3-beta-hydroxysteroid with NADP(+) reduction",
        )
