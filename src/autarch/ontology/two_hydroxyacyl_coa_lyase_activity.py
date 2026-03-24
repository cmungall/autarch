"""2-hydroxyacyl-CoA lyase activity.

Catalysis of the cleavage of a 2-hydroxyacyl-CoA to form formyl-CoA and an
aldehyde product.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_FORMYL_COA = "CHEBI:57376"


class TwoHydroxyacylCoALyaseActivity(ReactionClass):
    """2-hydroxyacyl-CoA lyase activity.

    Catalysis of the cleavage of a 2-hydroxyacyl-CoA to form formyl-CoA and an
    aldehyde product.
    """

    GO_ID = "GO:0106359"
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
        if not any(participant.chebi_id == CHEBI_FORMYL_COA for participant in right):
            return ClassificationResult(
                is_member=False,
                explanation="Requires formyl-CoA as a cleavage product",
            )

        hydroxy_thioesters = [
            participant
            for participant in left
            if participant.is_thioester() and participant.has_moiety(Moiety.HYDROXYL)
        ]
        if not hydroxy_thioesters:
            return ClassificationResult(
                is_member=False,
                explanation="Requires a hydroxylated acyl-CoA donor",
            )

        aldehydes = [
            participant
            for participant in right
            if participant.has_moiety(Moiety.ALDEHYDE)
        ]
        if not aldehydes:
            return ClassificationResult(
                is_member=False,
                explanation="Requires aldehyde formation alongside formyl-CoA",
            )

        return ClassificationResult(
            is_member=True,
            explanation="2-hydroxyacyl-CoA lyase activity: cleavage of a hydroxylated acyl-CoA to formyl-CoA and an aldehyde",
        )
