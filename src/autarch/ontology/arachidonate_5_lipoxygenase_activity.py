"""arachidonate 5-lipoxygenase activity.

Catalysis of oxygenation of arachidonate to 5-hydroperoxy-eicosatetraenoate or
its dehydration product leukotriene A4.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_ARACHIDONATE = "CHEBI:32395"
CHEBI_5_HPETE = "CHEBI:57450"
CHEBI_LEUKOTRIENE_A4 = "CHEBI:57463"


class Arachidonate5LipoxygenaseActivity(ReactionClass):
    """arachidonate 5-lipoxygenase activity.

    Catalysis of oxygenation of arachidonate to 5-hydroperoxy-eicosatetraenoate or
    its dehydration product leukotriene A4.
    """

    GO_ID = "GO:0004051"
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

        if {CHEBI_ARACHIDONATE, CHEBI_O2} <= left_ids and CHEBI_5_HPETE in right_ids:
            return ClassificationResult(
                is_member=True,
                explanation="arachidonate 5-lipoxygenase activity: oxygenation of arachidonate to 5-HPETE",
            )
        if {CHEBI_ARACHIDONATE, CHEBI_O2} <= left_ids and {CHEBI_LEUKOTRIENE_A4, CHEBI_H2O} <= right_ids:
            return ClassificationResult(
                is_member=True,
                explanation="arachidonate 5-lipoxygenase activity: oxygenation/dehydration of arachidonate to leukotriene A4",
            )
        if CHEBI_5_HPETE in left_ids and {CHEBI_LEUKOTRIENE_A4, CHEBI_H2O} <= right_ids:
            return ClassificationResult(
                is_member=True,
                explanation="arachidonate 5-lipoxygenase activity: dehydration of 5-HPETE to leukotriene A4",
            )

        return ClassificationResult(
            is_member=False,
            explanation="Requires the supported arachidonate to 5-HPETE/leukotriene A4 conversion set",
        )
