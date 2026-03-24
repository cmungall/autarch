"""sterol 12-alpha-hydroxylase activity.

Catalysis of 12-alpha hydroxylation of a sterol substrate using molecular oxygen
and reduced NADPH-hemoprotein reductase equivalents.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE = "CHEBI:57618"
CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE = "CHEBI:58210"
CHEBI_STEROL_DIOL = "CHEBI:28047"
CHEBI_STEROL_TRIOL = "CHEBI:16496"


class Sterol12AlphaHydroxylaseActivity(ReactionClass):
    """sterol 12-alpha-hydroxylase activity.

    Catalysis of 12-alpha hydroxylation of a sterol substrate using molecular oxygen
    and reduced NADPH-hemoprotein reductase equivalents.
    """

    GO_ID = "GO:0008397"
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
        if CHEBI_STEROL_DIOL not in left_ids or CHEBI_STEROL_TRIOL not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires the supported sterol diol to triol conversion",
            )
        if not {CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE, CHEBI_O2} <= left_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires reduced hemoprotein reductase and dioxygen",
            )
        if not {CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE, CHEBI_H2O, CHEBI_H_PLUS} <= right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires oxidized reductase, water, and hydron products",
            )
        return ClassificationResult(
            is_member=True,
            explanation="sterol 12-alpha-hydroxylase activity: monooxygenation of a sterol substrate",
        )
