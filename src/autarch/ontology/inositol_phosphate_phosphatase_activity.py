"""inositol phosphate phosphatase activity.

Catalysis of phosphate removal from a soluble inositol phosphate substrate by
hydrolysis.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_PHOSPHATE
from autarch.ontology.lipid_utils import (
    is_inositol_phosphate_like,
    phosphate_count,
)
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

PHOSPHATE_IDS = {CHEBI_PHOSPHATE, "CHEBI:16838", "CHEBI:18367"}


class InositolPhosphatePhosphataseActivity(ReactionClass):
    """inositol phosphate phosphatase activity.

    Catalysis of phosphate removal from a soluble inositol phosphate substrate
    by hydrolysis.
    """

    GO_ID = "GO:0052745"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward.is_member:
            return forward
        reverse = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
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

        left_core = [
            participant
            for participant in left_participants
            if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in right_participants
            if participant.chebi_id not in PHOSPHATE_IDS | {CHEBI_H_PLUS}
        ]

        substrates = [participant for participant in left_core if is_inositol_phosphate_like(participant)]
        products = [participant for participant in right_core if is_inositol_phosphate_like(participant)]
        if len(substrates) != 1 or len(products) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one soluble inositol phosphate substrate and one dephosphorylated inositol phosphate product",
            )

        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Unexpected additional organic participants",
            )

        if phosphate_count(products[0]) != phosphate_count(substrates[0]) - 1:
            return ClassificationResult(
                is_member=False,
                explanation="Inositol phosphate product does not lose exactly one phosphate group",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Inositol phosphate phosphatase activity: hydrolytic removal of one phosphate from a soluble inositol phosphate",
        )
