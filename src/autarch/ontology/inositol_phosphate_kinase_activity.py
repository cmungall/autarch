"""inositol phosphate kinase activity.

Catalysis of ATP-dependent phosphorylation of a soluble inositol phosphate.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_H_PLUS
from autarch.ontology.lipid_utils import (
    is_inositol_phosphate_like,
    phosphate_count,
)
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class InositolPhosphateKinaseActivity(ReactionClass):
    """inositol phosphate kinase activity.

    Catalysis of ATP-dependent phosphorylation of a soluble inositol
    phosphate.
    """

    GO_ID = "GO:0180030"
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
        left_ids = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        if CHEBI_ATP not in left_ids or CHEBI_ADP not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires ATP consumption with ADP production",
            )

        left_core = [
            participant
            for participant in left_participants
            if participant.chebi_id not in {CHEBI_ATP, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in right_participants
            if participant.chebi_id not in {CHEBI_ADP, CHEBI_H_PLUS}
        ]

        substrates = [participant for participant in left_core if is_inositol_phosphate_like(participant)]
        products = [participant for participant in right_core if is_inositol_phosphate_like(participant)]
        if len(substrates) != 1 or len(products) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one soluble inositol phosphate substrate and one more highly phosphorylated inositol product",
            )

        if phosphate_count(products[0]) != phosphate_count(substrates[0]) + 1:
            return ClassificationResult(
                is_member=False,
                explanation="Inositol phosphate kinase products must gain exactly one phosphate group",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Inositol phosphate kinase activity: ATP-dependent phosphorylation of a soluble inositol phosphate",
        )
