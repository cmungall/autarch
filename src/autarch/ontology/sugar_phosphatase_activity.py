"""sugar-phosphatase activity.

Catalysis of the hydrolysis of a phosphate monoester bond in a sugar phosphate.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_PHOSPHATE
from autarch.ontology.lipid_utils import phosphate_count
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.ontology.transport_utils import is_carbohydrate_or_derivative

PHOSPHATE_IDS = {CHEBI_PHOSPHATE, "CHEBI:16838", "CHEBI:18367"}


class SugarPhosphataseActivity(ReactionClass):
    """sugar-phosphatase activity.

    Catalysis of the hydrolysis of a phosphate monoester bond in a sugar phosphate.
    """

    GO_ID = "GO:0050308"
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

        substrates = [
            participant
            for participant in left_core
            if is_carbohydrate_or_derivative(participant) and phosphate_count(participant) >= 1
        ]
        products = [
            participant
            for participant in right_core
            if is_carbohydrate_or_derivative(participant)
        ]
        if not substrates or not products:
            return ClassificationResult(
                is_member=False,
                explanation="Expected carbohydrate-like phosphorylated substrate and carbohydrate-like product",
            )

        for substrate in substrates:
            for product in products:
                if phosphate_count(substrate) == phosphate_count(product) + 1:
                    return ClassificationResult(
                        is_member=True,
                        explanation="Sugar-phosphatase activity: hydrolytic removal of one phosphate from a carbohydrate-like substrate",
                    )

        return ClassificationResult(
            is_member=False,
            explanation="No carbohydrate-like substrate/product pair loses exactly one phosphate",
        )
