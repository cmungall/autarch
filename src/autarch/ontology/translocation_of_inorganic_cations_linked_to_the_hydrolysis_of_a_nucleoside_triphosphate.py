"""translocation of inorganic cations linked to the hydrolysis of a nucleoside triphosphate.

Precision-biased EC-level translocase classifier for ATP/GTP-coupled transport
of directly represented inorganic cations.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.translocase_linked_to_hydrolysis import (
    TranslocaseLinkedToHydrolysis,
)
from autarch.ontology.transport_utils import (
    is_inorganic_cation,
    reactive_transported_pairs,
)


class TranslocationOfInorganicCationsLinkedToTheHydrolysisOfANucleosideTriphosphate(
    TranslocaseLinkedToHydrolysis
):
    """translocation of inorganic cations linked to the hydrolysis of a nucleoside triphosphate.

    Precision-biased EC-level translocase classifier for ATP/GTP-coupled
    transport of directly represented inorganic cations.
    """

    GO_ID: ClassVar[Optional[str]] = None
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "7.2.2.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for ATP/GTP-coupled translocation of inorganic cations."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        transported_pairs = reactive_transported_pairs(reaction)
        if not transported_pairs:
            return ClassificationResult(
                is_member=False,
                explanation="No transported substrate beyond ATPase-coupling participants",
            )

        if not all(
            is_inorganic_cation(left_participant)
            for left_participant, _ in transported_pairs
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Transported substrate is not consistently an inorganic cation",
            )

        transport_info = [
            f"{left_participant.chebi_id or 'unknown'}: "
            f"{left_participant.location or 'unknown'}→{right_participant.location or 'unknown'}"
            for left_participant, right_participant in transported_pairs
        ]
        return ClassificationResult(
            is_member=True,
            explanation=(
                "Translocation of inorganic cations linked to nucleoside-triphosphate "
                f"hydrolysis: {', '.join(transport_info)}"
            ),
        )
