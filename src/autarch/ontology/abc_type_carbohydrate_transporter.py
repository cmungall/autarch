"""ABC-type carbohydrate transporter activity.

ATP-dependent transmembrane transport of carbohydrate substrates and closely
related carbohydrate derivatives.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.atpase_coupled_transmembrane_transporter import (
    ATPaseCoupledTransmembraneTransporter,
)
from autarch.ontology.transport_utils import (
    is_carbohydrate_or_derivative,
    reactive_transported_pairs,
)


class ABCTypeCarbohydrateTransporter(ATPaseCoupledTransmembraneTransporter):
    """ABC-type carbohydrate transporter activity.

    ATP-dependent transmembrane transport of carbohydrate substrates and closely
    related carbohydrate derivatives.
    """

    GO_ID: ClassVar[Optional[str]] = "GO:0043211"
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "7.5.2.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for ATP-dependent transport of carbohydrate substrates."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        transported_pairs = reactive_transported_pairs(reaction)
        if not transported_pairs:
            return ClassificationResult(
                is_member=False,
                explanation="No transported substrate beyond ATPase-coupling participants",
            )

        if not all(is_carbohydrate_or_derivative(left_participant) for left_participant, _ in transported_pairs):
            return ClassificationResult(
                is_member=False,
                explanation="Transported substrate is not consistently carbohydrate-like",
            )

        transport_info = [
            f"{left_participant.chebi_id or 'unknown'}: "
            f"{left_participant.location or 'unknown'}→{right_participant.location or 'unknown'}"
            for left_participant, right_participant in transported_pairs
        ]
        return ClassificationResult(
            is_member=True,
            explanation=f"ABC-type carbohydrate transport: {', '.join(transport_info)}",
        )
