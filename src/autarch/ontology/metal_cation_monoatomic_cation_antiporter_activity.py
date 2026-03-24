"""metal cation:monoatomic cation antiporter activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.antiporter_activity import AntiporterActivity
from autarch.ontology.transport_utils import (
    coupled_transported_pairs,
    is_inorganic_cation,
    is_transition_metal_ion,
)


class MetalCationMonoatomicCationAntiporterActivity(AntiporterActivity):
    """metal cation:monoatomic cation antiporter activity."""

    GO_ID = "GO:0140828"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        transported = coupled_transported_pairs(reaction)
        has_metal_cation = any(is_transition_metal_ion(left_participant) for left_participant, _ in transported)
        has_other_cation = any(
            is_inorganic_cation(left_participant) and not is_transition_metal_ion(left_participant)
            for left_participant, _ in transported
        )
        if not (has_metal_cation and has_other_cation):
            return ClassificationResult(
                is_member=False,
                explanation="No metal-cation/monoatomic-cation antiport detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Metal cation:monoatomic cation antiporter activity",
        )
