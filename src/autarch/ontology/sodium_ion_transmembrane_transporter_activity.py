"""sodium ion transmembrane transporter activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.metal_ion_transmembrane_transporter_activity import (
    MetalIonTransmembraneTransporterActivity,
)
from autarch.ontology.transport_utils import coupled_transported_pairs, is_sodium_ion


class SodiumIonTransmembraneTransporterActivity(MetalIonTransmembraneTransporterActivity):
    """sodium ion transmembrane transporter activity."""

    GO_ID = "GO:0015081"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        transported = coupled_transported_pairs(reaction)
        if not any(is_sodium_ion(left_participant) for left_participant, _ in transported):
            return ClassificationResult(
                is_member=False,
                explanation="No transported sodium ion detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Sodium ion transmembrane transporter activity",
        )
