"""transition metal ion transmembrane transporter activity.

Enables the transfer of transition metal ions from one side of a membrane to
the other.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.metal_ion_transmembrane_transporter_activity import (
    MetalIonTransmembraneTransporterActivity,
)
from autarch.ontology.transport_utils import coupled_transported_pairs, is_transition_metal_ion


class TransitionMetalIonTransmembraneTransporterActivity(
    MetalIonTransmembraneTransporterActivity
):
    """transition metal ion transmembrane transporter activity.

    Enables the transfer of transition metal ions from one side of a membrane
    to the other.
    """

    GO_ID = "GO:0046915"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent = super().check_membership_impl(reaction)
        if not parent.is_member:
            return parent

        transported = coupled_transported_pairs(reaction)
        substantive = [
            left_participant
            for left_participant, _ in transported
            if left_participant.chebi_id != "CHEBI:15378"
        ]
        if not substantive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive transported ion remains after removing hydron coupling",
            )
        if not all(is_transition_metal_ion(participant) for participant in substantive):
            return ClassificationResult(
                is_member=False,
                explanation="Transported ion is not consistently transition-metal-like",
            )
        return ClassificationResult(
            is_member=True,
            explanation="Transition metal ion transmembrane transporter activity",
        )
