"""ATPase-coupled monoatomic cation transmembrane transporter activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.atpase_coupled_transmembrane_transporter import (
    ATPaseCoupledTransmembraneTransporter,
)
from autarch.ontology.transport_utils import coupled_transported_pairs, is_metal_ion


class ATPaseCoupledMonoatomicCationTransmembraneTransporterActivity(
    ATPaseCoupledTransmembraneTransporter
):
    """ATPase-coupled monoatomic cation transmembrane transporter activity."""

    GO_ID = "GO:0019829"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        transported = coupled_transported_pairs(reaction)
        if not transported:
            return ClassificationResult(
                is_member=False,
                explanation="No transported substrate beyond coupling spectators",
            )
        if not all(is_metal_ion(left_participant) or left_participant.chebi_id == "CHEBI:15378" for left_participant, _ in transported):
            return ClassificationResult(
                is_member=False,
                explanation="Non-cation co-transported substrate detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="ATPase-coupled monoatomic cation transmembrane transporter activity",
        )
