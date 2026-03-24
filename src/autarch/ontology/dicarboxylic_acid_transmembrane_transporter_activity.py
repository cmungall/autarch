"""dicarboxylic acid transmembrane transporter activity.

Enables the transfer of dicarboxylic acids from one side of a membrane to the
other.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.carboxylic_acid_transmembrane_transporter_activity import (
    CarboxylicAcidTransmembraneTransporterActivity,
)
from autarch.ontology.transport_utils import (
    coupled_transported_pairs,
    is_dicarboxylic_acid,
)


class DicarboxylicAcidTransmembraneTransporterActivity(
    CarboxylicAcidTransmembraneTransporterActivity
):
    """dicarboxylic acid transmembrane transporter activity.

    Enables the transfer of dicarboxylic acids from one side of a membrane to
    the other.
    """

    GO_ID = "GO:0005310"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent = super().check_membership_impl(reaction)
        if not parent.is_member:
            return parent

        transported = coupled_transported_pairs(reaction)
        substantive = [
            left_participant
            for left_participant, _ in transported
            if left_participant.chebi_id not in {"CHEBI:15378", "CHEBI:29101"}
        ]
        if not substantive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive dicarboxylic-acid substrate detected",
            )
        if not all(is_dicarboxylic_acid(participant) for participant in substantive):
            return ClassificationResult(
                is_member=False,
                explanation="Transported substrate is not consistently dicarboxylic-acid-like",
            )
        return ClassificationResult(
            is_member=True,
            explanation="Dicarboxylic acid transmembrane transporter activity",
        )
