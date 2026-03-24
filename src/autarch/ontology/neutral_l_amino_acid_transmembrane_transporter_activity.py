"""neutral L-amino acid transmembrane transporter activity.

Enables the transfer of neutral L-amino acids from one side of a membrane to
the other.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.amino_acid_transmembrane_transporter_activity import (
    AminoAcidTransmembraneTransporterActivity,
)
from autarch.ontology.transport_utils import (
    coupled_transported_pairs,
    is_neutral_l_amino_acid,
)


class NeutralLAminoAcidTransmembraneTransporterActivity(
    AminoAcidTransmembraneTransporterActivity
):
    """neutral L-amino acid transmembrane transporter activity.

    Enables the transfer of neutral L-amino acids from one side of a membrane
    to the other.
    """

    GO_ID = "GO:0015175"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        if not reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Not a transport reaction",
            )

        transported = coupled_transported_pairs(reaction)
        substantive = [
            left_participant
            for left_participant, _ in transported
            if left_participant.chebi_id not in {"CHEBI:15378", "CHEBI:29101"}
        ]
        if not substantive:
            return ClassificationResult(
                is_member=False,
                explanation="No transported amino-acid substrate beyond coupling ions",
            )
        if not all(is_neutral_l_amino_acid(participant) for participant in substantive):
            return ClassificationResult(
                is_member=False,
                explanation="Transported substrate is not consistently a neutral L-amino acid",
            )
        return ClassificationResult(
            is_member=True,
            explanation="Neutral L-amino acid transmembrane transporter activity",
        )
