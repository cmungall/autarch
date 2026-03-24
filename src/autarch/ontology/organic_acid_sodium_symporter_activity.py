"""organic acid:sodium symporter activity.

Enables cotransport of an organic acid with sodium across a membrane.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE
from autarch.ontology.solute_sodium_symporter_activity import SoluteSodiumSymporterActivity
from autarch.ontology.transport_utils import (
    coupled_transported_pairs,
    is_carboxylic_acid_or_derivative,
)


class OrganicAcidSodiumSymporterActivity(SoluteSodiumSymporterActivity):
    """organic acid:sodium symporter activity.

    Enables cotransport of an organic acid with sodium across a membrane.
    """

    GO_ID = "GO:0005343"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent = super().check_membership_impl(reaction)
        if not parent.is_member:
            return parent

        transported = coupled_transported_pairs(reaction)
        organic_acids = [
            left_participant
            for left_participant, _ in transported
            if left_participant.chebi_id != "CHEBI:29101"
        ]
        if not organic_acids:
            return ClassificationResult(
                is_member=False,
                explanation="No organic-acid substrate accompanies sodium transport",
            )
        if not all(is_carboxylic_acid_or_derivative(participant) for participant in organic_acids):
            return ClassificationResult(
                is_member=False,
                explanation="Co-transported substrate is not consistently an organic acid",
            )
        return ClassificationResult(
            is_member=True,
            explanation="Organic acid:sodium symporter activity",
        )
