"""solute:sodium symporter activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.solute_monoatomic_cation_symporter_activity import (
    SoluteMonoatomicCationSymporterActivity,
)
from autarch.ontology.transport_utils import coupled_transported_pairs, is_sodium_ion


class SoluteSodiumSymporterActivity(SoluteMonoatomicCationSymporterActivity):
    """solute:sodium symporter activity."""

    GO_ID = "GO:0015370"

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
            explanation="Solute:sodium symporter activity",
        )
