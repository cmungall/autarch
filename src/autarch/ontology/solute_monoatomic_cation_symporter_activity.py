"""solute:monoatomic cation symporter activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.symporter_activity import SymporterActivity
from autarch.ontology.transport_utils import coupled_transported_pairs, is_inorganic_cation


class SoluteMonoatomicCationSymporterActivity(SymporterActivity):
    """solute:monoatomic cation symporter activity."""

    GO_ID = "GO:0015294"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        transported = coupled_transported_pairs(reaction)
        has_cation = any(is_inorganic_cation(left_participant) for left_participant, _ in transported)
        has_other_solute = any(not is_inorganic_cation(left_participant) for left_participant, _ in transported)
        if not (has_cation and has_other_solute):
            return ClassificationResult(
                is_member=False,
                explanation="No coupled cation-plus-solute symport detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Solute:monoatomic cation symporter activity",
        )
