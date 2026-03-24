"""solute:proton symporter activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.symporter_activity import SymporterActivity
from autarch.ontology.transport_utils import coupled_transported_pairs


class SoluteProtonSymporterActivity(SymporterActivity):
    """solute:proton symporter activity."""

    GO_ID = "GO:0015295"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        transported = coupled_transported_pairs(reaction)
        has_proton = any(left_participant.chebi_id == CHEBI_H_PLUS for left_participant, _ in transported)
        has_other_solute = any(left_participant.chebi_id != CHEBI_H_PLUS for left_participant, _ in transported)
        if not (has_proton and has_other_solute):
            return ClassificationResult(
                is_member=False,
                explanation="No coupled proton-plus-solute symport detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Solute:proton symporter activity",
        )
