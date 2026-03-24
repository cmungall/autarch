"""acyl-CoA hydrolase activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_COA
from autarch.ontology.thiolester_hydrolase import ThiolesterHydrolase


class AcylCoAHydrolase(ThiolesterHydrolase):
    """acyl-CoA hydrolase activity."""

    GO_ID = "GO:0016289"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Restrict thiolester hydrolysis to CoA-releasing branches."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_has_coa = any(participant.chebi_id == CHEBI_COA for participant in reaction.left_participants)
        right_has_coa = any(participant.chebi_id == CHEBI_COA for participant in reaction.right_participants)
        if not (left_has_coa ^ right_has_coa):
            return ClassificationResult(
                is_member=False,
                explanation="Requires CoA to appear on exactly one side of the hydrolysis",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Acyl-CoA hydrolase: thiolester hydrolysis that releases coenzyme A",
        )
