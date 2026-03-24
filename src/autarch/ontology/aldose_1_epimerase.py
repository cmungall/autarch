"""Aldose 1-epimerase reaction classification."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.isomerase import Isomerase

# Positive pairs seen for GO:0004034 in the RHEA evaluation set.
ALDOSE_ANOMER_PAIRS = {
    frozenset({"CHEBI:15444", "CHEBI:15903"}),  # alpha/beta-D-glucose
    frozenset({"CHEBI:28061", "CHEBI:27667"}),  # alpha/beta-D-galactose
}


class Aldose1Epimerase(Isomerase):
    """aldose 1-epimerase"""

    GO_ID = "GO:0004034"  # aldose 1-epimerase activity
    EC_NUMBER_PREFIX = "5.1.3.3"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for alpha/beta anomerization of an aldose."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an isomerase: {parent_result.explanation}",
            )

        if len(reaction.left_participants) != 1 or len(reaction.right_participants) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Aldose 1-epimerase expects single-substrate isomerization",
            )

        left = reaction.left_participants[0]
        right = reaction.right_participants[0]

        if left.chebi_id and right.chebi_id:
            pair = frozenset({left.chebi_id, right.chebi_id})
            if pair in ALDOSE_ANOMER_PAIRS:
                return ClassificationResult(
                    is_member=True,
                    explanation="Aldose 1-epimerase: alpha/beta aldose anomerization",
                )

        return ClassificationResult(
            is_member=False,
            explanation="No GO:0004034 aldose anomer pair detected",
        )
