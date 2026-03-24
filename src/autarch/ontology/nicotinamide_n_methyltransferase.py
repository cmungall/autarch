"""nicotinamide N-methyltransferase activity.

Catalysis of the reaction: nicotinamide + S-adenosyl-L-methionine = 1-methylnicotinamide + S-adenosyl-L-homocysteine.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase

CHEBI_NICOTINAMIDE = "CHEBI:17154"
CHEBI_1_METHYLNICOTINAMIDE = "CHEBI:16797"


class NicotinamideNMethyltransferase(Methyltransferase):
    """nicotinamide N-methyltransferase activity.

    Catalysis of the reaction: nicotinamide + S-adenosyl-L-methionine = 1-methylnicotinamide + S-adenosyl-L-homocysteine.
    """

    GO_ID = "GO:0008112"
    EC_NUMBER_PREFIX = "2.1.1.1"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [p for p in reaction.left_participants if p.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}]
        right_core = [p for p in reaction.right_participants if p.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one nicotinamide substrate and one methylated product",
            )

        if left_core[0].chebi_id != CHEBI_NICOTINAMIDE or right_core[0].chebi_id != CHEBI_1_METHYLNICOTINAMIDE:
            return ClassificationResult(
                is_member=False,
                explanation="Requires nicotinamide methylation to 1-methylnicotinamide",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Nicotinamide N-methyltransferase: SAM-dependent N-methylation of nicotinamide",
        )
