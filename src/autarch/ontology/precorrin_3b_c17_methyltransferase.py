"""precorrin-3B C17-methyltransferase activity.

Catalysis of the reaction: precorrin-3B + S-adenosyl-L-methionine = precorrin-4 + S-adenosyl-L-homocysteine + 3 H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase

CHEBI_PRECORRIN_3B = "CHEBI:77870"
CHEBI_PRECORRIN_4 = "CHEBI:57769"


class Precorrin3BC17Methyltransferase(Methyltransferase):
    """precorrin-3B C17-methyltransferase activity.

    Catalysis of the reaction: precorrin-3B + S-adenosyl-L-methionine = precorrin-4 + S-adenosyl-L-homocysteine + 3 H(+).
    """

    GO_ID = "GO:0030789"
    EC_NUMBER_PREFIX = "2.1.1.131"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [p for p in reaction.left_participants if p.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}]
        right_core = [p for p in reaction.right_participants if p.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one precorrin substrate and one methylated corrin product",
            )

        if left_core[0].chebi_id != CHEBI_PRECORRIN_3B or right_core[0].chebi_id != CHEBI_PRECORRIN_4:
            return ClassificationResult(
                is_member=False,
                explanation="Requires precorrin-3B methylation to precorrin-4",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Precorrin-3B C17-methyltransferase: SAM-dependent corrin-ring methylation at the precorrin-3B to precorrin-4 step",
        )
