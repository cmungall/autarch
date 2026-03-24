"""methionine S-methyltransferase activity.

Catalysis of the reaction: L-methionine + S-adenosyl-L-methionine = S-methyl-L-methionine + S-adenosyl-L-homocysteine.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase

CHEBI_METHIONINE = "CHEBI:57844"
CHEBI_S_METHYL_METHIONINE = "CHEBI:58252"


class MethionineSMethyltransferase(Methyltransferase):
    """methionine S-methyltransferase activity.

    Catalysis of the reaction: L-methionine + S-adenosyl-L-methionine = S-methyl-L-methionine + S-adenosyl-L-homocysteine.
    """

    GO_ID = "GO:0030732"
    EC_NUMBER_PREFIX = "2.1.1.12"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [p for p in reaction.left_participants if p.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}]
        right_core = [p for p in reaction.right_participants if p.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one methionine substrate and one methylated sulfur product",
            )

        if left_core[0].chebi_id != CHEBI_METHIONINE or right_core[0].chebi_id != CHEBI_S_METHYL_METHIONINE:
            return ClassificationResult(
                is_member=False,
                explanation="Requires sulfur methylation of L-methionine to S-methyl-L-methionine",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Methionine S-methyltransferase: SAM-dependent sulfur methylation of methionine",
        )
