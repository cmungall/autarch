"""tyramine N-methyltransferase activity.

Catalysis of the reaction: tyramine + S-adenosyl-L-methionine = N-methyltyramine + S-adenosyl-L-homocysteine + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase

CHEBI_TYRAMINE = "CHEBI:327995"
CHEBI_N_METHYLTYRAMINE = "CHEBI:58155"


class TyramineNMethyltransferase(Methyltransferase):
    """tyramine N-methyltransferase activity.

    Catalysis of the reaction: tyramine + S-adenosyl-L-methionine = N-methyltyramine + S-adenosyl-L-homocysteine + H(+).
    """

    GO_ID = "GO:0030738"
    EC_NUMBER_PREFIX = "2.1.1.27"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [p for p in reaction.left_participants if p.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}]
        right_core = [p for p in reaction.right_participants if p.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one tyramine substrate and one N-methylated amine product",
            )

        if left_core[0].chebi_id != CHEBI_TYRAMINE or right_core[0].chebi_id != CHEBI_N_METHYLTYRAMINE:
            return ClassificationResult(
                is_member=False,
                explanation="Requires N-methylation of tyramine to N-methyltyramine",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Tyramine N-methyltransferase: SAM-dependent N-methylation of tyramine",
        )
