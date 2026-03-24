"""10-hydroxydihydrosanguinarine 10-O-methyltransferase activity.

Catalysis of the reaction: 10-hydroxydihydrosanguinarine + S-adenosyl-L-methionine = dihydrochelirubine + S-adenosyl-L-homocysteine + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase

CHEBI_10_HYDROXYDIHYDROSANGUINARINE = "CHEBI:15878"
CHEBI_DIHYDROCHELIRUBINE = "CHEBI:17789"


class TenHydroxydihydrosanguinarine10OMethyltransferase(Methyltransferase):
    """10-hydroxydihydrosanguinarine 10-O-methyltransferase activity.

    Catalysis of the reaction: 10-hydroxydihydrosanguinarine + S-adenosyl-L-methionine = dihydrochelirubine + S-adenosyl-L-homocysteine + H(+).
    """

    GO_ID = "GO:0030779"
    EC_NUMBER_PREFIX = "2.1.1.119"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [p for p in reaction.left_participants if p.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}]
        right_core = [p for p in reaction.right_participants if p.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one alkaloid substrate and one methylated alkaloid product",
            )

        if (
            left_core[0].chebi_id != CHEBI_10_HYDROXYDIHYDROSANGUINARINE
            or right_core[0].chebi_id != CHEBI_DIHYDROCHELIRUBINE
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Requires 10-hydroxydihydrosanguinarine methylation to dihydrochelirubine",
            )

        return ClassificationResult(
            is_member=True,
            explanation="10-hydroxydihydrosanguinarine 10-O-methyltransferase: SAM-dependent O-methylation of the benzophenanthridine alkaloid substrate",
        )
