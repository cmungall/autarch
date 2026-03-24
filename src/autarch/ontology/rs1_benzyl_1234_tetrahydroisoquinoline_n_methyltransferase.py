"""(RS)-1-benzyl-1,2,3,4-tetrahydroisoquinoline N-methyltransferase activity.

Catalysis of the reaction: (RS)-1-benzyl-1,2,3,4-tetrahydroisoquinoline + S-adenosyl-L-methionine = N-methyl-(RS)-1-benzyl-1,2,3,4-tetrahydroisoquinoline + S-adenosyl-L-homocysteine + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase

CHEBI_BENZYL_TETRAHYDROISOQUINOLINE = "CHEBI:57902"
CHEBI_N_METHYL_BENZYL_TETRAHYDROISOQUINOLINE = "CHEBI:57598"


class RS1Benzyl1234TetrahydroisoquinolineNMethyltransferase(Methyltransferase):
    """(RS)-1-benzyl-1,2,3,4-tetrahydroisoquinoline N-methyltransferase activity.

    Catalysis of the reaction: (RS)-1-benzyl-1,2,3,4-tetrahydroisoquinoline + S-adenosyl-L-methionine = N-methyl-(RS)-1-benzyl-1,2,3,4-tetrahydroisoquinoline + S-adenosyl-L-homocysteine + H(+).
    """

    GO_ID = "GO:0030776"
    EC_NUMBER_PREFIX = "2.1.1.115"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [p for p in reaction.left_participants if p.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}]
        right_core = [p for p in reaction.right_participants if p.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one tetrahydroisoquinoline substrate and one methylated product",
            )
        if (
            left_core[0].chebi_id != CHEBI_BENZYL_TETRAHYDROISOQUINOLINE
            or right_core[0].chebi_id != CHEBI_N_METHYL_BENZYL_TETRAHYDROISOQUINOLINE
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Requires N-methylation of benzyl tetrahydroisoquinoline",
            )

        return ClassificationResult(
            is_member=True,
            explanation="(RS)-1-benzyl-1,2,3,4-tetrahydroisoquinoline N-methyltransferase: SAM-dependent N-methylation of a tetrahydroisoquinoline alkaloid",
        )
