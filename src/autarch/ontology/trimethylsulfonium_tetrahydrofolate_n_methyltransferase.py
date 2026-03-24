"""trimethylsulfonium-tetrahydrofolate N-methyltransferase activity.

Catalysis of the reaction: trimethylsulfonium + (6S)-5,6,7,8-tetrahydrofolate = dimethyl sulfide + (6S)-5-methyl-5,6,7,8-tetrahydrofolate + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_TRIMETHYLSULFONIUM = "CHEBI:17434"
CHEBI_DIMETHYL_SULFIDE = "CHEBI:17437"
CHEBI_THF = "CHEBI:57453"
CHEBI_METHYL_THF = "CHEBI:18608"


class TrimethylsulfoniumTetrahydrofolateNMethyltransferase(ReactionClass):
    """trimethylsulfonium-tetrahydrofolate N-methyltransferase activity.

    Catalysis of the reaction: trimethylsulfonium + (6S)-5,6,7,8-tetrahydrofolate = dimethyl sulfide + (6S)-5-methyl-5,6,7,8-tetrahydrofolate + H(+).
    """

    GO_ID = "GO:0047147"
    EC_NUMBER_PREFIX = "2.1.1.19"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebis = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        if left_chebis != {CHEBI_TRIMETHYLSULFONIUM, CHEBI_THF}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires trimethylsulfonium and tetrahydrofolate as the complete substrate set",
            )
        if right_chebis != {CHEBI_DIMETHYL_SULFIDE, CHEBI_METHYL_THF, CHEBI_H_PLUS}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires dimethyl sulfide, methyl-tetrahydrofolate, and hydron as the complete product set",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Trimethylsulfonium-tetrahydrofolate N-methyltransferase: folate-mediated methyl transfer from trimethylsulfonium",
        )
