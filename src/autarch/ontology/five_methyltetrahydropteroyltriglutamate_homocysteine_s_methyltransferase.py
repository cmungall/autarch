"""5-methyltetrahydropteroyltriglutamate-homocysteine S-methyltransferase activity.

Catalysis of the reaction: 5-methyltetrahydropteroyltri-L-glutamate + L-homocysteine = tetrahydropteroyltri-L-glutamate + L-methionine.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_METHYLTHF_TRIGLU = "CHEBI:58207"
CHEBI_THF_TRIGLU = "CHEBI:58140"
CHEBI_HOMOCYSTEINE = "CHEBI:58199"
CHEBI_METHIONINE = "CHEBI:57844"


class FiveMethyltetrahydropteroyltriglutamateHomocysteineSMethyltransferase(ReactionClass):
    """5-methyltetrahydropteroyltriglutamate-homocysteine S-methyltransferase activity.

    Catalysis of the reaction: 5-methyltetrahydropteroyltri-L-glutamate + L-homocysteine = tetrahydropteroyltri-L-glutamate + L-methionine.
    """

    GO_ID = "GO:0003871"
    EC_NUMBER_PREFIX = "2.1.1.14"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebis = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        if left_chebis != {CHEBI_METHYLTHF_TRIGLU, CHEBI_HOMOCYSTEINE}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires methylated folate donor and homocysteine as the complete substrate set",
            )
        if right_chebis != {CHEBI_THF_TRIGLU, CHEBI_METHIONINE}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires demethylated folate and methionine as the complete product set",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Methyltetrahydropteroyltriglutamate-homocysteine S-methyltransferase: folate-mediated sulfur methyl transfer to homocysteine",
        )
