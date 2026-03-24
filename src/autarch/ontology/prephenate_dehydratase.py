"""prephenate dehydratase activity.

Catalysis of the reaction: prephenate + H(+) = 3-phenylpyruvate + CO2 + H2O.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_CO2, CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_PREPHENATE = "CHEBI:29934"
CHEBI_PHENYLPYRUVATE = "CHEBI:18005"


class PrephenateDehydratase(ReactionClass):
    """prephenate dehydratase activity.

    Catalysis of the reaction: prephenate + H(+) = 3-phenylpyruvate + CO2 + H2O.
    """

    GO_ID = "GO:0004664"
    EC_NUMBER_PREFIX = "4.2.1.51"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebis = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        if left_chebis != {CHEBI_PREPHENATE, CHEBI_H_PLUS}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires prephenate and hydron as the complete substrate set",
            )
        if right_chebis != {CHEBI_PHENYLPYRUVATE, CHEBI_CO2, CHEBI_H2O}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires phenylpyruvate, carbon dioxide, and water as the complete product set",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Prephenate dehydratase: dehydration and decarboxylation of prephenate to phenylpyruvate",
        )
