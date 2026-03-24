"""glycerol dehydratase activity.

Catalysis of the reaction: glycerol = 3-hydroxypropanal + H2O.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_GLYCEROL = "CHEBI:17522"
CHEBI_HYDROXYPROPANAL = "CHEBI:17871"


class GlycerolDehydratase(ReactionClass):
    """glycerol dehydratase activity.

    Catalysis of the reaction: glycerol = 3-hydroxypropanal + H2O.
    """

    GO_ID = "GO:0046405"
    EC_NUMBER_PREFIX = "4.2.1.30"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebis = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        if left_chebis != {CHEBI_GLYCEROL}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires glycerol as the sole substantive substrate",
            )
        if right_chebis != {CHEBI_HYDROXYPROPANAL, CHEBI_H2O}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires 3-hydroxypropanal and water as the complete product set",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Glycerol dehydratase: dehydration of glycerol to 3-hydroxypropanal",
        )
