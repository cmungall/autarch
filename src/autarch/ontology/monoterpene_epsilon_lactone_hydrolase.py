"""monoterpene epsilon-lactone hydrolase.

Catalysis of the reaction: a monoterpene epsilon-lactone + H2O = the
corresponding hydroxyheptanoate + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class MonoterpeneEpsilonLactoneHydrolase(ReactionClass):
    """monoterpene epsilon-lactone hydrolase.

    Catalysis of the reaction: a monoterpene epsilon-lactone + H2O = the
    corresponding hydroxyheptanoate + H(+).
    """

    EC_NUMBER_PREFIX = "3.1.1.83"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:50238": "CHEBI:64224",
        "CHEBI:228": "CHEBI:64225",
        "CHEBI:233": "CHEBI:64226",
        "CHEBI:50250": "CHEBI:64221",
        "CHEBI:50243": "CHEBI:64223",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_H2O not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires water as reactant")
        if CHEBI_H_PLUS not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydron product")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_H2O]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_H_PLUS]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one lactone substrate and one hydroxyacid product")

        substrate = left_core[0]
        product = right_core[0]
        if self.SUBSTRATE_TO_PRODUCT.get(substrate) != product:
            return ClassificationResult(is_member=False, explanation="Organic branch does not match a supported monoterpene epsilon-lactone hydrolase pair")

        return ClassificationResult(is_member=True, explanation="Monoterpene epsilon-lactone hydrolase: water-dependent ring opening of a supported monoterpene epsilon-lactone")
