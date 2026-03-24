"""monocyclic monoterpene ketone monooxygenase.

Catalysis of NADPH-dependent monooxygenation of supported monocyclic
monoterpene ketones to the corresponding oxepanone products.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class MonocyclicMonoterpeneKetoneMonooxygenase(ReactionClass):
    """monocyclic monoterpene ketone monooxygenase.

    Catalysis of NADPH-dependent monooxygenation of supported monocyclic
    monoterpene ketones to the corresponding oxepanone products.
    """

    EC_NUMBER_PREFIX = "1.14.13.105"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCTS = {
        "CHEBI:15410": {"CHEBI:50250"},
        "CHEBI:23733": {"CHEBI:50238", "CHEBI:64229"},
        "CHEBI:154": {"CHEBI:64230"},
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        required_left = {CHEBI_NADPH, CHEBI_O2, CHEBI_H_PLUS}
        if not required_left.issubset(left_ids):
            return ClassificationResult(is_member=False, explanation="Requires NADPH, dioxygen, and hydron substrates")
        if CHEBI_NADP_PLUS not in right_ids or CHEBI_H2O not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires NADP(+) and water products")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in required_left]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_NADP_PLUS, CHEBI_H2O}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one monoterpene ketone substrate and one oxygenated product")

        if right_core[0] not in self.SUBSTRATE_TO_PRODUCTS.get(left_core[0], set()):
            return ClassificationResult(is_member=False, explanation="Organic branch does not match a supported monoterpene ketone monooxygenase pair")

        return ClassificationResult(is_member=True, explanation="Monocyclic monoterpene ketone monooxygenase: NADPH-dependent monooxygenation of a supported monoterpene ketone")
