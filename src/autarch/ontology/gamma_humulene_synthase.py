"""gamma-humulene synthase.

Catalysis of sesquiterpene cyclization of farnesyl diphosphate to supported
gamma-humulene synthase product branches with diphosphate release.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_FARNESYL_DIPHOSPHATE = "CHEBI:175763"
CHEBI_DIPHOSPHATE = "CHEBI:33019"


class GammaHumuleneSynthase(ReactionClass):
    """gamma-humulene synthase.

    Catalysis of sesquiterpene cyclization of farnesyl diphosphate to supported
    gamma-humulene synthase product branches with diphosphate release.
    """

    EC_NUMBER_PREFIX = "4.2.3.56"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    PRODUCTS = {
        "CHEBI:49290",
        "CHEBI:49231",
        "CHEBI:6530",
        "CHEBI:49210",
        "CHEBI:49224",
        "CHEBI:49214",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_FARNESYL_DIPHOSPHATE not in left_ids or CHEBI_DIPHOSPHATE not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires farnesyl diphosphate substrate and diphosphate product",
            )

        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_FARNESYL_DIPHOSPHATE]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_DIPHOSPHATE]
        if left_core:
            return ClassificationResult(
                is_member=False,
                explanation="Expected farnesyl diphosphate to be the only organic substrate",
            )
        if len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one sesquiterpene cyclization product",
            )
        if right_core[0] not in self.PRODUCTS:
            return ClassificationResult(
                is_member=False,
                explanation="Product branch does not match a supported gamma-humulene synthase product",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Gamma-humulene synthase: cyclization of farnesyl diphosphate to a supported sesquiterpene product",
        )
