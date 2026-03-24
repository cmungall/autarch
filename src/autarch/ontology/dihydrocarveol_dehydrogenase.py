"""dihydrocarveol dehydrogenase.

Catalysis of the reaction: dihydrocarveol + NAD(+) = dihydrocarvone + NADH + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class DihydrocarveolDehydrogenase(ReactionClass):
    """dihydrocarveol dehydrogenase.

    Catalysis of the reaction: dihydrocarveol + NAD(+) = dihydrocarvone + NADH + H(+).
    """

    EC_NUMBER_PREFIX = "1.1.1.296"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCTS = {
        "CHEBI:50215": {"CHEBI:23733"},
        "CHEBI:149": {"CHEBI:154"},
        "CHEBI:50235": {"CHEBI:168"},
        "CHEBI:158": {"CHEBI:168"},
        "CHEBI:152": {"CHEBI:154"},
        "CHEBI:153": {"CHEBI:155"},
        "CHEBI:150": {"CHEBI:155"},
        "CHEBI:50233": {"CHEBI:166"},
        "CHEBI:50232": {"CHEBI:166"},
    }
    PRODUCT_TO_SUBSTRATES = {
        "CHEBI:23733": {"CHEBI:50215"},
        "CHEBI:154": {"CHEBI:149", "CHEBI:152"},
        "CHEBI:168": {"CHEBI:50235", "CHEBI:158"},
        "CHEBI:155": {"CHEBI:153", "CHEBI:150"},
        "CHEBI:166": {"CHEBI:50233", "CHEBI:50232"},
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
            oxidation=True,
        )
        if forward.is_member:
            return forward

        reverse = self._check_direction(
            left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
            oxidation=False,
        )
        if reverse.is_member:
            return reverse

        return forward

    def _check_direction(self, left_ids: list[str], right_ids: list[str], oxidation: bool) -> ClassificationResult:
        if oxidation:
            if CHEBI_NAD_PLUS not in left_ids or CHEBI_NADH not in right_ids:
                return ClassificationResult(
                    is_member=False,
                    explanation="Requires NAD(+) reduction to NADH",
                )
            left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_NAD_PLUS]
            right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_NADH, CHEBI_H_PLUS}]
            matcher = self.SUBSTRATE_TO_PRODUCTS
            explanation = "Dihydrocarveol dehydrogenase: NAD-dependent oxidation of a dihydrocarveol stereoisomer"
        else:
            if CHEBI_NADH not in left_ids or CHEBI_NAD_PLUS not in right_ids:
                return ClassificationResult(
                    is_member=False,
                    explanation="Requires NADH oxidation to NAD(+) for reverse orientation",
                )
            left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_NADH]
            right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_NAD_PLUS, CHEBI_H_PLUS}]
            matcher = self.PRODUCT_TO_SUBSTRATES
            explanation = "Dihydrocarveol dehydrogenase: reverse reduction of a dihydrocarvone stereoisomer"

        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one organic substrate and one organic product",
            )

        substrate = left_core[0]
        product = right_core[0]
        if substrate not in matcher or product not in matcher[substrate]:
            return ClassificationResult(
                is_member=False,
                explanation="Organic branch does not match a supported dihydrocarveol/dihydrocarvone pair",
            )

        return ClassificationResult(is_member=True, explanation=explanation)
