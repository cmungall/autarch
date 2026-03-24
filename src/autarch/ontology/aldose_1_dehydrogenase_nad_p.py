"""aldose 1-dehydrogenase [NAD(P)(+)].

Catalysis of the reaction: aldopyranose + NAD(P)(+) = aldono-1,5-lactone +
NAD(P)H + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class Aldose1DehydrogenaseNADP(ReactionClass):
    """aldose 1-dehydrogenase [NAD(P)(+)].

    Catalysis of the reaction: aldopyranose + NAD(P)(+) = aldono-1,5-lactone +
    NAD(P)H + H(+).
    """

    EC_NUMBER_PREFIX = "1.1.1.359"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:53455": "CHEBI:15867",
        "CHEBI:4167": "CHEBI:16217",
        "CHEBI:140379": "CHEBI:140380",
        "CHEBI:46987": "CHEBI:17100",
        "CHEBI:4139": "CHEBI:15945",
    }
    PRODUCT_TO_SUBSTRATE = {product: substrate for substrate, product in SUBSTRATE_TO_PRODUCT.items()}
    COFACTOR_PAIRS = {
        CHEBI_NAD_PLUS: CHEBI_NADH,
        CHEBI_NADP_PLUS: CHEBI_NADPH,
    }
    REVERSE_COFACTOR_PAIRS = {value: key for key, value in COFACTOR_PAIRS.items()}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
            substrate_to_product=self.SUBSTRATE_TO_PRODUCT,
            cofactor_pairs=self.COFACTOR_PAIRS,
            explanation="Aldose 1-dehydrogenase [NAD(P)(+)]: oxidation of an aldopyranose to an aldono-1,5-lactone",
        )
        if forward.is_member:
            return forward

        reverse = self._check_direction(
            left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
            substrate_to_product=self.PRODUCT_TO_SUBSTRATE,
            cofactor_pairs=self.REVERSE_COFACTOR_PAIRS,
            explanation="Aldose 1-dehydrogenase [NAD(P)(+)]: reverse reduction of an aldono-1,5-lactone",
        )
        if reverse.is_member:
            return reverse

        return forward

    def _check_direction(
        self,
        left_ids: list[str],
        right_ids: list[str],
        substrate_to_product: dict[str, str],
        cofactor_pairs: dict[str, str],
        explanation: str,
    ) -> ClassificationResult:
        matched_pair = next(
            ((left_cofactor, right_cofactor) for left_cofactor, right_cofactor in cofactor_pairs.items() if left_cofactor in left_ids and right_cofactor in right_ids),
            None,
        )
        if matched_pair is None:
            return ClassificationResult(
                is_member=False,
                explanation="Requires a matching NAD(P) cofactor conversion",
            )

        left_core = [chebi_id for chebi_id in left_ids if chebi_id != matched_pair[0]]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {matched_pair[1], CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one aldose-like substrate and one lactone product",
            )

        substrate = left_core[0]
        product = right_core[0]
        if substrate_to_product.get(substrate) != product:
            return ClassificationResult(
                is_member=False,
                explanation="Organic branch does not match a supported aldose/lactone pair",
            )

        return ClassificationResult(is_member=True, explanation=explanation)
