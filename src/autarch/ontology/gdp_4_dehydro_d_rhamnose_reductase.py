"""GDP-4-dehydro-D-rhamnose reductase activity.

Catalysis of reversible interconversion between GDP-4-dehydro-D-rhamnose and
GDP-6-deoxy-alpha-D-talose or GDP-alpha-D-rhamnose with NAD(P)-linked redox.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_GDP_4_DEHYDRO_D_RHAMNOSE = "CHEBI:57964"


class GDP4DehydroDRhamnoseReductase(ReactionClass):
    """GDP-4-dehydro-D-rhamnose reductase activity.

    Catalysis of reversible interconversion between GDP-4-dehydro-D-rhamnose and
    GDP-6-deoxy-alpha-D-talose or GDP-alpha-D-rhamnose with NAD(P)-linked redox.
    """

    GO_ID = "GO:0042356"
    EC_NUMBER_PREFIX = "1.1.1.187"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    REDUCED_TO_OXIDIZED = {
        "CHEBI:57638": CHEBI_GDP_4_DEHYDRO_D_RHAMNOSE,
        "CHEBI:58224": CHEBI_GDP_4_DEHYDRO_D_RHAMNOSE,
    }
    COFACTOR_PAIRS = {
        CHEBI_NADP_PLUS: CHEBI_NADPH,
        CHEBI_NAD_PLUS: CHEBI_NADH,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_redox_pair(
            reaction,
            substrate_to_product=self.REDUCED_TO_OXIDIZED,
            cofactor_pairs=self.COFACTOR_PAIRS,
            explanation="GDP-4-dehydro-D-rhamnose reductase: NAD(P)-linked interconversion of GDP-6-deoxy sugars and GDP-4-dehydro-D-rhamnose",
        )


def _check_reversible_redox_pair(
    reaction: Reaction,
    substrate_to_product: dict[str, str],
    cofactor_pairs: dict[str, str],
    explanation: str,
) -> ClassificationResult:
    forward = _check_redox_direction(
        left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
        right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
        substrate_to_product=substrate_to_product,
        cofactor_pairs=cofactor_pairs,
        explanation=explanation,
    )
    if forward.is_member:
        return forward

    reverse = _check_redox_direction(
        left_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
        right_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
        substrate_to_product={product: substrate for substrate, product in substrate_to_product.items()},
        cofactor_pairs={product: substrate for substrate, product in cofactor_pairs.items()},
        explanation=f"{explanation} (reverse reaction orientation)",
    )
    if reverse.is_member:
        return reverse
    return forward


def _check_redox_direction(
    left_ids: list[str],
    right_ids: list[str],
    substrate_to_product: dict[str, str],
    cofactor_pairs: dict[str, str],
    explanation: str,
) -> ClassificationResult:
    matched_pair = next(
        ((left, right) for left, right in cofactor_pairs.items() if left in left_ids and right in right_ids),
        None,
    )
    if matched_pair is None:
        return ClassificationResult(is_member=False, explanation="Requires matched NAD(P) cofactor conversion")

    left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {matched_pair[0], CHEBI_H_PLUS}]
    right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {matched_pair[1], CHEBI_H_PLUS}]
    if len(left_core) != 1 or len(right_core) != 1:
        return ClassificationResult(is_member=False, explanation="Expected one organic substrate and one product")

    substrate = left_core[0]
    product = right_core[0]
    if substrate_to_product.get(substrate) != product:
        return ClassificationResult(is_member=False, explanation="Organic branch does not match a supported redox pair")

    return ClassificationResult(is_member=True, explanation=explanation)
