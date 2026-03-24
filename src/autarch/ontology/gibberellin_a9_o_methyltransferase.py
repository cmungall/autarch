"""gibberellin A9 O-methyltransferase.

Catalysis of SAM-dependent methyl ester formation for supported gibberellin
substrates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_SAH, CHEBI_SAM
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class GibberellinA9OMethyltransferase(ReactionClass):
    """gibberellin A9 O-methyltransferase.

    Catalysis of SAM-dependent methyl ester formation for supported gibberellin
    substrates.
    """

    EC_NUMBER_PREFIX = "2.1.1.275"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:73251": "CHEBI:73252",
        "CHEBI:58590": "CHEBI:73253",
        "CHEBI:73255": "CHEBI:73256",
        "CHEBI:58526": "CHEBI:73257",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_sam_transfer_pairs(
            reaction,
            substrate_to_product=self.SUBSTRATE_TO_PRODUCT,
            explanation="Gibberellin A9 O-methyltransferase: SAM-dependent methyl ester formation for a supported gibberellin substrate",
        )


def _check_sam_transfer_pairs(
    reaction: Reaction,
    substrate_to_product: dict[str, str],
    explanation: str,
) -> ClassificationResult:
    left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
    right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

    if CHEBI_SAM not in left_ids or CHEBI_SAH not in right_ids:
        return ClassificationResult(is_member=False, explanation="Requires SAM donor and SAH coproduct")

    left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_SAM]
    right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_SAH]
    if len(left_core) != 1 or len(right_core) != 1:
        return ClassificationResult(is_member=False, explanation="Expected one gibberellin substrate and one methyl ester product")

    if substrate_to_product.get(left_core[0]) != right_core[0]:
        return ClassificationResult(is_member=False, explanation="Gibberellin branch does not match a supported methyltransferase pair")

    return ClassificationResult(is_member=True, explanation=explanation)
