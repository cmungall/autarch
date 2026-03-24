"""methylxanthine N(1)-demethylase.

Catalysis of N(1)-demethylation of methylxanthines with dioxygen and NAD(P)H,
producing formaldehyde.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_FORMALDEHYDE = "CHEBI:16842"


class MethylxanthineN1Demethylase(ReactionClass):
    """methylxanthine N(1)-demethylase.

    Catalysis of N(1)-demethylation of methylxanthines with dioxygen and NAD(P)H,
    producing formaldehyde.
    """

    EC_NUMBER_PREFIX = "1.14.13.178"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:25858": "CHEBI:48991",
        "CHEBI:27732": "CHEBI:28946",
        "CHEBI:28177": "CHEBI:62208",
    }
    COFACTOR_PAIRS = {
        CHEBI_NADH: CHEBI_NAD_PLUS,
        CHEBI_NADPH: CHEBI_NADP_PLUS,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_O2 not in left_ids or CHEBI_H_PLUS not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires dioxygen and hydron substrates")
        if CHEBI_FORMALDEHYDE not in right_ids or CHEBI_H2O not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires formaldehyde and water products")

        matched_pair = next(((left, right) for left, right in self.COFACTOR_PAIRS.items() if left in left_ids and right in right_ids), None)
        if matched_pair is None:
            return ClassificationResult(is_member=False, explanation="Requires NAD(P)H to NAD(P)(+) cofactor conversion")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {matched_pair[0], CHEBI_O2, CHEBI_H_PLUS}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {matched_pair[1], CHEBI_FORMALDEHYDE, CHEBI_H2O}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one methylxanthine substrate and one demethylated xanthine product")

        substrate = left_core[0]
        product = right_core[0]
        if self.SUBSTRATE_TO_PRODUCT.get(substrate) != product:
            return ClassificationResult(is_member=False, explanation="Organic branch does not match a supported methylxanthine N(1)-demethylase pair")

        return ClassificationResult(is_member=True, explanation="Methylxanthine N(1)-demethylase: NAD(P)H-dependent oxidative demethylation of a methylxanthine")
