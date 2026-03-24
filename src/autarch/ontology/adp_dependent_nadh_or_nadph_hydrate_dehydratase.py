"""ADP-dependent NAD(P)H-hydrate dehydratase activity.

Catalysis of the reaction: NAD(P)HX + ADP = NAD(P)H + AMP + phosphate + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_AMP, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_NADHX = "CHEBI:64074"
CHEBI_NADPHX = "CHEBI:64076"
CHEBI_POLYPHOSPHATE = "CHEBI:16838"

SUBSTRATE_PRODUCT_PAIRS = {
    CHEBI_NADHX: CHEBI_NADH,
    CHEBI_NADPHX: CHEBI_NADPH,
}


class ADPDependentNADHOrNADPHHydrateDehydratase(ReactionClass):
    """ADP-dependent NAD(P)H-hydrate dehydratase activity.

    Catalysis of the reaction: NAD(P)HX + ADP = NAD(P)H + AMP + phosphate + H(+).
    """

    GO_ID = "GO:0052855"
    EC_NUMBER_PREFIX = "4.2.1.136"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebis = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        if CHEBI_ADP not in left_chebis or CHEBI_AMP not in right_chebis or CHEBI_POLYPHOSPHATE not in right_chebis:
            return ClassificationResult(is_member=False, explanation="Requires ADP consumption with AMP and phosphate release")

        substrate = next((chebi for chebi in SUBSTRATE_PRODUCT_PAIRS if chebi in left_chebis), None)
        if substrate is None:
            return ClassificationResult(is_member=False, explanation="No NAD(P)HX hydrate substrate detected")
        if SUBSTRATE_PRODUCT_PAIRS[substrate] not in right_chebis:
            return ClassificationResult(is_member=False, explanation="Missing the corresponding restored NAD(P)H product")

        substantive_left = left_chebis - {CHEBI_ADP, substrate}
        substantive_right = right_chebis - {CHEBI_AMP, CHEBI_POLYPHOSPHATE, SUBSTRATE_PRODUCT_PAIRS[substrate], CHEBI_H_PLUS}
        if substantive_left or substantive_right:
            return ClassificationResult(is_member=False, explanation="Contains additional substantive participants outside ADP-dependent NAD(P)HX repair chemistry")

        return ClassificationResult(
            is_member=True,
            explanation="ADP-dependent NAD(P)H-hydrate dehydratase: ADP-coupled dehydration of NAD(P)HX back to NAD(P)H",
        )
