"""estradiol 17-beta-dehydrogenase [NAD(P)+] activity.

Catalysis of reversible interconversion between 17beta-estradiol and estrone
with NAD(+) or NADP(+) as cofactor.
"""

from autarch.ontology.gdp_4_dehydro_d_rhamnose_reductase import _check_reversible_redox_pair
from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class Estradiol17BetaDehydrogenaseNADOrNADPActivity(ReactionClass):
    """estradiol 17-beta-dehydrogenase [NAD(P)+] activity."""

    GO_ID = "GO:0004303"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_redox_pair(
            reaction,
            substrate_to_product={"CHEBI:16469": "CHEBI:17263"},
            cofactor_pairs={CHEBI_NADP_PLUS: CHEBI_NADPH, CHEBI_NAD_PLUS: CHEBI_NADH},
            explanation="estradiol 17-beta-dehydrogenase [NAD(P)+] activity: reversible estradiol/estrone interconversion with NAD(P) coupling",
        )
