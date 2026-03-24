"""quinate/shikimate dehydrogenase [NAD(P)(+)].

Catalysis of reversible NAD(P)-linked interconversion between quinate or
shikimate substrates and the corresponding dehydro products.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.gdp_4_dehydro_d_rhamnose_reductase import _check_reversible_redox_pair
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class QuinateShikimateDehydrogenaseNADP(ReactionClass):
    """quinate/shikimate dehydrogenase [NAD(P)(+)].

    Catalysis of reversible NAD(P)-linked interconversion between quinate or
    shikimate substrates and the corresponding dehydro products.
    """

    EC_NUMBER_PREFIX = "1.1.1.282"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    REDUCED_TO_OXIDIZED = {
        "CHEBI:36208": "CHEBI:16630",
        "CHEBI:29751": "CHEBI:32364",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_redox_pair(
            reaction,
            substrate_to_product=self.REDUCED_TO_OXIDIZED,
            cofactor_pairs={CHEBI_NADP_PLUS: CHEBI_NADPH, CHEBI_NAD_PLUS: CHEBI_NADH},
            explanation="Quinate/shikimate dehydrogenase [NAD(P)(+)]: NAD(P)-linked oxidation of quinate or shikimate substrates",
        )
