"""6,7-dihydropteridine reductase activity.

Catalysis of reversible NAD(P)-linked interconversion between 6,7-dihydropteridines
and the corresponding tetrahydropteridines.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.gdp_4_dehydro_d_rhamnose_reductase import _check_reversible_redox_pair
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class SixSevenDihydropteridineReductase(ReactionClass):
    """6,7-dihydropteridine reductase activity.

    Catalysis of reversible NAD(P)-linked interconversion between 6,7-dihydropteridines
    and the corresponding tetrahydropteridines.
    """

    GO_ID = "GO:0004155"
    EC_NUMBER_PREFIX = "1.5.1.34"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    REDUCED_TO_OXIDIZED = {
        "CHEBI:28889": "CHEBI:30156",
        "CHEBI:59560": "CHEBI:43120",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_redox_pair(
            reaction,
            substrate_to_product=self.REDUCED_TO_OXIDIZED,
            cofactor_pairs={CHEBI_NADP_PLUS: CHEBI_NADPH, CHEBI_NAD_PLUS: CHEBI_NADH},
            explanation="6,7-dihydropteridine reductase: NAD(P)-linked interconversion of dihydropteridines and tetrahydropteridines",
        )
