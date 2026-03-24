"""1-piperideine-2-carboxylate/1-pyrroline-2-carboxylate reductase[NAD(P)H].

Catalysis of reversible NAD(P)-linked interconversion between proline or
pipecolate substrates and the corresponding cyclic imino-acid products.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.gdp_4_dehydro_d_rhamnose_reductase import _check_reversible_redox_pair
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class OnePiperideine2CarboxylateOnePyrroline2CarboxylateReductaseNADHOrNADPH(ReactionClass):
    """1-piperideine-2-carboxylate/1-pyrroline-2-carboxylate reductase[NAD(P)H].

    Catalysis of reversible NAD(P)-linked interconversion between proline or
    pipecolate substrates and the corresponding cyclic imino-acid products.
    """

    EC_NUMBER_PREFIX = "1.5.1.1"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    REDUCED_TO_OXIDIZED = {
        "CHEBI:61185": "CHEBI:77631",
        "CHEBI:60039": "CHEBI:39785",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_redox_pair(
            reaction,
            substrate_to_product=self.REDUCED_TO_OXIDIZED,
            cofactor_pairs={CHEBI_NADP_PLUS: CHEBI_NADPH, CHEBI_NAD_PLUS: CHEBI_NADH},
            explanation="1-piperideine-2-carboxylate/1-pyrroline-2-carboxylate reductase[NAD(P)H]: NAD(P)-linked oxidation of proline or pipecolate substrates",
        )
