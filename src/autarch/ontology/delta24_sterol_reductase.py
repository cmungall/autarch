"""Delta24-sterol reductase activity.

Catalysis of reversible NADP-linked interconversion between Delta(24)-reduced
sterols and their Delta(24)-unsaturated counterparts.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_NADPH, CHEBI_NADP_PLUS
from autarch.ontology.gdp_4_dehydro_d_rhamnose_reductase import _check_reversible_redox_pair
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class Delta24SterolReductase(ReactionClass):
    """Delta24-sterol reductase activity.

    Catalysis of reversible NADP-linked interconversion between Delta(24)-reduced
    sterols and their Delta(24)-unsaturated counterparts.
    """

    GO_ID = "GO:0050614"
    EC_NUMBER_PREFIX = "1.3.1.72"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    REDUCED_TO_OXIDIZED = {
        "CHEBI:134072": "CHEBI:147457",
        "CHEBI:16113": "CHEBI:15889",
        "CHEBI:15889": "CHEBI:147458",
        "CHEBI:16114": "CHEBI:15890",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_redox_pair(
            reaction,
            substrate_to_product=self.REDUCED_TO_OXIDIZED,
            cofactor_pairs={CHEBI_NADP_PLUS: CHEBI_NADPH},
            explanation="Delta24-sterol reductase: NADP-linked interconversion of Delta(24)-reduced and Delta(24)-unsaturated sterols",
        )
