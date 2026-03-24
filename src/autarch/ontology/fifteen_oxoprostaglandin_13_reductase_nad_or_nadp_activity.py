"""15-oxoprostaglandin 13-reductase [NAD(P)+] activity.

Catalysis of oxidation of 13,14-dihydro-15-oxoprostaglandin E2 to
15-oxoprostaglandin E2 with NAD(+) or NADP(+) as acceptor.
"""

from autarch.ontology.gdp_4_dehydro_d_rhamnose_reductase import _check_reversible_redox_pair
from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class FifteenOxoprostaglandin13ReductaseNADOrNADPActivity(ReactionClass):
    """15-oxoprostaglandin 13-reductase [NAD(P)+] activity."""

    GO_ID = "GO:0047522"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_redox_pair(
            reaction,
            substrate_to_product={"CHEBI:57402": "CHEBI:57400"},
            cofactor_pairs={CHEBI_NADP_PLUS: CHEBI_NADPH, CHEBI_NAD_PLUS: CHEBI_NADH},
            explanation="15-oxoprostaglandin 13-reductase [NAD(P)+] activity: reversible 13,14-double-bond oxidation state interconversion of 15-oxoprostaglandin E2",
        )
