"""dihydrolipoyl dehydrogenase (NADH).

Catalysis of reversible NAD-linked interconversion between supported lipoyl-lysyl
protein redox states.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_NADH, CHEBI_NAD_PLUS
from autarch.ontology.gdp_4_dehydro_d_rhamnose_reductase import _check_reversible_redox_pair
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_DIHYDROLIPOYL_PROTEIN = "CHEBI:83100"
CHEBI_LIPOYL_PROTEIN = "CHEBI:83099"


class DihydrolipoylDehydrogenaseNADH(ReactionClass):
    """dihydrolipoyl dehydrogenase (NADH).

    Catalysis of reversible NAD-linked interconversion between supported
    lipoyl-lysyl protein redox states.

    The GO benchmark for this term currently includes neighboring NAD-linked
    oxidoreductases, so this classifier intentionally stays chemistry-first
    rather than broadening into those unrelated branches.
    """

    GO_ID = "GO:0004148"
    EC_BROAD_XREFS = ["1.4.1.27", "1.2.1.25", "1.2.1.105", "1.2.1.104"]
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_redox_pair(
            reaction,
            substrate_to_product={CHEBI_DIHYDROLIPOYL_PROTEIN: CHEBI_LIPOYL_PROTEIN},
            cofactor_pairs={CHEBI_NAD_PLUS: CHEBI_NADH},
            explanation="Dihydrolipoyl dehydrogenase (NADH): NAD-linked interconversion of lipoyl-lysyl protein redox states",
        )
