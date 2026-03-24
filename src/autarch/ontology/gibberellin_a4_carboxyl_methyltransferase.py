"""gibberellin A4 carboxyl methyltransferase.

Catalysis of SAM-dependent methyl ester formation for supported gibberellin
substrates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.gibberellin_a9_o_methyltransferase import _check_sam_transfer_pairs
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class GibberellinA4CarboxylMethyltransferase(ReactionClass):
    """gibberellin A4 carboxyl methyltransferase.

    Catalysis of SAM-dependent methyl ester formation for supported gibberellin
    substrates.
    """

    EC_NUMBER_PREFIX = "2.1.1.276"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:73251": "CHEBI:73252",
        "CHEBI:58590": "CHEBI:73253",
        "CHEBI:73255": "CHEBI:73256",
        "CHEBI:73258": "CHEBI:73260",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_sam_transfer_pairs(
            reaction,
            substrate_to_product=self.SUBSTRATE_TO_PRODUCT,
            explanation="Gibberellin A4 carboxyl methyltransferase: SAM-dependent methyl ester formation for a supported gibberellin substrate",
        )
