"""rRNA small subunit pseudouridine methyltransferase Nep1.

Catalysis of SAM-dependent methylation of pseudouridine residues in small-subunit
rRNA to N(1)-methylpseudouridine.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.ontology.trna_cytidine_5_methyltransferase import _check_reversible_sam_exact_pair

CHEBI_PSEUDOURIDINE_IN_RRNA = "CHEBI:65314"
CHEBI_METHYLPSEUDOURIDINE_IN_RRNA = "CHEBI:74890"


class RRNASmallSubunitPseudouridineMethyltransferaseNep1(ReactionClass):
    """rRNA small subunit pseudouridine methyltransferase Nep1.

    Catalysis of SAM-dependent methylation of pseudouridine residues in small-subunit
    rRNA to N(1)-methylpseudouridine.
    """

    EC_NUMBER_PREFIX = "2.1.1.260"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_sam_exact_pair(
            reaction,
            substrate=CHEBI_PSEUDOURIDINE_IN_RRNA,
            product=CHEBI_METHYLPSEUDOURIDINE_IN_RRNA,
            explanation="rRNA small subunit pseudouridine methyltransferase Nep1: SAM-dependent methylation of pseudouridine in small-subunit rRNA",
        )
