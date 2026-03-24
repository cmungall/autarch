"""N-terminal methionine N(alpha)-acetyltransferase NatB.

Catalysis of acetyl transfer from acetyl-CoA to selected methionyl protein N termini.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.n_terminal_amino_acid_n_alpha_acetyltransferase_nat_a import _check_n_terminal_acetyltransferase
from autarch.ontology.reaction import IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE, ReactionClass


class NTerminalMethionineNAlphaAcetyltransferaseNatB(ReactionClass):
    """N-terminal methionine N(alpha)-acetyltransferase NatB.

    Catalysis of acetyl transfer from acetyl-CoA to selected methionyl protein N termini.
    """

    EC_NUMBER_PREFIX = "2.3.1.254"
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:133045": "CHEBI:133063",
        "CHEBI:133356": "CHEBI:133358",
        "CHEBI:133359": "CHEBI:133360",
        "CHEBI:133361": "CHEBI:133362",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_n_terminal_acetyltransferase(reaction, self.SUBSTRATE_TO_PRODUCT, "NatB")
