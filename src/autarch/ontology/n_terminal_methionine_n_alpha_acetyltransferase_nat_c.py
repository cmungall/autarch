"""N-terminal methionine N(alpha)-acetyltransferase NatC.

Catalysis of acetyl transfer from acetyl-CoA to selected methionyl protein
N termini.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.n_terminal_amino_acid_n_alpha_acetyltransferase_nat_a import _check_n_terminal_acetyltransferase
from autarch.ontology.reaction import IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE, ReactionClass


class NTerminalMethionineNAlphaAcetyltransferaseNatC(ReactionClass):
    """N-terminal methionine N(alpha)-acetyltransferase NatC.

    Catalysis of acetyl transfer from acetyl-CoA to selected methionyl protein
    N termini.
    """

    EC_NUMBER_PREFIX = "2.3.1.256"
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:133377": "CHEBI:133378",
        "CHEBI:133387": "CHEBI:133386",
        "CHEBI:133383": "CHEBI:133382",
        "CHEBI:133380": "CHEBI:133379",
        "CHEBI:133385": "CHEBI:133384",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_n_terminal_acetyltransferase(reaction, self.SUBSTRATE_TO_PRODUCT, "NatC")
