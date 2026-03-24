"""N-terminal amino-acid N(alpha)-acetyltransferase NatA.

Catalysis of acetyl transfer from acetyl-CoA to the N terminus of selected
protein substrates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_ACETYL_COA = "CHEBI:57288"
CHEBI_COA = "CHEBI:57287"


class NTerminalAminoAcidNAlphaAcetyltransferaseNatA(ReactionClass):
    """N-terminal amino-acid N(alpha)-acetyltransferase NatA.

    Catalysis of acetyl transfer from acetyl-CoA to the N terminus of selected
    protein substrates.
    """

    EC_NUMBER_PREFIX = "2.3.1.255"
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:64739": "CHEBI:133375",
        "CHEBI:65250": "CHEBI:133372",
        "CHEBI:64718": "CHEBI:83683",
        "CHEBI:64723": "CHEBI:133369",
        "CHEBI:64738": "CHEBI:83690",
        "CHEBI:64741": "CHEBI:133371",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_n_terminal_acetyltransferase(
            reaction,
            self.SUBSTRATE_TO_PRODUCT,
            "NatA",
        )


def _check_n_terminal_acetyltransferase(
    reaction: Reaction,
    substrate_to_product: dict[str, str],
    label: str,
) -> ClassificationResult:
    left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
    right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

    forward = _check_n_terminal_direction(left_ids, right_ids, substrate_to_product, label)
    if forward.is_member:
        return forward

    reverse = _check_n_terminal_direction(right_ids, left_ids, substrate_to_product, label)
    if reverse.is_member:
        return ClassificationResult(is_member=True, explanation=f"{reverse.explanation} (reverse reaction orientation)")

    return forward


def _check_n_terminal_direction(left_ids: list[str], right_ids: list[str], substrate_to_product: dict[str, str], label: str) -> ClassificationResult:
    if CHEBI_ACETYL_COA not in left_ids or CHEBI_COA not in right_ids:
        return ClassificationResult(is_member=False, explanation="Requires acetyl-CoA donor and CoA product")

    left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_ACETYL_COA]
    right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_COA, CHEBI_H_PLUS}]
    if len(left_core) != 1 or len(right_core) != 1:
        return ClassificationResult(is_member=False, explanation="Expected one protein N-terminus substrate and one acetylated protein product")

    substrate = left_core[0]
    product = right_core[0]
    if substrate_to_product.get(substrate) != product:
        return ClassificationResult(is_member=False, explanation=f"Protein branch does not match a supported {label} substrate/product pair")

    return ClassificationResult(is_member=True, explanation=f"N-terminal amino-acid N(alpha)-acetyltransferase {label}: acetyl transfer to a supported protein N terminus")
