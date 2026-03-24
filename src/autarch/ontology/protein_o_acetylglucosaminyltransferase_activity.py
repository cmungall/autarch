"""protein O-acetylglucosaminyltransferase activity.

Catalysis of transfer of N-acetylglucosamine from UDP-GlcNAc onto seryl or
threonyl residues in protein substrates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_UDP
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_UDP_GLCNAC = "CHEBI:57705"
SUBSTRATE_TO_PRODUCT = {
    "CHEBI:29999": "CHEBI:90838",
    "CHEBI:30013": "CHEBI:90840",
}


class ProteinOAcetylglucosaminyltransferaseActivity(ReactionClass):
    """protein O-acetylglucosaminyltransferase activity."""

    GO_ID = "GO:0097363"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            left_ids=[p.chebi_id for p in reaction.left_participants if p.chebi_id],
            right_ids=[p.chebi_id for p in reaction.right_participants if p.chebi_id],
            explanation="protein O-acetylglucosaminyltransferase activity: transfer of GlcNAc from UDP-GlcNAc to a protein seryl/threonyl residue",
        )
        if forward.is_member:
            return forward
        reverse = self._check_direction(
            left_ids=[p.chebi_id for p in reaction.right_participants if p.chebi_id],
            right_ids=[p.chebi_id for p in reaction.left_participants if p.chebi_id],
            explanation="protein O-acetylglucosaminyltransferase activity: transfer of GlcNAc from UDP-GlcNAc to a protein seryl/threonyl residue (reverse reaction orientation)",
        )
        if reverse.is_member:
            return reverse
        return forward

    def _check_direction(self, left_ids: list[str], right_ids: list[str], explanation: str) -> ClassificationResult:
        if CHEBI_UDP_GLCNAC not in left_ids or CHEBI_UDP not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires UDP-GlcNAc donor and UDP coproduct")
        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_UDP_GLCNAC]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_UDP, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one protein acceptor and one O-GlcNAcylated protein product")
        if SUBSTRATE_TO_PRODUCT.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="Protein branch does not match a supported O-GlcNAc transfer pair")
        return ClassificationResult(is_member=True, explanation=explanation)
