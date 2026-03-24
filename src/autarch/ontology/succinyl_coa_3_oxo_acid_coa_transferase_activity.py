"""succinyl-CoA:3-oxo-acid CoA-transferase activity.

Catalysis of CoA transfer from succinyl-CoA to a 3-oxo acid acceptor.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_SUCCINYL_COA = "CHEBI:57292"
CHEBI_SUCCINATE = "CHEBI:30031"


class SuccinylCoA3OxoAcidCoATransferaseActivity(ReactionClass):
    """succinyl-CoA:3-oxo-acid CoA-transferase activity."""

    GO_ID = "GO:0008260"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:13705": "CHEBI:57286",
        "CHEBI:35973": "CHEBI:90726",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            left_ids=[p.chebi_id for p in reaction.left_participants if p.chebi_id],
            right_ids=[p.chebi_id for p in reaction.right_participants if p.chebi_id],
            explanation="succinyl-CoA:3-oxo-acid CoA-transferase activity: CoA transfer from succinyl-CoA to a 3-oxo acid acceptor",
        )
        if forward.is_member:
            return forward
        reverse = self._check_direction(
            left_ids=[p.chebi_id for p in reaction.right_participants if p.chebi_id],
            right_ids=[p.chebi_id for p in reaction.left_participants if p.chebi_id],
            explanation="succinyl-CoA:3-oxo-acid CoA-transferase activity: CoA transfer from succinyl-CoA to a 3-oxo acid acceptor (reverse reaction orientation)",
        )
        if reverse.is_member:
            return reverse
        return forward

    def _check_direction(self, left_ids: list[str], right_ids: list[str], explanation: str) -> ClassificationResult:
        if CHEBI_SUCCINYL_COA not in left_ids or CHEBI_SUCCINATE not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires succinyl-CoA on the substrate side and succinate on the product side")
        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_SUCCINYL_COA]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_SUCCINATE]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one 3-oxo acid acceptor and one 3-oxoacyl-CoA product")
        if self.SUBSTRATE_TO_PRODUCT.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="3-oxo acid branch does not match a supported succinyl-CoA transfer pair")
        return ClassificationResult(is_member=True, explanation=explanation)
