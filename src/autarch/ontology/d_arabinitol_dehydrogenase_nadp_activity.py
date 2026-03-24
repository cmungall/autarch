"""D-arabinitol dehydrogenase (NADP+) activity.

Catalysis of oxidation of D-arabinitol to D-ribulose or D-xylulose with NADP(+)
as acceptor.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_NADPH, CHEBI_NADP_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_D_ARABINITOL = "CHEBI:18333"
PRODUCT_IDS = {"CHEBI:17173", "CHEBI:17140"}


class DArabinitolDehydrogenaseNADPActivity(ReactionClass):
    """D-arabinitol dehydrogenase (NADP+) activity."""

    GO_ID = "GO:0033709"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            left_ids=[p.chebi_id for p in reaction.left_participants if p.chebi_id],
            right_ids=[p.chebi_id for p in reaction.right_participants if p.chebi_id],
            explanation="D-arabinitol dehydrogenase (NADP+) activity: NADP-linked oxidation of D-arabinitol to a pentulose",
        )
        if forward.is_member:
            return forward
        reverse = self._check_direction(
            left_ids=[p.chebi_id for p in reaction.right_participants if p.chebi_id],
            right_ids=[p.chebi_id for p in reaction.left_participants if p.chebi_id],
            explanation="D-arabinitol dehydrogenase (NADP+) activity: NADP-linked oxidation of D-arabinitol to a pentulose (reverse reaction orientation)",
        )
        if reverse.is_member:
            return reverse
        return forward

    def _check_direction(
        self,
        left_ids: list[str],
        right_ids: list[str],
        explanation: str,
    ) -> ClassificationResult:
        if CHEBI_D_ARABINITOL not in left_ids or not (PRODUCT_IDS & set(right_ids)):
            return ClassificationResult(is_member=False, explanation="Requires D-arabinitol substrate and a supported pentulose product")
        if CHEBI_NADP_PLUS not in left_ids or CHEBI_NADPH not in right_ids or CHEBI_H_PLUS not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires NADP(+)/NADPH coupling with hydron release")
        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {CHEBI_D_ARABINITOL, CHEBI_NADP_PLUS}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_NADPH, CHEBI_H_PLUS}]
        if len(left_core) != 0 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Unexpected side products for D-arabinitol dehydrogenase chemistry")
        if right_core[0] not in PRODUCT_IDS:
            return ClassificationResult(is_member=False, explanation="Product does not match a supported D-arabinitol oxidation product")
        return ClassificationResult(is_member=True, explanation=explanation)
