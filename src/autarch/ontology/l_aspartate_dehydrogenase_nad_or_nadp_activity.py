"""L-aspartate dehydrogenase [NAD(P)+] activity.

Catalysis of oxidative deamination of L-aspartate to oxaloacetate with
NAD(+) or NADP(+) as acceptor.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS, CHEBI_NH4
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_L_ASPARTATE = "CHEBI:29991"
CHEBI_OXALOACETATE = "CHEBI:16452"


class LAspartateDehydrogenaseNADOrNADPActivity(ReactionClass):
    """L-aspartate dehydrogenase [NAD(P)+] activity."""

    GO_ID = "GO:0033735"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    COFACTOR_PAIRS = {
        CHEBI_NADP_PLUS: CHEBI_NADPH,
        CHEBI_NAD_PLUS: CHEBI_NADH,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        for left_ids, right_ids, suffix in [
            ([p.chebi_id for p in reaction.left_participants if p.chebi_id], [p.chebi_id for p in reaction.right_participants if p.chebi_id], ""),
            ([p.chebi_id for p in reaction.right_participants if p.chebi_id], [p.chebi_id for p in reaction.left_participants if p.chebi_id], " (reverse reaction orientation)"),
        ]:
            if CHEBI_L_ASPARTATE not in left_ids or CHEBI_OXALOACETATE not in right_ids:
                continue
            if CHEBI_H2O not in left_ids or CHEBI_NH4 not in right_ids or CHEBI_H_PLUS not in right_ids:
                continue
            if any(oxidized in left_ids and reduced in right_ids for oxidized, reduced in self.COFACTOR_PAIRS.items()):
                return ClassificationResult(
                    is_member=True,
                    explanation=f"L-aspartate dehydrogenase [NAD(P)+] activity: oxidative deamination of L-aspartate to oxaloacetate{suffix}",
                )
        return ClassificationResult(is_member=False, explanation="Requires L-aspartate/oxaloacetate with water, ammonium, and matched NAD(P) cofactor conversion")
