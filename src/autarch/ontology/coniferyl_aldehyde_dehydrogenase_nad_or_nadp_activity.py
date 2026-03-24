"""coniferyl-aldehyde dehydrogenase [NAD(P)+] activity.

Catalysis of oxidation of coniferyl aldehyde to ferulate with NAD(+) or NADP(+)
as acceptor.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_CONIFERYL_ALDEHYDE = "CHEBI:16547"
CHEBI_FERULATE = "CHEBI:29749"


class ConiferylAldehydeDehydrogenaseNADOrNADPActivity(ReactionClass):
    """coniferyl-aldehyde dehydrogenase [NAD(P)+] activity."""

    GO_ID = "GO:0050269"
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
            if CHEBI_CONIFERYL_ALDEHYDE not in left_ids or CHEBI_FERULATE not in right_ids:
                continue
            if CHEBI_H2O not in left_ids or right_ids.count(CHEBI_H_PLUS) != 2:
                continue
            if any(oxidized in left_ids and reduced in right_ids for oxidized, reduced in self.COFACTOR_PAIRS.items()):
                return ClassificationResult(
                    is_member=True,
                    explanation=f"coniferyl-aldehyde dehydrogenase [NAD(P)+] activity: oxidation of coniferyl aldehyde to ferulate{suffix}",
                )
        return ClassificationResult(is_member=False, explanation="Requires coniferyl aldehyde/ferulate with water and matched NAD(P) cofactor conversion")
