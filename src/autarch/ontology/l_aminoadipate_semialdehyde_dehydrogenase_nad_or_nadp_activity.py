"""L-aminoadipate-semialdehyde dehydrogenase [NAD(P)+] activity.

Catalysis of oxidation of L-aminoadipate-semialdehyde to L-2-aminoadipate with
NAD(+) or NADP(+) as acceptor.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_L_ALLYSINE = "CHEBI:58321"
CHEBI_L_2_AMINOADIPATE = "CHEBI:58672"


class LAminoadipateSemialdehydeDehydrogenaseNADOrNADPActivity(ReactionClass):
    """L-aminoadipate-semialdehyde dehydrogenase [NAD(P)+] activity."""

    GO_ID = "GO:0004043"
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
            if CHEBI_L_ALLYSINE not in left_ids or CHEBI_L_2_AMINOADIPATE not in right_ids:
                continue
            if CHEBI_H2O not in left_ids or right_ids.count(CHEBI_H_PLUS) != 2:
                continue
            if any(reduced in right_ids and oxidized in left_ids for oxidized, reduced in self.COFACTOR_PAIRS.items()):
                return ClassificationResult(
                    is_member=True,
                    explanation=f"L-aminoadipate-semialdehyde dehydrogenase [NAD(P)+] activity: oxidation of L-aminoadipate-semialdehyde to L-2-aminoadipate{suffix}",
                )
        return ClassificationResult(is_member=False, explanation="Requires L-aminoadipate-semialdehyde/L-2-aminoadipate with water and matched NAD(P) cofactor conversion")
