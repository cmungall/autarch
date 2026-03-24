"""nitroquinoline-N-oxide reductase [NAD(P)H] activity.

Catalysis of reduction of nitroquinoline N-oxide to the corresponding
hydroxyamino compound with NADH or NADPH.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_HYDROXYAMINO_QUINOLINE_N_OXIDE = "CHEBI:28469"
CHEBI_NITROQUINOLINE_N_OXIDE = "CHEBI:16907"


class NitroquinolineNOxideReductaseNADHOrNADPHActivity(ReactionClass):
    """nitroquinoline-N-oxide reductase [NAD(P)H] activity."""

    GO_ID = "GO:0050465"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    COFACTOR_PAIRS = {
        CHEBI_NADPH: CHEBI_NADP_PLUS,
        CHEBI_NADH: CHEBI_NAD_PLUS,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        for left_ids, right_ids, suffix in [
            ([p.chebi_id for p in reaction.left_participants if p.chebi_id], [p.chebi_id for p in reaction.right_participants if p.chebi_id], " (reverse reaction orientation)"),
            ([p.chebi_id for p in reaction.right_participants if p.chebi_id], [p.chebi_id for p in reaction.left_participants if p.chebi_id], ""),
        ]:
            if CHEBI_NITROQUINOLINE_N_OXIDE not in left_ids or CHEBI_HYDROXYAMINO_QUINOLINE_N_OXIDE not in right_ids:
                continue
            if CHEBI_H2O not in right_ids or left_ids.count(CHEBI_H_PLUS) != 2:
                continue
            if any(left_ids.count(reduced) == 2 and right_ids.count(oxidized) == 2 for reduced, oxidized in self.COFACTOR_PAIRS.items()):
                return ClassificationResult(
                    is_member=True,
                    explanation=f"nitroquinoline-N-oxide reductase [NAD(P)H] activity: reduction of nitroquinoline N-oxide with two nicotinamide equivalents{suffix}",
                )
        return ClassificationResult(is_member=False, explanation="Requires nitroquinoline N-oxide/hydroxyaminoquinoline N-oxide, water, and two matched nicotinamide equivalents")
