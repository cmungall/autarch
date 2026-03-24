"""limonene 1,2-monooxygenase [NAD(P)H] activity.

Catalysis of the epoxidation of limonene to limonene 1,2-epoxide with NADH or
NADPH as reductant.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_LIMONENE_R = "CHEBI:15382"
CHEBI_LIMONENE_S = "CHEBI:15383"
CHEBI_LIMONENE_EPOXIDE = "CHEBI:16431"


class Limonene12MonooxygenaseNADHOrNADPHActivity(ReactionClass):
    """limonene 1,2-monooxygenase [NAD(P)H] activity.

    Catalysis of the epoxidation of limonene to limonene 1,2-epoxide with NADH or
    NADPH as reductant.
    """

    GO_ID = "GO:0052601"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    COFACTOR_PAIRS = {
        CHEBI_NADPH: CHEBI_NADP_PLUS,
        CHEBI_NADH: CHEBI_NAD_PLUS,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward
        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(is_member=True, explanation=f"{reverse.explanation} (reverse reaction orientation)")
        return forward

    def _check_direction(self, left: list[Participant], right: list[Participant]) -> ClassificationResult:
        left_ids = {p.chebi_id for p in left if p.chebi_id}
        right_ids = {p.chebi_id for p in right if p.chebi_id}
        if not ({CHEBI_LIMONENE_R, CHEBI_LIMONENE_S} & left_ids):
            return ClassificationResult(is_member=False, explanation="Requires limonene as substrate")
        if CHEBI_LIMONENE_EPOXIDE not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires limonene 1,2-epoxide as product")
        if not {CHEBI_O2, CHEBI_H_PLUS} <= left_ids or CHEBI_H2O not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires dioxygen, hydron, and water in the monooxygenase stoichiometry")
        for reduced, oxidized in self.COFACTOR_PAIRS.items():
            if reduced in left_ids and oxidized in right_ids:
                return ClassificationResult(
                    is_member=True,
                    explanation="limonene 1,2-monooxygenase [NAD(P)H] activity: limonene epoxidation coupled to nicotinamide oxidation",
                )
        return ClassificationResult(is_member=False, explanation="Requires NADH/NADPH oxidation to the matching nicotinamide product")
