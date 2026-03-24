"""glyceollin synthase activity.

Catalysis of oxidative cyclization of prenylated pterocarpans to glyceollins.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE = "CHEBI:57618"
CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE = "CHEBI:58210"
SUBSTRATE_IDS = {"CHEBI:50036", "CHEBI:50118"}
PRODUCT_IDS = {"CHEBI:16470", "CHEBI:52127", "CHEBI:52086"}


class GlyceollinSynthaseActivity(ReactionClass):
    """glyceollin synthase activity.

    Catalysis of oxidative cyclization of prenylated pterocarpans to glyceollins.
    """

    GO_ID = "GO:0033769"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward
        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(is_member=True, explanation=f"{reverse.explanation} (reverse reaction orientation)")
        return forward

    def _check_direction(self, left: list[Participant], right: list[Participant]) -> ClassificationResult:
        left_ids = [p.chebi_id for p in left if p.chebi_id]
        right_ids = [p.chebi_id for p in right if p.chebi_id]
        left_set = set(left_ids)
        right_set = set(right_ids)
        if not (SUBSTRATE_IDS & left_set):
            return ClassificationResult(is_member=False, explanation="Requires a supported prenylated pterocarpan substrate")
        if not (PRODUCT_IDS & right_set):
            return ClassificationResult(is_member=False, explanation="Requires glyceollin product formation")
        if not {CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE, CHEBI_O2} <= left_set:
            return ClassificationResult(is_member=False, explanation="Requires reduced hemoprotein reductase and dioxygen")
        if not {CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE, CHEBI_H2O, CHEBI_H_PLUS} <= right_set:
            return ClassificationResult(is_member=False, explanation="Requires oxidized reductase, water, and hydron products")
        return ClassificationResult(
            is_member=True,
            explanation="glyceollin synthase activity: oxidative cyclization of prenylated pterocarpans",
        )
