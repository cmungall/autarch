"""sterol 14-demethylase activity.

Catalysis of sterol 14-demethylation using molecular oxygen and reduced
NADPH-hemoprotein reductase equivalents.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE = "CHEBI:57618"
CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE = "CHEBI:58210"
CHEBI_FORMATE = "CHEBI:15740"
SUBSTRATE_IDS = {"CHEBI:17791", "CHEBI:16521"}
PRODUCT_IDS = {"CHEBI:30109", "CHEBI:17813"}


class Sterol14DemethylaseActivity(ReactionClass):
    """sterol 14-demethylase activity.

    Catalysis of sterol 14-demethylation using molecular oxygen and reduced
    NADPH-hemoprotein reductase equivalents.
    """

    GO_ID = "GO:0008398"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward

        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse reaction orientation)",
            )

        return forward

    def _check_direction(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in left if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in right if participant.chebi_id]
        left_set = set(left_ids)
        right_set = set(right_ids)
        if not (SUBSTRATE_IDS & left_set):
            return ClassificationResult(is_member=False, explanation="Requires a supported sterol substrate")
        if not (PRODUCT_IDS & right_set):
            return ClassificationResult(is_member=False, explanation="Requires a supported demethylated sterol product")
        if left_ids.count(CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE) < 3 or right_ids.count(CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE) < 3:
            return ClassificationResult(is_member=False, explanation="Requires the reduced/oxidized hemoprotein reductase cofactor set")
        if left_ids.count(CHEBI_O2) < 3 or right_ids.count(CHEBI_H2O) < 3:
            return ClassificationResult(is_member=False, explanation="Requires multi-oxygen demethylation stoichiometry")
        if CHEBI_FORMATE not in right_set:
            return ClassificationResult(is_member=False, explanation="Requires formate release")
        if right_ids.count(CHEBI_H_PLUS) < 1:
            return ClassificationResult(is_member=False, explanation="Requires hydron release")
        return ClassificationResult(
            is_member=True,
            explanation="sterol 14-demethylase activity: oxygen-dependent sterol demethylation with formate release",
        )
