"""steroid hydroxylase activity.

Catalysis of a reaction in which one atom of molecular oxygen is incorporated
into a steroid substrate.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.biochemical_context_utils import is_steroid_like
from autarch.ontology.lipid_utils import carbon_count, oxygen_count
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class SteroidHydroxylaseActivity(ReactionClass):
    """steroid hydroxylase activity.

    Catalysis of a reaction in which one atom of molecular oxygen is
    incorporated into a steroid substrate.
    """

    GO_ID = "GO:0008395"
    EC_BROAD_XREFS = ["1.-.-.-"]
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward
        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(is_member=True, explanation=f"{reverse.explanation} (reverse reaction orientation)")
        return forward

    @staticmethod
    def _check_direction(left_participants: list[Participant], right_participants: list[Participant]) -> ClassificationResult:
        left_ids = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        if CHEBI_O2 not in left_ids or CHEBI_H2O not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires O2 consumption with H2O production")

        left_core = [participant for participant in left_participants if participant.chebi_id not in {CHEBI_O2, CHEBI_H_PLUS}]
        right_core = [participant for participant in right_participants if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}]
        left_steroids = [participant for participant in left_core if is_steroid_like(participant)]
        right_steroids = [participant for participant in right_core if is_steroid_like(participant)]
        if not left_steroids or not right_steroids:
            return ClassificationResult(is_member=False, explanation="Requires steroid-like substrate and product")

        for substrate in left_steroids:
            for product in right_steroids:
                if carbon_count(substrate) == carbon_count(product) and oxygen_count(product) >= oxygen_count(substrate) + 1:
                    return ClassificationResult(
                        is_member=True,
                        explanation="Steroid hydroxylase activity: oxygen-dependent hydroxylation of a steroid-like scaffold",
                    )

        return ClassificationResult(is_member=False, explanation="No steroid-like hydroxylation pattern detected")
