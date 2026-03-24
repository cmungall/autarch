"""deacetylase activity.

Catalysis of the hydrolytic removal of an acetyl group with acetate release.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.lipid_utils import carbon_count
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

ACETATE_IDS = {"CHEBI:30089", "CHEBI:40480"}


class DeacetylaseActivity(ReactionClass):
    """deacetylase activity.

    Catalysis of the hydrolytic removal of an acetyl group with acetate
    release.
    """

    GO_ID = "GO:0019213"
    EC_BROAD_XREFS = ["3.5.1.-", "3.1.1.-", "2.3.1.-"]
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
        if not any(participant.chebi_id == CHEBI_H2O for participant in left_participants):
            return ClassificationResult(is_member=False, explanation="Requires water as a hydrolytic substrate")
        if not any(participant.chebi_id in ACETATE_IDS for participant in right_participants):
            return ClassificationResult(is_member=False, explanation="Requires acetate release")

        left_core = [participant for participant in left_participants if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}]
        right_core = [participant for participant in right_participants if participant.chebi_id not in ACETATE_IDS | {CHEBI_H_PLUS}]
        if not left_core or not right_core:
            return ClassificationResult(is_member=False, explanation="Requires a conserved non-acetate substrate/product branch")
        if not any(carbon_count(participant) >= 2 for participant in right_core):
            return ClassificationResult(is_member=False, explanation="No substantive deacetylated product detected")

        return ClassificationResult(is_member=True, explanation="Deacetylase activity: hydrolytic acetate release from an acetylated substrate")
