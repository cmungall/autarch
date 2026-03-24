"""demethylase activity.

Catalysis of methyl-group removal, typically with release of formaldehyde or
formate.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_CO2, CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_O2
from autarch.ontology.lipid_utils import carbon_count
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_FORMALDEHYDE = "CHEBI:16842"
FORMATE_IDS = {"CHEBI:15740", "CHEBI:29947"}
CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
CHEBI_SUCCINATE = "CHEBI:15741"


class DemethylaseActivity(ReactionClass):
    """demethylase activity.

    Catalysis of methyl-group removal, typically with release of formaldehyde
    or formate.
    """

    GO_ID = "GO:0032451"
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
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        left_ids = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        if CHEBI_FORMALDEHYDE not in right_ids and not (FORMATE_IDS & right_ids):
            return ClassificationResult(is_member=False, explanation="Requires formaldehyde or formate release")

        oxidative_context = (
            CHEBI_O2 in left_ids
            or CHEBI_2_OXOGLUTARATE in left_ids
            or CHEBI_NADH in left_ids
            or CHEBI_NADPH in left_ids
        )
        if not oxidative_context:
            return ClassificationResult(is_member=False, explanation="Requires oxygenase or reductant context for oxidative demethylation")

        left_core = [
            participant
            for participant in left_participants
            if participant.chebi_id not in {CHEBI_O2, CHEBI_2_OXOGLUTARATE, CHEBI_NADH, CHEBI_NADPH, CHEBI_H_PLUS}
            and carbon_count(participant) >= 2
        ]
        right_core = [
            participant
            for participant in right_participants
            if participant.chebi_id not in {CHEBI_FORMALDEHYDE, CHEBI_SUCCINATE, CHEBI_CO2, CHEBI_H2O, CHEBI_H_PLUS}
            and participant.chebi_id not in FORMATE_IDS
            and carbon_count(participant) >= 2
        ]
        if not left_core or not right_core:
            return ClassificationResult(is_member=False, explanation="Missing conserved organic substrate/product branch")

        return ClassificationResult(is_member=True, explanation="Demethylase activity: oxidative one-carbon removal with formaldehyde/formate release")
