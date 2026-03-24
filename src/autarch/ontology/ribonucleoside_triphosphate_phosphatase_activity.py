"""ribonucleoside triphosphate phosphatase activity.

Catalysis of the hydrolysis of a ribonucleoside triphosphate to the
corresponding diphosphate and phosphate.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_PHOSPHATE
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

PHOSPHATE_IDS = {CHEBI_PHOSPHATE, "CHEBI:16838", "CHEBI:18367"}
TRIPHOSPHATE_TO_DIPHOSPHATE = {
    "CHEBI:30616": "CHEBI:456216",  # ATP -> ADP
    "CHEBI:37565": "CHEBI:58189",  # GTP -> GDP
    "CHEBI:61557": "CHEBI:57930",  # CTP -> CDP
    "CHEBI:61402": "CHEBI:58280",  # UTP -> UDP
    "CHEBI:37563": "CHEBI:58069",  # ITP -> IDP
}


class RibonucleosideTriphosphatePhosphataseActivity(ReactionClass):
    """ribonucleoside triphosphate phosphatase activity.

    Catalysis of the hydrolysis of a ribonucleoside triphosphate to the
    corresponding diphosphate and phosphate.
    """

    GO_ID = "GO:0017111"
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
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        if not any(participant.chebi_id == CHEBI_H2O for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires water as a hydrolytic substrate",
            )
        if not any(participant.chebi_id in PHOSPHATE_IDS for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires phosphate release",
            )

        left_ids = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        for triphosphate, diphosphate in TRIPHOSPHATE_TO_DIPHOSPHATE.items():
            if triphosphate in left_ids and diphosphate in right_ids:
                return ClassificationResult(
                    is_member=True,
                    explanation="Ribonucleoside triphosphate phosphatase activity: hydrolytic conversion of an NTP to the matching NDP",
                )

        return ClassificationResult(
            is_member=False,
            explanation="No supported ribonucleoside triphosphate to diphosphate conversion detected",
        )
