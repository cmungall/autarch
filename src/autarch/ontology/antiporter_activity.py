"""antiporter activity."""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.transport_utils import TRANSPORT_SPECTATOR_IDS, transported_pairs


class AntiporterActivity(ReactionClass):
    """antiporter activity."""

    GO_ID = "GO:0015297"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    COUPLING_SPECTATOR_IDS = TRANSPORT_SPECTATOR_IDS - {CHEBI_H_PLUS}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        if not reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="No transported substrate between locations",
            )

        transported = [
            (left_participant, right_participant)
            for left_participant, right_participant in transported_pairs(reaction)
            if left_participant.chebi_id not in self.COUPLING_SPECTATOR_IDS
        ]
        if len(transported) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Antiport requires at least two transported substrates",
            )

        directions = {
            (left_participant.location, right_participant.location)
            for left_participant, right_participant in transported
        }
        if len(directions) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="No oppositely directed transport detected",
            )

        direction_text = ", ".join(
            f"{from_location}->{to_location}" for from_location, to_location in sorted(directions)
        )
        return ClassificationResult(
            is_member=True,
            explanation=f"Antiporter activity: oppositely directed transport ({direction_text})",
        )
