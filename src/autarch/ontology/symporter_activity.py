"""symporter activity."""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.transport_utils import TRANSPORT_SPECTATOR_IDS, transported_pairs


class SymporterActivity(ReactionClass):
    """symporter activity."""

    GO_ID = "GO:0015293"
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
                explanation="Symport requires at least two transported substrates",
            )

        directions = {
            (left_participant.location, right_participant.location)
            for left_participant, right_participant in transported
        }
        if len(directions) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Transport directions are not aligned",
            )

        direction = next(iter(directions))
        return ClassificationResult(
            is_member=True,
            explanation=f"Symporter activity: coupled transport in one direction {direction[0]}->{direction[1]}",
        )
