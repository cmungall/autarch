"""sesquiterpene synthase activity.

Catalysis of the reaction: trans,trans-farnesyl diphosphate = a sesquiterpene +
diphosphate.
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_DIPHOSPHATE, CHEBI_H2O
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_FARNESYL_DIPHOSPHATE = "CHEBI:175763"
CHEBI_ALT_FARNESYL_DIPHOSPHATE = "CHEBI:57223"


class SesquiterpeneSynthaseActivity(ReactionClass):
    """sesquiterpene synthase activity.

    Catalysis of the reaction: trans,trans-farnesyl diphosphate = a
    sesquiterpene + diphosphate.
    """

    GO_ID = "GO:0010334"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    DONOR_IDS = {CHEBI_FARNESYL_DIPHOSPHATE, CHEBI_ALT_FARNESYL_DIPHOSPHATE}

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
        if any(participant.chebi_id == CHEBI_H2O for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Uses water as substrate - not sesquiterpene synthase chemistry",
            )

        left_ids = [participant.chebi_id for participant in left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in right_participants if participant.chebi_id]
        if not any(chebi_id in self.DONOR_IDS for chebi_id in left_ids):
            return ClassificationResult(
                is_member=False,
                explanation="Requires a farnesyl diphosphate donor",
            )
        if CHEBI_DIPHOSPHATE not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires diphosphate release",
            )

        left_core = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.DONOR_IDS
        ]
        right_core = [
            participant
            for participant in right_participants
            if participant.chebi_id != CHEBI_DIPHOSPHATE
        ]
        if any(participant.chebi_id not in {None, CHEBI_H2O} for participant in left_core):
            return ClassificationResult(
                is_member=False,
                explanation="Expected farnesyl diphosphate to be the only substantive organic substrate",
            )
        if len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one substantive sesquiterpene product",
            )
        if right_core[0].is_thioester() or right_core[0].is_phosphorylated():
            return ClassificationResult(
                is_member=False,
                explanation="Product retains a high-energy transfer group rather than a sesquiterpene scaffold",
            )
        return ClassificationResult(
            is_member=True,
            explanation="Sesquiterpene synthase activity: cyclization of farnesyl diphosphate with diphosphate release",
        )
