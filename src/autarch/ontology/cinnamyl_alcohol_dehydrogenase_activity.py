"""cinnamyl-alcohol dehydrogenase activity.

Catalysis of the oxidation of cinnamyl alcohol to cinnamaldehyde with NADP as
acceptor.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_NADPH, CHEBI_NADP_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_CINNAMYL_ALCOHOL = "CHEBI:33227"
CHEBI_CINNAMALDEHYDE = "CHEBI:16731"


class CinnamylAlcoholDehydrogenaseActivity(ReactionClass):
    """cinnamyl-alcohol dehydrogenase activity.

    Catalysis of the oxidation of cinnamyl alcohol to cinnamaldehyde with NADP as
    acceptor.
    """

    GO_ID = "GO:0045551"
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
        left_ids = {participant.chebi_id for participant in left if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right if participant.chebi_id}
        if not {CHEBI_CINNAMYL_ALCOHOL, CHEBI_NADP_PLUS} <= left_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires cinnamyl alcohol and NADP(+) on the substrate side",
            )
        if not {CHEBI_CINNAMALDEHYDE, CHEBI_NADPH} <= right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires cinnamaldehyde and NADPH as products",
            )
        if CHEBI_H_PLUS not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires hydron release",
            )
        return ClassificationResult(
            is_member=True,
            explanation="cinnamyl-alcohol dehydrogenase activity: oxidation of cinnamyl alcohol with NADP(+) reduction",
        )
