"""hydroxycinnamoyltransferase activity.

Catalysis of the transfer of a hydroxycinnamoyl group from a coenzyme A donor
to an acceptor substrate.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_COA, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

HYDROXYCINNAMOYL_COA_DONORS = {
    "CHEBI:57355",  # 4-coumaroyl-CoA
    "CHEBI:57393",  # sinapoyl-CoA
    "CHEBI:87136",  # caffeoyl-CoA
    "CHEBI:87305",  # feruloyl-CoA
}


class HydroxycinnamoyltransferaseActivity(ReactionClass):
    """hydroxycinnamoyltransferase activity.

    Catalysis of the transfer of a hydroxycinnamoyl group from a coenzyme A donor
    to an acceptor substrate.
    """

    GO_ID = "GO:0050734"
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
        donor = next(
            (participant for participant in left_participants if participant.chebi_id in HYDROXYCINNAMOYL_COA_DONORS),
            None,
        )
        if donor is None:
            return ClassificationResult(
                is_member=False,
                explanation="No supported hydroxycinnamoyl-CoA donor detected",
            )
        if not any(participant.chebi_id == CHEBI_COA for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires released coenzyme A",
            )

        left_core = [
            participant
            for participant in left_participants
            if participant is not donor and participant.chebi_id not in {CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in right_participants
            if participant.chebi_id not in {CHEBI_COA, CHEBI_H_PLUS}
        ]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one non-donor acceptor substrate and one hydroxycinnamoylated product",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Hydroxycinnamoyltransferase activity: hydroxycinnamoyl-CoA donor releases CoA during acyl transfer",
        )
