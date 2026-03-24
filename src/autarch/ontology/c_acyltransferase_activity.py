"""C-acyltransferase activity.

Catalysis of the transfer of an acyl group to a carbon acceptor.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_COA
from autarch.ontology.acyltransferase_acyl_groups_converted_into_alkyl_on_transfer import (
    AcyltransferaseAcylGroupsConvertedIntoAlkylOnTransfer,
)
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
)
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.ontology.succinyltransferase_activity import SuccinyltransferaseActivity


class CAcyltransferaseActivity(ReactionClass):
    """C-acyltransferase activity.

    Catalysis of the transfer of an acyl group to a carbon acceptor.
    """

    GO_ID = "GO:0016408"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    CHILD_CLASSES = (
        AcyltransferaseAcylGroupsConvertedIntoAlkylOnTransfer,
        SuccinyltransferaseActivity,
    )

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        aggregate = aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "C-acyltransferase activity",
        )
        if aggregate.is_member:
            return aggregate

        forward = self._check_direct_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward

        reverse = self._check_direct_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse reaction orientation)",
            )
        return aggregate

    def _check_direct_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        left_thioesters = [participant for participant in left_participants if participant.is_thioester()]
        right_thioesters = [participant for participant in right_participants if participant.is_thioester()]
        if len(left_thioesters) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="No supported multi-acyl C-acyltransferase donor pattern detected",
            )
        if not any(participant.chebi_id == CHEBI_COA for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires released coenzyme A",
            )
        if not right_thioesters:
            return ClassificationResult(
                is_member=False,
                explanation="Requires an acyl-thioester product",
            )

        return ClassificationResult(
            is_member=True,
            explanation="C-acyltransferase activity: condensation of acyl-thioester donors with CoA release",
        )
