"""UDP-glycosyltransferase activity."""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_UDP
from autarch.ontology.glycosyltransferase import Glycosyltransferase
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE


class UDPGlycosyltransferaseActivity(Glycosyltransferase):
    """UDP-glycosyltransferase activity."""

    GO_ID = "GO:0008194"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    CONCEPT_PHRASE = "UDP-glycosyltransferase activity"
    DONOR_IDS = frozenset(
        {
            "CHEBI:18066",
            "CHEBI:18307",
            "CHEBI:58885",
            "CHEBI:66914",
            "CHEBI:57705",
            "CHEBI:58052",
            "CHEBI:57632",
            "CHEBI:67138",
            "CHEBI:68623",
            "CHEBI:83836",
            "CHEBI:70731",
        }
    )
    PRODUCT_IDS = frozenset({CHEBI_UDP})

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_specific_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward.is_member:
            return forward

        reverse = self._check_specific_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse reaction orientation)",
            )

        return forward

    def _check_specific_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        if not any(participant.chebi_id in self.DONOR_IDS for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No UDP-sugar donor detected",
            )
        if not any(participant.chebi_id in self.PRODUCT_IDS for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No released UDP detected",
            )

        result = super().check_membership_impl(
            Reaction(
                left_participants=left_participants,
                right_participants=right_participants,
            )
        )
        if not result.is_member:
            return result

        return ClassificationResult(
            is_member=True,
            explanation="UDP-glycosyltransferase activity: UDP-sugar donor transfers a glycosyl group",
        )
