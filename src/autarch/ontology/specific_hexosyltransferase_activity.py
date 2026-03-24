"""Helpers for donor-specific hexosyltransferase GO classes."""

from __future__ import annotations

from typing import ClassVar

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.ontology.hexosyltransferase import Hexosyltransferase
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE


class SpecificHexosyltransferaseActivity(Hexosyltransferase):
    """Base class for donor-restricted hexosyltransferase activities."""

    EXCLUDE_FROM_DISCOVERY: ClassVar[bool] = True
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    DONOR_IDS: ClassVar[frozenset[str]]
    PRODUCT_IDS: ClassVar[frozenset[str]]
    CONCEPT_PHRASE: ClassVar[str]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not glycosyl transfer",
            )

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
                explanation=f"No {self.CONCEPT_PHRASE} donor detected",
            )

        if not any(participant.chebi_id in self.PRODUCT_IDS for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation=f"No released carrier for {self.CONCEPT_PHRASE} detected",
            )

        result = self._check_direction(left_participants, right_participants)
        if not result.is_member:
            return result

        return ClassificationResult(
            is_member=True,
            explanation=f"{self.CONCEPT_PHRASE}: donor-restricted hexosyl transfer",
        )
