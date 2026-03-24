"""zeaxanthin epoxidase activity.

Catalysis of the epoxidation of zeaxanthin or antheraxanthin using molecular
oxygen and reduced ferredoxin.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_REDUCED_FERREDOXIN = "CHEBI:33738"
CHEBI_OXIDIZED_FERREDOXIN = "CHEBI:33737"
CHEBI_ZEAXANTHIN = "CHEBI:27547"
CHEBI_ANTHERAXANTHIN = "CHEBI:27867"
CHEBI_VIOLAXANTHIN = "CHEBI:35288"


class ZeaxanthinEpoxidaseActivity(ReactionClass):
    """zeaxanthin epoxidase activity.

    Catalysis of the epoxidation of zeaxanthin or antheraxanthin using molecular
    oxygen and reduced ferredoxin.
    """

    GO_ID = "GO:0052662"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    PRODUCTS_BY_SUBSTRATE = {
        CHEBI_ZEAXANTHIN: {CHEBI_ANTHERAXANTHIN, CHEBI_VIOLAXANTHIN},
        CHEBI_ANTHERAXANTHIN: {CHEBI_VIOLAXANTHIN},
    }

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
        left_ids = [participant.chebi_id for participant in left if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in right if participant.chebi_id]
        left_set = set(left_ids)
        right_set = set(right_ids)

        substrate = next((chebi_id for chebi_id in self.PRODUCTS_BY_SUBSTRATE if chebi_id in left_set), None)
        if substrate is None:
            return ClassificationResult(
                is_member=False,
                explanation="Requires zeaxanthin or antheraxanthin as the epoxidation substrate",
            )
        if not (self.PRODUCTS_BY_SUBSTRATE[substrate] & right_set):
            return ClassificationResult(
                is_member=False,
                explanation="Requires antheraxanthin or violaxanthin as the epoxidized product",
            )
        if left_ids.count(CHEBI_REDUCED_FERREDOXIN) < 2 or right_ids.count(CHEBI_OXIDIZED_FERREDOXIN) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Requires reduced and oxidized ferredoxin cofactors",
            )
        if CHEBI_O2 not in left_set or CHEBI_H2O not in right_set:
            return ClassificationResult(
                is_member=False,
                explanation="Requires dioxygen consumption and water formation",
            )
        if left_ids.count(CHEBI_H_PLUS) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Requires hydrons for the epoxidation chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Zeaxanthin epoxidase activity: ferredoxin-coupled epoxidation of zeaxanthin/antheraxanthin",
        )
