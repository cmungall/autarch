"""glycerophospholipase activity.

Catalysis of hydrolysis of glycerophospholipids to smaller lipid fragments and
phosphate-containing alcohol products.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.lipid_utils import (
    carbon_count,
    is_glycerophospholipid_like,
    is_lipid_like,
    oxygen_count,
    phosphate_count,
)
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class GlycerophospholipaseActivity(ReactionClass):
    """glycerophospholipase activity.

    Catalysis of hydrolysis of glycerophospholipids to smaller lipid fragments
    and phosphate-containing alcohol products.
    """

    GO_ID = "GO:0004620"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward.is_member:
            return forward
        reverse = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
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

        left_core = [
            participant
            for participant in left_participants
            if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in right_participants
            if participant.chebi_id != CHEBI_H_PLUS
        ]

        substrates = [participant for participant in left_core if is_glycerophospholipid_like(participant)]
        if len(substrates) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one glycerophospholipid substrate",
            )

        if len([participant for participant in right_core if carbon_count(participant) >= 2]) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Expected at least two substantive organic hydrolysis products",
            )

        has_lipid_fragment = any(
            (
                is_lipid_like(participant)
                or (carbon_count(participant) >= 5 and oxygen_count(participant) >= 3)
            )
            and participant not in substrates
            for participant in right_core
        )
        has_small_phosphorylated_fragment = any(
            phosphate_count(participant) >= 1
            and carbon_count(participant) <= 8
            and oxygen_count(participant) >= 3
            for participant in right_core
        )
        has_fatty_acid_release = any(
            participant.has_moiety(Moiety.CARBOXYL) and carbon_count(participant) >= 8
            for participant in right_core
        )
        if not (has_lipid_fragment and (has_small_phosphorylated_fragment or has_fatty_acid_release)):
            return ClassificationResult(
                is_member=False,
                explanation="Products do not match glycerophospholipid hydrolysis chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Glycerophospholipase activity: hydrolysis of a glycerophospholipid into smaller lipid fragments",
        )
