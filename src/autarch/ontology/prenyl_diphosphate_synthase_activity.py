"""prenyl diphosphate synthase activity.

Catalysis of the chain-elongation reaction between prenyl diphosphates,
releasing diphosphate and forming a longer prenyl diphosphate.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_DIPHOSPHATE, CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.lipid_utils import carbon_count, is_prenyl_diphosphate_like
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class PrenylDiphosphateSynthaseActivity(ReactionClass):
    """prenyl diphosphate synthase activity.

    Catalysis of the chain-elongation reaction between prenyl diphosphates,
    releasing diphosphate and forming a longer prenyl diphosphate.
    """

    GO_ID = "GO:0120531"
    EC_BROAD_XREFS = ["2.5.1.-"]
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    LEFT_SPECTATOR_IDS = {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_IDS = {CHEBI_DIPHOSPHATE, CHEBI_H_PLUS}

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
        if any(participant.chebi_id == CHEBI_H2O for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Water-reactive chemistry is not prenyl diphosphate synthesis",
            )

        if not any(participant.chebi_id == CHEBI_DIPHOSPHATE for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires diphosphate release",
            )

        left_core = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_IDS
        ]
        right_core = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_IDS
        ]

        left_prenyl = [participant for participant in left_core if is_prenyl_diphosphate_like(participant)]
        right_prenyl = [participant for participant in right_core if is_prenyl_diphosphate_like(participant)]
        if len(left_prenyl) != 2 or len(right_prenyl) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected two prenyl diphosphate substrates and one longer prenyl diphosphate product",
            )
        if len(left_core) != 2 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Unexpected non-prenyl participants in chain-elongation chemistry",
            )

        substrate_carbons = sorted(carbon_count(participant) for participant in left_prenyl)
        product_carbons = carbon_count(right_prenyl[0])
        if 5 not in substrate_carbons:
            return ClassificationResult(
                is_member=False,
                explanation="Prenyl diphosphate synthases require an isopentenyl-sized C5 donor",
            )
        if product_carbons != sum(substrate_carbons):
            return ClassificationResult(
                is_member=False,
                explanation="Product carbon count does not match prenyl chain elongation",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Prenyl diphosphate synthase activity: C5 prenyl extension with diphosphate release",
        )
