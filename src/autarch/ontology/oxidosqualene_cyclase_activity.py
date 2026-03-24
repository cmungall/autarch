"""oxidosqualene cyclase activity.

Catalysis of the cyclization of (S)-2,3-epoxysqualene to form a triterpene.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_OXIDOSQUALENE = "CHEBI:15441"


class OxidosqualeneCyclaseActivity(ReactionClass):
    """oxidosqualene cyclase activity.

    Catalysis of the cyclization of (S)-2,3-epoxysqualene to form a triterpene.
    """

    GO_ID = "GO:0031559"
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
        left: list[Participant],
        right: list[Participant],
    ) -> ClassificationResult:
        if not any(participant.chebi_id == CHEBI_OXIDOSQUALENE for participant in left):
            return ClassificationResult(
                is_member=False,
                explanation="Requires (S)-2,3-epoxysqualene as the substrate scaffold",
            )

        right_products = [
            participant
            for participant in right
            if participant.chebi_id != CHEBI_H2O
        ]
        cyclized_products = [
            participant
            for participant in right_products
            if participant.chebi_id != CHEBI_OXIDOSQUALENE and self._is_cyclized_triterpene(participant)
        ]
        if not cyclized_products:
            return ClassificationResult(
                is_member=False,
                explanation="Requires a cyclized triterpene-like product",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Oxidosqualene cyclase activity: cyclization of (S)-2,3-epoxysqualene to a triterpene scaffold",
        )

    @staticmethod
    def _is_cyclized_triterpene(participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        ring_count = mol.GetRingInfo().NumRings()
        return carbon_count >= 25 and ring_count >= 4
