"""S-acyltransferase activity.

Catalysis of the transfer of an acyl group to a sulfur atom on the acceptor molecule.
"""

from __future__ import annotations

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_COA, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_PANTETHEINE_4P = "CHEBI:61723"


class SAcyltransferaseActivity(ReactionClass):
    """S-acyltransferase activity.

    Catalysis of the transfer of an acyl group to a sulfur atom on the acceptor molecule.
    """

    GO_ID = "GO:0016417"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    ACYL_SULFUR_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[S;X1,X2,-]")

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
        left_core = [participant for participant in left if participant.chebi_id != CHEBI_H_PLUS]
        right_core = [participant for participant in right if participant.chebi_id != CHEBI_H_PLUS]

        donors = [participant for participant in left_core if participant.is_thioester()]
        sulfur_acceptors = [
            participant
            for participant in left_core
            if not participant.is_thioester() and self._has_sulfur(participant)
        ]
        released_carriers = [
            participant
            for participant in right_core
            if participant.chebi_id in {CHEBI_COA, CHEBI_PANTETHEINE_4P}
        ]
        sulfur_products = [
            participant
            for participant in right_core
            if participant.chebi_id not in {CHEBI_COA, CHEBI_PANTETHEINE_4P}
            if self._is_acyl_sulfur_product(participant)
        ]

        if not donors:
            return ClassificationResult(
                is_member=False,
                explanation="Requires an activated thioester acyl donor",
            )
        if not sulfur_acceptors:
            return ClassificationResult(
                is_member=False,
                explanation="Requires a sulfur-containing acceptor on the donor side",
            )
        if not released_carriers:
            return ClassificationResult(
                is_member=False,
                explanation="Requires CoA or phosphopantetheine carrier release",
            )
        if not sulfur_products:
            return ClassificationResult(
                is_member=False,
                explanation="Requires formation of a new sulfur-acyl product",
            )

        return ClassificationResult(
            is_member=True,
            explanation="S-acyltransferase activity: transfer of an activated acyl group to a sulfur acceptor with carrier release",
        )

    @staticmethod
    def _has_sulfur(participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        return any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())

    @classmethod
    def _is_acyl_sulfur_product(cls, participant: Participant) -> bool:
        if participant.is_thioester():
            return True
        mol = participant.get_mol()
        if mol is None or cls.ACYL_SULFUR_PATTERN is None:
            return False
        return mol.HasSubstructMatch(cls.ACYL_SULFUR_PATTERN)
