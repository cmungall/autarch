"""hydrolase acting on acid anhydrides in phosphorus-containing anhydrides.

Catalysis of the hydrolysis of any acid anhydride which contains phosphorus.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.hydrolase import Hydrolase


class HydrolaseActingOnAcidAnhydridesInPhosphorusContainingAnhydrides(Hydrolase):
    """hydrolase acting on acid anhydrides in phosphorus-containing anhydrides.

    Catalysis of the hydrolysis of any acid anhydride which contains
    phosphorus.
    """

    GO_ID = "GO:0016818"
    EC_NUMBER_PREFIX = "3.6.1.-"
    PHOSPHORIC_ANHYDRIDE_PATTERN = Chem.MolFromSmarts("[P][O][P]")
    ACYL_PHOSPHATE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])O[P]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolysis of phosphorus-containing acid anhydrides."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not phosphorus-anhydride hydrolysis",
            )

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
        """Evaluate one hydrolytic direction."""
        if not any(participant.chebi_id == CHEBI_H2O for participant in left):
            return ClassificationResult(
                is_member=False,
                explanation="No water substrate - not hydrolytic phosphorus-anhydride cleavage",
            )

        left_reactive = [participant for participant in left if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}]
        right_reactive = [participant for participant in right if participant.chebi_id != CHEBI_H_PLUS]
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive participants remain after removing water and hydrons",
            )

        left_anhydrides = self._anhydride_count(left_reactive)
        right_anhydrides = self._anhydride_count(right_reactive)
        if left_anhydrides == 0:
            return ClassificationResult(
                is_member=False,
                explanation="No phosphorus-containing acid anhydride detected on the substrate side",
            )
        if len(left_reactive) > 1 and not all(self._contains_anhydride(participant) for participant in left_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Additional non-anhydride substrates indicate ATP-coupled synthesis or transfer rather than direct anhydride hydrolysis",
            )
        if right_anhydrides >= left_anhydrides:
            return ClassificationResult(
                is_member=False,
                explanation="Anhydride count does not decrease across the hydrolysis reaction",
            )

        left_p = sum(self._phosphorus_count(participant) for participant in left_reactive)
        right_p = sum(self._phosphorus_count(participant) for participant in right_reactive)
        if left_p == 0 or left_p != right_p:
            return ClassificationResult(
                is_member=False,
                explanation="Phosphorus inventory is inconsistent with hydrolysis of a phosphorus-containing anhydride",
            )

        return ClassificationResult(
            is_member=True,
            explanation=(
                "Hydrolase acting on a phosphorus-containing acid anhydride: water-mediated cleavage "
                "reduces the number of P-containing anhydride bonds"
            ),
        )

    @classmethod
    def _anhydride_count(cls, participants: list[Participant]) -> int:
        return cls._pattern_count(participants, cls.PHOSPHORIC_ANHYDRIDE_PATTERN) + cls._pattern_count(
            participants,
            cls.ACYL_PHOSPHATE_PATTERN,
        )

    @staticmethod
    def _pattern_count(participants: list[Participant], pattern: Chem.Mol | None) -> int:
        if pattern is None:
            return 0
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            total += len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)
        return total

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    @classmethod
    def _contains_anhydride(cls, participant: Participant) -> bool:
        return cls._pattern_count([participant], cls.PHOSPHORIC_ANHYDRIDE_PATTERN) > 0 or cls._pattern_count(
            [participant],
            cls.ACYL_PHOSPHATE_PATTERN,
        ) > 0
