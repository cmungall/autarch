"""thiolester hydrolase.

Catalysis of the reaction: RCO-SR' + H2O = RCOOH + HSR'. This reaction is the
hydrolysis of a thiolester bond, an ester formed from a carboxylic acid and a
thiol (i.e., RCO-SR'), such as that found in acetyl-coenzyme A.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_O2,
)
from autarch.ontology.hydrolase import Hydrolase


class ThiolesterHydrolase(Hydrolase):
    """thiolester hydrolase.

    Catalysis of the reaction: RCO-SR' + H2O = RCOOH + HSR'. This reaction is
    the hydrolysis of a thiolester bond, an ester formed from a carboxylic acid
    and a thiol (i.e., RCO-SR'), such as that found in acetyl-coenzyme A.
    """

    GO_ID = "GO:0016790"
    EC_NUMBER_PREFIX = "3.1.2.-"
    LEFT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H_PLUS}
    THIOESTER_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    CARBOXYL_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[O;H1,X1-]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolysis of a thiolester bond."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not thiolester hydrolysis",
            )

        forward_result = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward_result.is_member:
            return forward_result

        reverse_result = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
        if reverse_result.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse_result.explanation} (reverse reaction orientation)",
            )

        return forward_result

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one reaction orientation."""
        parent_result = super().check_membership_impl(
            Reaction(left_participants=left_participants, right_participants=right_participants)
        )
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase in this reaction orientation: {parent_result.explanation}",
            )

        if any(participant.chebi_id in {CHEBI_ATP, CHEBI_O2, CHEBI_H2O2} for participant in left_participants + right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="ATP or oxygen chemistry indicates a different reaction class",
            )

        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
        ]

        if len(left_reactive) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Thiolester hydrolases act on a single substantive thiolester substrate",
            )

        left_thiolester_count = self._substructure_count(left_reactive, self.THIOESTER_PATTERN)
        if left_thiolester_count == 0:
            return ClassificationResult(
                is_member=False,
                explanation="No thiolester substrate detected",
            )

        right_thiolester_count = self._substructure_count(right_reactive, self.THIOESTER_PATTERN)
        if right_thiolester_count >= left_thiolester_count:
            return ClassificationResult(
                is_member=False,
                explanation="No thiolester bond is cleaved",
            )

        left_carboxyl_count = self._substructure_count(left_reactive, self.CARBOXYL_PATTERN)
        right_carboxyl_count = self._substructure_count(right_reactive, self.CARBOXYL_PATTERN)
        if right_carboxyl_count <= left_carboxyl_count:
            return ClassificationResult(
                is_member=False,
                explanation="Hydrolysis does not generate a carboxylate product from the thiolester substrate",
            )

        has_sulfur_product = any(
            self._atom_count(participant, 16) > 0 and self._substructure_count([participant], self.THIOESTER_PATTERN) == 0
            for participant in right_reactive
        )
        if not has_sulfur_product:
            return ClassificationResult(
                is_member=False,
                explanation="No released sulfur-containing product consistent with thiolester cleavage",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Thiolester hydrolase: hydrolytic cleavage of a thiolester bond to carboxylate and sulfur-containing products",
        )

    @staticmethod
    def _substructure_count(participants: list[Participant], pattern: Chem.Mol | None) -> int:
        """Count pattern matches across participants."""
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
    def _atom_count(participant: Participant, atomic_num: int) -> int:
        """Count atoms of a given element."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_num)
