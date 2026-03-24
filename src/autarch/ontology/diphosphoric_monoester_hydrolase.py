"""diphosphoric monoester hydrolase.

Catalysis of the hydrolysis of a diphosphoester, releasing a diphosphate and a
free hydroxyl group.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_DIPHOSPHATE, CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.hydrolase import Hydrolase


class DiphosphoricMonoesterHydrolase(Hydrolase):
    """diphosphoric monoester hydrolase.

    Catalysis of the hydrolysis of a diphosphoester, releasing a diphosphate
    and a free hydroxyl group.
    """

    GO_ID = "GO:0016794"
    EC_NUMBER_PREFIX = "3.1.7.-"
    DIPHOSPHOESTER_PATTERN = Chem.MolFromSmarts(
        "[#6,#7,#8][OX2][PX4](=[OX1])([OX2,OX1-])[OX2][PX4](=[OX1])([OX2,OX1-])[OX2,OX1-]"
    )
    PRIMARY_DIPHOSPHATE_PATTERN = Chem.MolFromSmarts(
        "[CH2][OX2][PX4](=[OX1])([OX2,OX1-])[OX2][PX4](=[OX1])([OX2,OX1-])[OX2,OX1-]"
    )
    PRIMARY_ALCOHOL_PATTERN = Chem.MolFromSmarts("[CH2][OX2H]")
    PHOSPHATE_ESTER_PATTERN = Chem.MolFromSmarts("[#6][OX2][PX4]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolysis of an organic diphosphoester."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not diphosphoester hydrolysis",
            )

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
                explanation="Diphosphoric monoester hydrolases consume water",
            )

        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id != CHEBI_H_PLUS
        ]
        substrates = [
            participant
            for participant in left_reactive
            if self._carbon_count(participant) > 0 and self._phosphorus_count(participant) >= 2
        ]
        if len(substrates) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Requires one substantive organic diphosphoester substrate",
            )

        substrate = substrates[0]
        if not self._matches(substrate, self.DIPHOSPHOESTER_PATTERN):
            return ClassificationResult(
                is_member=False,
                explanation="Substrate lacks an organic diphosphoester linkage",
            )

        diphosphate_products = [
            participant
            for participant in right_reactive
            if participant.chebi_id == CHEBI_DIPHOSPHATE
            or (
                self._phosphorus_count(participant) == 2 and self._carbon_count(participant) == 0
            )
        ]
        if not diphosphate_products:
            return ClassificationResult(
                is_member=False,
                explanation="Hydrolysis does not release diphosphate",
            )

        remaining_products = [
            participant for participant in right_reactive if participant not in diphosphate_products
        ]
        organic_products = [
            participant for participant in remaining_products if self._carbon_count(participant) > 0
        ]
        if not organic_products:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive organic hydrolysis product remains after diphosphate release",
            )

        if self._carbon_total(organic_products) != self._carbon_count(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="Organic product carbon skeleton does not match the diphosphoester substrate",
            )

        if self._phosphorus_total(organic_products) != self._phosphorus_count(substrate) - 2:
            return ClassificationResult(
                is_member=False,
                explanation="Organic products do not lose exactly one diphosphate group",
            )

        remaining_p = self._phosphorus_total(organic_products)
        if remaining_p == 0:
            if not self._matches(substrate, self.PRIMARY_DIPHOSPHATE_PATTERN) or not any(
                self._matches(product, self.PRIMARY_ALCOHOL_PATTERN) for product in organic_products
            ):
                return ClassificationResult(
                    is_member=False,
                    explanation="Dephosphorylated products in this branch should arise from cleavage of a primary diphosphoester to a primary alcohol",
                )
        elif self._phosphate_ester_count(substrate) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Phosphate-retaining products require a substrate with multiple phosphate ester attachment sites",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Diphosphoric monoester hydrolase: water hydrolyzes an organic diphosphoester to release diphosphate",
        )

    @staticmethod
    def _matches(participant: Participant, pattern: Chem.Mol | None) -> bool:
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return False
        return mol.HasSubstructMatch(pattern)

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    @classmethod
    def _carbon_total(cls, participants: list[Participant]) -> int:
        return sum(cls._carbon_count(participant) * max(participant.count, 1) for participant in participants)

    @classmethod
    def _phosphorus_total(cls, participants: list[Participant]) -> int:
        return sum(
            cls._phosphorus_count(participant) * max(participant.count, 1) for participant in participants
        )

    def _phosphate_ester_count(self, participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return len(mol.GetSubstructMatches(self.PHOSPHATE_ESTER_PATTERN))
