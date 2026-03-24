"""oxo-acid-lyase.

Catalysis of the cleavage of a C-C bond by other means than by hydrolysis or
oxidation, of a 3-hydroxy acid.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_CO2, CHEBI_H2O, CHEBI_H2O2, CHEBI_O2
from autarch.ontology.lyase import Lyase


class OxoAcidLyase(Lyase):
    """oxo-acid-lyase.

    Catalysis of the cleavage of a C-C bond by other means than by hydrolysis
    or oxidation, of a 3-hydroxy acid.
    """

    GO_ID = "GO:0016833"
    EC_NUMBER_PREFIX = "4.1.3.-"

    CARBOXYL_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[O;H1,X1-]")
    HYDROXYL_PATTERN = Chem.MolFromSmarts("[O;H1][C;!$(C=O)]")
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for non-hydrolytic cleavage of hydroxy-carboxylate substrates."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not oxo-acid lyase chemistry",
            )

        if any(participant.chebi_id in {CHEBI_H2O, CHEBI_O2} for participant in reaction.left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Hydrolytic or oxidative cosubstrates indicate a different reaction class",
            )

        if any(participant.chebi_id in {CHEBI_CO2, CHEBI_H2O2} for participant in reaction.right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Carbon dioxide or peroxide products indicate a different lyase/oxidoreductase branch",
            )

        left_reactive = [participant for participant in reaction.left_participants if self._has_carbon(participant)]
        right_reactive = [participant for participant in reaction.right_participants if self._has_carbon(participant)]
        if len(left_reactive) != 1 or len(right_reactive) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Oxo-acid lyases split one hydroxy-acid-like substrate into multiple carbon-containing products",
            )

        substrate = left_reactive[0]
        if not (substrate.is_carboxylic_acid() or substrate.is_thioester()):
            return ClassificationResult(
                is_member=False,
                explanation="No oxo-acid or thioester substrate detected",
            )

        if not self._has_hydroxyl(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="Oxo-acid lyases act on hydroxy-acid substrates",
            )

        if self._carboxyl_count([substrate]) < 2 and not substrate.is_thioester():
            return ClassificationResult(
                is_member=False,
                explanation="Precision-biased oxo-acid lyase rule requires a multi-acid substrate or a hydroxy-thioester",
            )

        if self._carboxyl_count(right_reactive) < 1 or self._carbonyl_count(right_reactive) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Products do not show the expected oxo-acid fragmentation signature",
            )

        if self._carbon_count(right_reactive) != self._carbon_count([substrate]):
            return ClassificationResult(
                is_member=False,
                explanation="Carbon skeleton is not conserved across the cleavage",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Oxo-acid lyase: non-hydrolytic C-C bond cleavage of a hydroxy-acid substrate",
        )

    @classmethod
    def _has_hydroxyl(cls, participant: Participant) -> bool:
        """Check for a non-carbonyl hydroxyl group."""
        mol = participant.get_mol()
        if mol is None or cls.HYDROXYL_PATTERN is None:
            return False
        return mol.HasSubstructMatch(cls.HYDROXYL_PATTERN)

    @classmethod
    def _carboxyl_count(cls, participants: list[Participant]) -> int:
        """Count carboxyl groups across participants."""
        return cls._pattern_count(participants, cls.CARBOXYL_PATTERN)

    @classmethod
    def _carbonyl_count(cls, participants: list[Participant]) -> int:
        """Count carbonyl groups across participants."""
        return cls._pattern_count(participants, cls.CARBONYL_PATTERN)

    @classmethod
    def _pattern_count(cls, participants: list[Participant], pattern: Chem.Mol | None) -> int:
        """Count SMARTS matches across participants."""
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
    def _has_carbon(participant: Participant) -> bool:
        """Check whether a participant has at least one carbon atom."""
        mol = participant.get_mol()
        if mol is None:
            return False
        return any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())

    @staticmethod
    def _carbon_count(participants: list[Participant]) -> int:
        """Count carbon atoms across participants."""
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            total += sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        return total
