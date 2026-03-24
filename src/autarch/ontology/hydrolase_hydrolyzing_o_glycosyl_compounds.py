"""hydrolase activity, hydrolyzing O-glycosyl compounds.

Catalysis of the hydrolysis of any O-glycosyl bond.
"""

from typing import ClassVar, Optional

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.hydrolase import Hydrolase


class HydrolaseHydrolyzingOGlycosylCompounds(Hydrolase):
    """hydrolase activity, hydrolyzing O-glycosyl compounds.

    Catalysis of the hydrolysis of any O-glycosyl bond.
    """

    GO_ID = "GO:0004553"
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "3.2.1.-"
    SUGAR_CHEBIS = {
        "CHEBI:4167",
        "CHEBI:15903",
        "CHEBI:16551",
        "CHEBI:27667",
        "CHEBI:17218",
        "CHEBI:18079",
        "CHEBI:47013",
        "CHEBI:28009",
        "CHEBI:35418",
    }
    HYDROXYL_PATTERN = Chem.MolFromSmarts("[OX2H]")
    EXOCYCLIC_O_GLYCOSIDE_PATTERN = Chem.MolFromSmarts("[C;R][O;X2;!R][C]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolysis of an O-glycosidic bond."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase: {parent_result.explanation}",
            )

        has_water = any(
            participant.chebi_id == CHEBI_H2O for participant in reaction.left_participants
        )
        if not has_water:
            return ClassificationResult(
                is_member=False,
                explanation="No water reactant - not O-glycosyl hydrolysis",
            )

        left_core = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id != CHEBI_H_PLUS
        ]

        has_o_glycoside_substrate = any(self._is_o_glycosidic_substrate(participant) for participant in left_core)
        if not has_o_glycoside_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="No O-glycosidic substrate detected",
            )

        has_sugar_like_product = any(self._is_sugar_like(participant) for participant in right_core)
        if not has_sugar_like_product:
            return ClassificationResult(
                is_member=False,
                explanation="No sugar-like hydrolysis product detected",
            )

        organic_products = [participant for participant in right_core if self._is_organic_product(participant)]
        if len(organic_products) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="O-glycosyl hydrolysis should yield two organic product fragments rather than simple deprotection products",
            )

        return ClassificationResult(
            is_member=True,
            explanation="O-glycosyl hydrolase: water-dependent cleavage of an O-glycosidic bond",
        )

    @classmethod
    def _is_sugar_like(cls, participant: Participant) -> bool:
        if participant.chebi_id in cls.SUGAR_CHEBIS:
            return True
        mol = participant.get_mol()
        if mol is None or cls.HYDROXYL_PATTERN is None:
            return False
        oxygen_count = sum(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
        hydroxyl_count = len(mol.GetSubstructMatches(cls.HYDROXYL_PATTERN))
        has_ring_oxygen = any(atom.GetAtomicNum() == 8 and atom.IsInRing() for atom in mol.GetAtoms())
        return oxygen_count >= 4 and hydroxyl_count >= 2 and has_ring_oxygen

    @classmethod
    def _is_o_glycosidic_substrate(cls, participant: Participant) -> bool:
        """Check for an O-glycosidic linkage without using labels."""
        if not participant.smiles:
            return False
        if participant.has_moiety(Moiety.O_GLYCOSIDE):
            return True
        mol = participant.get_mol()
        if mol is None or cls.EXOCYCLIC_O_GLYCOSIDE_PATTERN is None:
            return False
        if not cls._is_sugar_like(participant):
            return False
        return mol.HasSubstructMatch(cls.EXOCYCLIC_O_GLYCOSIDE_PATTERN)

    @staticmethod
    def _is_organic_product(participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
        return carbon_count >= 2
