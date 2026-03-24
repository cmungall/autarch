"""carboxylic ester hydrolase.

Catalysis of the hydrolysis of a carboxylic ester bond.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
)
from autarch.ontology.reaction import ReactionClass


class CarboxylicEsterHydrolase(ReactionClass):
    """carboxylic ester hydrolase.

    Catalysis of the hydrolysis of a carboxylic ester bond.
    """

    GO_ID = "GO:0052689"  # carboxylic ester hydrolase activity
    EC_NUMBER_PREFIX = "3.1.1.-"
    SPECTATOR_LEFT_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    SPECTATOR_RIGHT_CHEBIS = {CHEBI_H_PLUS}
    CARBOXYL_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[O;H1,X1-]")
    HYDROXYL_PATTERN = Chem.MolFromSmarts("[O;H1][#6;!$(C=O)]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolysis of a carboxylic ester or lactone bond."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not carboxylic ester hydrolysis",
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
        """Evaluate one reaction orientation for ester hydrolysis."""
        has_water = any(participant.chebi_id == CHEBI_H2O for participant in left_participants)
        if not has_water:
            return ClassificationResult(
                is_member=False,
                explanation="No water reactant found",
            )

        if any(participant.chebi_id == CHEBI_O2 for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Uses O2 - oxidative chemistry, not ester hydrolysis",
            )

        if any(participant.chebi_id == CHEBI_H2O2 for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Produces H2O2 - oxidative chemistry, not ester hydrolysis",
            )

        if any(participant.chebi_id == CHEBI_ATP for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Uses ATP - not simple ester hydrolysis",
            )

        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.SPECTATOR_LEFT_CHEBIS
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.SPECTATOR_RIGHT_CHEBIS
        ]

        if len(left_reactive) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Carboxylic ester hydrolases act on a single substantive ester substrate",
            )

        if any(participant.chebi_id in {CHEBI_NH3, CHEBI_NH4} for participant in right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Produces ammonia/ammonium - likely amidase chemistry",
            )

        if any(participant.is_thioester() for participant in left_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Thioester substrate - carbon-sulfur hydrolase chemistry, not carboxylic ester hydrolysis",
            )

        left_ester_bonds = sum(
            self._carboxylic_ester_bond_count(participant) * max(participant.count, 1)
            for participant in left_reactive
        )
        if left_ester_bonds == 0:
            return ClassificationResult(
                is_member=False,
                explanation="No carboxylic ester substrate detected",
            )

        right_ester_bonds = sum(
            self._carboxylic_ester_bond_count(participant) * max(participant.count, 1)
            for participant in right_reactive
        )
        if right_ester_bonds >= left_ester_bonds:
            return ClassificationResult(
                is_member=False,
                explanation="No carboxylic ester bond is cleaved",
            )

        left_carboxyl_count = sum(
            self._substructure_count(participant, self.CARBOXYL_PATTERN) * max(participant.count, 1)
            for participant in left_reactive
        )
        right_carboxyl_count = sum(
            self._substructure_count(participant, self.CARBOXYL_PATTERN) * max(participant.count, 1)
            for participant in right_reactive
        )
        if right_carboxyl_count > left_carboxyl_count:
            return ClassificationResult(
                is_member=True,
                explanation="Carboxylic ester hydrolase: ester bond cleavage decreases ester count and increases carboxylate products",
            )

        left_hydroxyl_count = sum(
            self._substructure_count(participant, self.HYDROXYL_PATTERN) * max(participant.count, 1)
            for participant in left_reactive
        )
        right_hydroxyl_count = sum(
            self._substructure_count(participant, self.HYDROXYL_PATTERN) * max(participant.count, 1)
            for participant in right_reactive
        )
        produces_co2 = any(
            participant.chebi_id == CHEBI_CO2 or participant.smiles == "O=C=O"
            for participant in right_reactive
        )
        if produces_co2 and right_hydroxyl_count > left_hydroxyl_count:
            return ClassificationResult(
                is_member=True,
                explanation="Carboxylic ester hydrolase: ester hydrolysis followed by decarboxylation",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No carboxylate increase or decarboxylative ester-hydrolysis signature detected",
        )

    @staticmethod
    def _substructure_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        if pattern is None:
            return 0
        mol = participant.get_mol()
        if mol is None:
            return 0
        return len(mol.GetSubstructMatches(pattern))

    @staticmethod
    def _carboxylic_ester_bond_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0

        ester_bonds = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue

            double_bonded_oxygens = []
            single_bonded_oxygens = []
            non_oxygen_neighbors = []
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8:
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        double_bonded_oxygens.append(neighbor)
                    else:
                        single_bonded_oxygens.append(neighbor)
                else:
                    non_oxygen_neighbors.append(neighbor)

            if not double_bonded_oxygens or not single_bonded_oxygens:
                continue
            if non_oxygen_neighbors and any(neighbor.GetAtomicNum() != 6 for neighbor in non_oxygen_neighbors):
                continue

            for oxygen in single_bonded_oxygens:
                has_carbon_neighbor = any(
                    neighbor.GetIdx() != atom.GetIdx() and neighbor.GetAtomicNum() == 6
                    for neighbor in oxygen.GetNeighbors()
                )
                if has_carbon_neighbor:
                    ester_bonds += 1

        return ester_bonds
