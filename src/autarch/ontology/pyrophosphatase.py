"""pyrophosphatase.

Catalysis of the hydrolysis of a phosphorus-oxygen-phosphorus bond in a
phosphoric anhydride substrate.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.hydrolase import Hydrolase


class Pyrophosphatase(Hydrolase):
    """pyrophosphatase.

    Catalysis of the hydrolysis of a phosphorus-oxygen-phosphorus bond in a
    phosphoric anhydride substrate.
    """

    GO_ID = "GO:0016462"  # pyrophosphatase activity
    EC_NUMBER_PREFIX = "3.6.1.-"
    SPECTATOR_LEFT_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    SPECTATOR_RIGHT_CHEBIS = {CHEBI_H_PLUS}
    PHOSPHORIC_ANHYDRIDE_PATTERN = Chem.MolFromSmarts("[P][O][P]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolysis of a single phosphoric anhydride substrate."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not pyrophosphatase",
            )

        water_count = sum(
            participant.count
            for participant in reaction.left_participants
            if participant.chebi_id == CHEBI_H2O
        )
        if water_count == 0:
            return ClassificationResult(
                is_member=False,
                explanation="No water substrate - pyrophosphatases are hydrolases",
            )
        if water_count != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Pyrophosphatase chemistry hydrolyzes a single phosphoric anhydride bond with one water equivalent",
            )

        left_reactive = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in self.SPECTATOR_LEFT_CHEBIS
        ]
        right_reactive = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in self.SPECTATOR_RIGHT_CHEBIS
        ]

        if len(left_reactive) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Pyrophosphatase chemistry should hydrolyze a single phosphoric anhydride substrate",
            )

        substrate = left_reactive[0]
        substrate_p_count = self._atom_count(substrate, 15)
        if substrate_p_count < 2 or not self._has_phosphoric_anhydride(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="Substrate lacks a phosphoric anhydride (P-O-P) bond",
            )
        if self._largest_anhydride_phosphorus_component_size(substrate) != substrate_p_count:
            return ClassificationResult(
                is_member=False,
                explanation="Substrate contains disconnected phosphorus domains beyond a single anhydride chain",
            )

        phosphorus_products = [
            participant
            for participant in right_reactive
            if self._atom_count(participant, 15) > 0
        ]
        phosphorus_product_slots = sum(max(participant.count, 1) for participant in phosphorus_products)
        if phosphorus_product_slots < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Hydrolysis does not yield multiple phosphorus-bearing products",
            )

        has_nonphosphorus_product = any(
            self._atom_count(participant, 15) == 0 for participant in right_reactive
        )
        if has_nonphosphorus_product:
            return ClassificationResult(
                is_member=False,
                explanation="Additional non-phosphorus products indicate chemistry beyond phosphoric anhydride hydrolysis",
            )

        if self._is_organophosphorus(substrate) and not any(
            self._is_organophosphorus(participant) for participant in phosphorus_products
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Organic phosphorus substrate loses phosphate to an inorganic leaving group instead of phosphoric anhydride hydrolysis",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Pyrophosphatase: hydrolytic cleavage of a phosphoric anhydride into multiple phosphorus-bearing products",
        )

    @staticmethod
    def _atom_count(participant: Participant, atomic_num: int) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_num)

    @classmethod
    def _is_organophosphorus(cls, participant: Participant) -> bool:
        return cls._atom_count(participant, 15) > 0 and cls._atom_count(participant, 6) > 0

    @classmethod
    def _has_phosphoric_anhydride(cls, participant: Participant) -> bool:
        if cls.PHOSPHORIC_ANHYDRIDE_PATTERN is None:
            return False
        mol = participant.get_mol()
        return mol is not None and mol.HasSubstructMatch(cls.PHOSPHORIC_ANHYDRIDE_PATTERN)

    @staticmethod
    def _largest_anhydride_phosphorus_component_size(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0

        phosphorus_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
        if not phosphorus_atoms:
            return 0

        adjacency: dict[int, set[int]] = {idx: set() for idx in phosphorus_atoms}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 8:
                continue
            phosphorus_neighbors = [
                neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 15
            ]
            if len(phosphorus_neighbors) < 2:
                continue
            for index, phosphorus_idx in enumerate(phosphorus_neighbors):
                for other_idx in phosphorus_neighbors[index + 1 :]:
                    adjacency[phosphorus_idx].add(other_idx)
                    adjacency[other_idx].add(phosphorus_idx)

        visited: set[int] = set()
        largest_component = 0
        for phosphorus_idx in phosphorus_atoms:
            if phosphorus_idx in visited:
                continue
            stack = [phosphorus_idx]
            component_size = 0
            while stack:
                current_idx = stack.pop()
                if current_idx in visited:
                    continue
                visited.add(current_idx)
                component_size += 1
                stack.extend(adjacency[current_idx] - visited)
            largest_component = max(largest_component, component_size)
        return largest_component
