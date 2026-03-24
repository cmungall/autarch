"""3-hydroxyacyl-CoA dehydratase activity.

Catalysis of the dehydration of a 3-hydroxyacyl-CoA thioester to the
corresponding enoyl-CoA.
"""

from collections import Counter

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.lyase import Lyase
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE


class ThreeHydroxyacylCoADehydratase(Lyase):
    """3-hydroxyacyl-CoA dehydratase activity.

    Catalysis of the dehydration of a 3-hydroxyacyl-CoA thioester to the
    corresponding enoyl-CoA.
    """

    GO_ID = "GO:0018812"
    EC_NUMBER_PREFIX = "4.2.1.-"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for reversible dehydration of a 3-hydroxyacyl-CoA thioester."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not 3-hydroxyacyl-CoA dehydratase chemistry",
            )

        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward

        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse hydration orientation)",
            )

        return forward

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        if not any(participant.chebi_id == CHEBI_H2O for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No water product found for the dehydration branch",
            )

        left_core = [
            participant
            for participant in left_participants
            if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS} and participant.smiles is not None
        ]
        right_core = [
            participant
            for participant in right_participants
            if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS} and participant.smiles is not None
        ]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected a single CoA thioester substrate/product pair",
            )

        substrate = left_core[0]
        product = right_core[0]
        if not self._is_hydroxyacyl_coa(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="No 3-hydroxyacyl-CoA substrate detected",
            )
        if not self._is_enoyl_coa(product):
            return ClassificationResult(
                is_member=False,
                explanation="No enoyl-CoA product detected",
            )
        if self._core_formula(substrate) != self._core_formula(product):
            return ClassificationResult(
                is_member=False,
                explanation="Substrate and product do not preserve the acyl-CoA scaffold",
            )

        return ClassificationResult(
            is_member=True,
            explanation="3-hydroxyacyl-CoA dehydratase: dehydration of a hydroxyacyl-CoA thioester to an enoyl-CoA",
        )

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    @classmethod
    def _is_coa_like(cls, participant: Participant) -> bool:
        return bool(participant.is_thioester() and cls._phosphorus_count(participant) >= 2)

    @classmethod
    def _is_hydroxyacyl_coa(cls, participant: Participant) -> bool:
        positions = cls._acyl_chain_positions(participant)
        if positions is None:
            return False
        _, alpha_carbon, beta_carbon = positions
        mol = participant.get_mol()
        if mol is None:
            return False
        alpha_beta_bond = mol.GetBondBetweenAtoms(alpha_carbon.GetIdx(), beta_carbon.GetIdx())
        if alpha_beta_bond is None or alpha_beta_bond.GetBondTypeAsDouble() != 1:
            return False
        return cls._has_beta_hydroxyl(beta_carbon) and cls._is_fatty_acyl_fragment(
            beta_carbon, excluded_atom_ids={alpha_carbon.GetIdx()}
        )

    @classmethod
    def _is_enoyl_coa(cls, participant: Participant) -> bool:
        positions = cls._acyl_chain_positions(participant)
        if positions is None:
            return False
        _, alpha_carbon, beta_carbon = positions
        mol = participant.get_mol()
        if mol is None:
            return False
        alpha_beta_bond = mol.GetBondBetweenAtoms(alpha_carbon.GetIdx(), beta_carbon.GetIdx())
        return alpha_beta_bond is not None and alpha_beta_bond.GetBondTypeAsDouble() == 2

    @classmethod
    def _acyl_chain_positions(cls, participant: Participant):
        mol = participant.get_mol()
        if mol is None or not cls._is_coa_like(participant):
            return None
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            has_thioester_sulfur = any(
                neighbor.GetAtomicNum() == 16
                and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 1
                for neighbor in atom.GetNeighbors()
            )
            has_carbonyl_oxygen = any(
                neighbor.GetAtomicNum() == 8
                and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2
                for neighbor in atom.GetNeighbors()
            )
            if not (has_thioester_sulfur and has_carbonyl_oxygen):
                continue
            alpha_candidates = [
                neighbor
                for neighbor in atom.GetNeighbors()
                if neighbor.GetAtomicNum() == 6
            ]
            if len(alpha_candidates) != 1:
                continue
            alpha_carbon = alpha_candidates[0]
            beta_candidates = [
                neighbor
                for neighbor in alpha_carbon.GetNeighbors()
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != atom.GetIdx()
            ]
            if len(beta_candidates) != 1:
                continue
            return atom, alpha_carbon, beta_candidates[0]
        return None

    @staticmethod
    def _has_beta_hydroxyl(beta_carbon) -> bool:
        for neighbor in beta_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() != 8:
                continue
            bond = beta_carbon.GetOwningMol().GetBondBetweenAtoms(beta_carbon.GetIdx(), neighbor.GetIdx())
            if bond is not None and bond.GetBondTypeAsDouble() == 1:
                return True
        return False

    @classmethod
    def _is_fatty_acyl_fragment(cls, start_atom, excluded_atom_ids: set[int]) -> bool:
        fragment_atoms = cls._collect_fragment_atoms(start_atom, excluded_atom_ids)
        fragment_ids = {atom.GetIdx() for atom in fragment_atoms}
        if not any(atom.GetAtomicNum() in {0, 6} and atom.GetIdx() != start_atom.GetIdx() for atom in fragment_atoms):
            return False

        oxygen_count = 0
        for atom in fragment_atoms:
            atomic_num = atom.GetAtomicNum()
            if atomic_num in {0, 6}:
                pass
            elif atomic_num == 8 and atom.GetIdx() != start_atom.GetIdx():
                oxygen_count += 1
            else:
                return False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in fragment_ids or neighbor.GetIdx() in excluded_atom_ids:
                    continue
                if neighbor.GetAtomicNum() != 1:
                    return False
        return oxygen_count <= 1

    @staticmethod
    def _collect_fragment_atoms(start_atom, excluded_atom_ids: set[int]):
        fragment_atoms = []
        stack = [start_atom]
        seen = set(excluded_atom_ids)
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in seen:
                continue
            seen.add(atom.GetIdx())
            fragment_atoms.append(atom)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in seen or neighbor.GetAtomicNum() not in {0, 6, 8}:
                    continue
                bond = atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond is None or bond.GetBondTypeAsDouble() != 1:
                    continue
                stack.append(neighbor)
        return fragment_atoms

    @staticmethod
    def _core_formula(participant: Participant) -> Counter[int]:
        mol = participant.get_mol()
        if mol is None:
            return Counter()
        return Counter(
            atom.GetAtomicNum()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() not in {0, 1, 8}
        )
